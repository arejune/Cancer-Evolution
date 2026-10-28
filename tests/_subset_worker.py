"""Runs a fixed subset of one experiment and prints its outputs as JSON. Used by tests/regression.py.

  python tests/_subset_worker.py EXPERIMENT

For agent-based runs it also records the simulator's random-number state at the end of each run, by
wrapping the simulator class so that each instance registers itself (this draws no random numbers).
"""
import json, os, sys, hashlib
from pathlib import Path

EXP = sys.argv[1]
REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "src")]
T = "threshold_transfer."
MOD = dict(sim=T + "model.simulator", main=T + "experiments.transfer_main", arms=T + "experiments.transfer_arms",
           ctmc=T + "model.ctmc", rmap=T + "experiments.recal_map", rsweep=T + "experiments.recal_sweep",
           log=T + "experiments.statesuff_log", policy=T + "experiments.policy", grid=T + "threshold.grid",
           pilot=T + "threshold.construction")
import importlib
import numpy as np

_LAST = []


def _patch_simulators():
    base = importlib.import_module(MOD["sim"]).TumorSimulator
    if getattr(base, "_tt_wrapped", False): return

    class Wrapped(base):
        _tt_wrapped = True

        def __init__(self, *a, **k):
            super().__init__(*a, **k); _LAST.append(self)
    for key in ("sim", "pilot", "rmap"):     # every module that builds a simulator
        importlib.import_module(MOD[key]).TumorSimulator = Wrapped


def _state():
    return json.dumps(_LAST[-1]._rng.bit_generator.state, sort_keys=True, default=str)


def _init():
    _patch_simulators()


def jsonable(x):
    if isinstance(x, dict): return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray): return jsonable(x.tolist())
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.floating, float)): return repr(float(x))
    if isinstance(x, np.bool_): return bool(x)
    return x


# per-run jobs, at module level so that worker processes can import them
def j_main(a):
    r = importlib.import_module(MOD["main"]).job(a); return [r, _state()]


def j_transit(a):
    r = importlib.import_module(MOD["rmap"]).transit(a); return [r, _state()]


def j_sweep(a):
    r = importlib.import_module(MOD["rsweep"]).job(a); return [r, _state()]


def j_log(a):
    return importlib.import_module(MOD["log"]).job_logged(a)


def j_policy(a):
    PL = importlib.import_module(MOD["policy"])
    if not PL._M: PL._M.update(PL.load_models(check=True))
    seed, arm, tau, setting, thr_size = a
    r = PL.run_arm(seed, arm, tau, setting, thr_size)
    r["eps"] = [{k: v for k, v in e.items() if k != "feats"} for e in r["eps"]]
    return r


def pool_map(f, args):
    from multiprocessing import Pool
    with Pool(int(os.environ.get("TT_WORKERS", "7")), initializer=_init) as pool:
        return pool.map(f, args, chunksize=4)


def ctmc_cells(cells):
    """CTMC batches at full R = 500; records the generator state after each batch."""
    C = importlib.import_module(MOD["ctmc"])
    gens = []
    orig = np.random.default_rng
    def rec(*a, **k):
        g = orig(*a, **k); gens.append(g); return g
    np.random.default_rng = rec
    out = []
    try:
        for kw in cells:
            r = C.batch(**kw)
            out.append(dict(prem=r["prem"], t_prog=np.where(np.isnan(r["t_prog"]), -1.0, r["t_prog"]),
                            burden=r["burden"], counts=r["counts"],
                            rng=json.dumps(gens[-1].bit_generator.state, sort_keys=True, default=str)))
    finally:
        np.random.default_rng = orig
    return out


def experiment(name):
    G = importlib.import_module(MOD["grid"])
    if name == "main_abm":            # 2 cells x 100 seeds; the m = 1.00 cell uses the counted simulator
        args = [("AT", 7, 1.00, s) for s in range(100)] + [("AT", 28, 0.95, s) for s in range(100)]
        return pool_map(j_main, args)
    if name == "main_arms":           # ODE cell + CTMC-nominal and CTMC-pk at two cells, R = 500 each
        A = importlib.import_module(MOD["arms"])
        out = {}
        for t, m in ((7, 1.00), (56, 1.20)):
            thr = m * A.NS[t]
            out[f"ode|{t}|{m:.2f}"] = A.ode_cell(t, thr)
            cells = [dict(p=A.p, R=G.R, block_seed=G.BLOCK_ABM, policy="AT", tau=t, thr=thr,
                          horizon_days=float(G.T_END), pk=pk) for pk in (None, A.PK)]
            r = ctmc_cells(cells)
            out[f"ctmc_nom|{t}|{m:.2f}"], out[f"ctmc_pk|{t}|{m:.2f}"] = r
        return out
    if name == "recal_map":           # 2 starting sizes x 100 seeds
        args = [(n, s) for n in (2600, 3000) for s in range(100)]
        return pool_map(j_transit, args)
    if name == "recal_sweep":         # 2 cells x 100 seeds, thresholds from the archived N*_recal
        res = importlib.import_module(T + "paths").RESULTS
        z = np.load(Path(res) / "p_t_recal_sweep.npz"); NSR = dict(zip(z["taus"].tolist(), z["nstar_recal"]))
        args = [(t, m, m * NSR[t], s) for t, m in ((14, 0.85), (56, 0.90)) for s in range(100)]
        return pool_map(j_sweep, args)
    if name == "recal_ctmc":          # CTMC under the re-fitted law, 2 cells, R = 500; and the fine grid
        C = importlib.import_module(MOD["ctmc"])
        res = importlib.import_module(T + "paths").RESULTS
        z = np.load(Path(res) / "p_t_recal_sweep.npz"); NSR = dict(zip(z["taus"].tolist(), z["nstar_recal"]))
        law = (float(z["c_refit"]), float(z["phi_refit"]))
        p = importlib.import_module(T + "model.ode").ODE2Params(mu_a=0.0, mu_escape=0.0)
        cells = [dict(p=p, R=G.R, block_seed=1_900_000, policy="AT", tau=28, thr=0.80 * NSR[28], nmax=C.N_MAX, horizon_days=750.0, growth=law),
                 dict(p=p, R=G.R, block_seed=2_000_000, policy="AT", tau=7, thr=0.9750 * NSR[7], nmax=C.N_MAX, horizon_days=750.0, growth=law),
                 dict(p=p, R=G.R, block_seed=2_000_000, policy="AT", tau=56, thr=1.0250 * NSR[56], nmax=C.N_MAX, horizon_days=750.0, growth=law)]
        return ctmc_cells(cells)
    if name == "statesuff_log":       # one training and one test cell, 100 seeds each
        args = [("train", 28, 2810, s) for s in range(100)] + [("test", 7, 3180, s) for s in range(100)]
        return pool_map(j_log, args)
    if name == "policy":              # scoring seeds: history-and-image rule at tau = 7 and size-only at tau = 28, 100 each
        PL = importlib.import_module(MOD["policy"])
        sel = PL.selected()
        args = [(2_500_000 + s, "CLIN", 7, sel["CLIN|7"]["selected"], sel["SIZE|7"]["selected"]) for s in range(100)] + \
               [(2_500_000 + s, "SIZE", 28, sel["SIZE|28"]["selected"], sel["SIZE|28"]["selected"]) for s in range(100)]
        return pool_map(j_policy, args)
    raise KeyError(name)


if __name__ == "__main__":
    _patch_simulators()
    out = jsonable(experiment(EXP))
    blob = json.dumps(out, sort_keys=True)
    print(json.dumps(dict(experiment=EXP, sha256=hashlib.sha256(blob.encode()).hexdigest(), out=out)))
