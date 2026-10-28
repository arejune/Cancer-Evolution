"""Main run in the agent-based model: the size-only policy at theta = m * N*(tau) on the grid, plus continuous
treatment (CT), 500 seeds per cell from seed 1,200,000.

Each file data/results/p_t_cells/<cell>__<first seed>.npz holds CHUNK seeds and is skipped if it exists, so an
interrupted run resumes. When all files exist they are joined into p_t_main_raw.npz. Division counters are
stored for the cells at m = 1.00 and CT only (zeros elsewhere).
"""
import os
for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
import sys, time, json, contextlib
import numpy as np
from multiprocessing import Pool
from threshold_transfer.paths import RESULTS, WORKERS
from threshold_transfer.threshold import construction as P, grid as G
from threshold_transfer.model.ode import nstar

CELLDIR = RESULTS / "p_t_cells"
_NSF = RESULTS / "p_t_nstar.json"
if _NSF.exists():
    NS = {int(k): float(v) for k, v in json.load(open(_NSF)).items()}
else:
    NS = {t: nstar(float(t)) for t in G.TAUS}
    RESULTS.mkdir(parents=True, exist_ok=True)
    json.dump({str(k): v for k, v in NS.items()}, open(_NSF, "w"))
GRID = G.grid(NS)
DT = [('pol', 'U3'), ('tau', 'i4'), ('m', 'f8'), ('seed', 'i4'), ('prem', 'i4'),
      ('t_prog', 'f8'), ('burden', 'f8'), ('nmax', 'f8')]


def job(a):
    """One run, a = (policy, tau, m, seed index). Undefined tau, m and t_prog are stored as -1."""
    pol, tau, m, s = a
    keep_led = (pol == "CT") or (abs(m - 1.00) < 1e-9)
    kw = dict(tau=tau, thr=m * NS[tau]) if pol == "AT" else {}
    r = P.run(G.BLOCK_ABM + s, pol, N0=G.N0, horizon=G.T_END, **kw)
    led = r["led"].ravel().tolist() if keep_led else [0.0] * 12
    return (pol, -1 if tau is None else tau, -1.0 if m is None else m, s,
            r["prem"], -1.0 if np.isnan(r["t_prog"]) else r["t_prog"],
            r["burden"], r["N_max"], led)


def cells_in_order():
    """All cells, m >= 0.90 and CT first."""
    t1 = [("AT", t, m) for t in G.TAUS for m in GRID[t] if m >= 0.90]
    t1.append(("CT", None, None))
    t2 = [("AT", t, m) for t in G.TAUS for m in GRID[t] if m < 0.90]
    return t1 + t2


def cell_file(pol, tau, m, c):
    return CELLDIR / f"{'CT' if pol == 'CT' else f'AT_{tau}_{m:.2f}'}__{c:03d}.npz"


def main(CHUNK=25, NW=WORKERS):
    CELLDIR.mkdir(parents=True, exist_ok=True)
    cells = cells_in_order()
    units = [(pol, tau, m, c) for (pol, tau, m) in cells for c in range(0, G.R, CHUNK)]
    todo = [u for u in units if not cell_file(*u).exists()]
    print(f"{len(todo)} of {len(units)} files to run ({CHUNK} seeds each)", flush=True)
    t_all = time.time()
    pool = Pool(NW) if NW > 0 else None
    mapper = pool.map if pool else (lambda f, xs, chunksize=1: [f(x) for x in xs])
    with (pool if pool else contextlib.nullcontext()):
        for n, (pol, tau, m, c) in enumerate(todo, 1):
            f = cell_file(pol, tau, m, c)
            if f.exists():
                continue
            t0 = time.time()
            out = mapper(job, [(pol, tau, m, s) for s in range(c, c + CHUNK)], chunksize=4)
            rec = np.array([o[:8] for o in out], dtype=DT)
            led = np.array([o[8] for o in out], dtype=float)
            np.savez(f, rec=rec, led=led)
            print(f"[{n}/{len(todo)}] {f.name}  {time.time() - t0:.1f} s", flush=True)

    if any(not cell_file(*u).exists() for u in units):
        print("incomplete", flush=True); sys.exit(2)
    recs, leds = [], []
    for pol, tau, m in cells:
        for c in range(0, G.R, CHUNK):
            z = np.load(cell_file(pol, tau, m, c))
            recs.append(z['rec']); leds.append(z['led'])
    rec = np.concatenate(recs); led = np.concatenate(leds)
    np.savez(RESULTS / "p_t_main_raw.npz", rec=rec, led=led,
             nstar=np.array([NS[t] for t in G.TAUS]), taus=np.array(G.TAUS))
    days = float(np.where(rec['t_prog'] > 0, rec['t_prog'], G.T_END).sum())
    wall = time.time() - t_all
    print(f"wrote p_t_main_raw.npz: {len(rec)} runs, {wall:.1f} s", flush=True)
    json.dump(dict(wall_s=wall, units_run=len(todo), n=len(rec), workers=NW, sim_days=days),
              open(RESULTS / "p_t_main_cost.json", "w"))
