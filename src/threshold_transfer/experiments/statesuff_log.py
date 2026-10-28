"""Prediction study: the size-only policy with the tumour's state recorded at every decision.

Seed blocks: gate 2,300,000 (10 seeds), training set 2,100,000 and test set 2,200,000 (500 seeds each), at five
thresholds per interval (BAND). At each decision the run records cell counts, the mean drug level over
occupied sites C and the shape of the tumour (rim cells B, resistant rim cells B_R, connected components K,
radius of gyration Rg, maximum radius Rmax, distance to the lattice edge). Reading the state draws no random
numbers, so every run equals construction.run with the same seed; the gate checks this.
Output: data/results/p_s_cells/{train,test}/<tau>_<thr>__<first seed>.npz (25 seeds per file), gate.json.
"""
import os
for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
import sys, time, json
import numpy as np
from scipy import ndimage
from multiprocessing import Pool
from threshold_transfer.paths import RESULTS, WORKERS
from threshold_transfer.threshold import construction as P

N0 = 3001.0
T_END = 750
PROG = P.PROG_MUL * N0
TAUS = (7, 14, 28, 56)
BAND = {7: (3060, 3120, 3180, 3240, 3300),
        14: (2900, 2950, 3000, 3050, 3100),
        28: (2670, 2740, 2810, 2880, 2950),
        56: (2660, 2740, 2820, 2900, 2980)}
CELLS = [(t, thr) for t in TAUS for thr in BAND[t]]
BLOCKS = {"gate": 2_300_000, "train": 2_100_000, "test": 2_200_000}
NSEED = {"gate": 10, "train": 500, "test": 500}
SHARD = 25
NW = WORKERS
OUT = RESULTS / "p_s_cells"
SENS, RA, ESC = 1, 2, 3
_S3 = np.ones((3, 3), dtype=int)

EP_DT = [('tau', 'i2'), ('thr', 'i2'), ('s', 'i2'), ('e', 'i2'), ('day', 'i2'),
         ('N', 'i4'), ('S', 'i4'), ('R', 'i4'), ('Esc', 'i4'), ('on', 'i1'),
         ('C', 'f8'), ('B', 'i4'), ('B_R', 'i4'), ('K', 'i2'), ('Rg', 'f8'),
         ('Rmax', 'f8'), ('d_wall', 'i2'), ('y', 'i1'), ('trunc', 'i1')]
RUN_DT = [('tau', 'i2'), ('thr', 'i2'), ('s', 'i2'), ('prem', 'i4'), ('t_prog', 'f8'),
          ('burden', 'f8'), ('N_end', 'f8'), ('N_max', 'f8')]


def epoch_reads(sim):
    """State at a decision: counts, mean drug level over occupied sites and shape features."""
    st = sim.get_state()
    g = np.array(sim.grid, copy=True)
    occ = g != 0
    Nx, Ny = g.shape
    ii, jj = np.nonzero(occ)
    C = float(np.asarray(sim.C_A)[occ].mean())
    r0, r1, c0, c1 = ii.min(), ii.max(), jj.min(), jj.max()
    d_wall = int(min(r0, c0, Nx - 1 - r1, Ny - 1 - c1))
    # crop with a 1-site margin; sites outside the lattice count as empty
    sub = np.pad(occ[max(r0-1, 0):r1+2, max(c0-1, 0):c1+2], 1, constant_values=False)
    gs = np.pad(g[max(r0-1, 0):r1+2, max(c0-1, 0):c1+2], 1, constant_values=0)
    h, w = sub.shape
    nb = np.zeros((h - 2, w - 2), dtype=np.int16)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                nb += sub[1+di:h-1+di, 1+dj:w-1+dj]
    core = sub[1:-1, 1:-1]
    bnd = core & (nb < 8)
    B = int(bnd.sum())
    B_R = int((bnd & (gs[1:-1, 1:-1] == RA)).sum())
    K = int(ndimage.label(occ, structure=_S3)[1])
    ci, cj = ii.mean(), jj.mean()
    d2 = (ii - ci) ** 2 + (jj - cj) ** 2
    return dict(N=int(st["true_total_volume"]), S=int(st["n_S"]), R=int(st["n_Ra"]),
                Esc=int(st["n_Escape"]), C=C, B=B, B_R=B_R, K=K,
                Rg=float(np.sqrt(d2.mean())), Rmax=float(np.sqrt(d2.max())), d_wall=d_wall)


def run_logged(seed, tau, thr, log=True, horizon=T_END):
    """construction.run(policy="AT") plus the state at each decision.

    y = 1 on the last decision if progression happens in the hold after it; trunc = 1 if the last interval is
    cut by the horizon. rng0 and rng1 are the simulator's generator state before and after the run.
    """
    sim, _, _ = P.build(seed)
    rng0 = json.dumps(sim._rng.bit_generator.state, sort_keys=True, default=str)
    prog = P.PROG_MUL * N0
    burden = 0.0; on = True; prem = 0; t_prog = np.nan
    Ntrace_max = 0.0; below_at_epoch = False
    eps = []
    for k in range(1, horizon + 1):
        if (k - 1) % tau == 0:
            V = float(sim.get_state()["true_total_volume"])
            on = (V >= thr)
            below_at_epoch = (V < thr)
            if log:
                r = epoch_reads(sim)
                assert r["N"] == int(V)
                r.update(e=(k - 1) // tau, day=k - 1, on=int(on), y=0, trunc=0)
                eps.append(r)
        sim.step(u_A=P.DOSE if on else 0.0, u_B=0.0)
        st = sim.get_state(); N = float(st["true_total_volume"])
        burden += N
        if N > Ntrace_max: Ntrace_max = N
        if N >= prog:
            if np.isnan(t_prog): t_prog = float(k)
            if below_at_epoch: prem = 1
            break
        if N <= 0: break
    if log and eps:
        last = eps[-1]
        if prem == 1:
            last["y"] = 1                      # progression inside this hold
        elif np.isnan(t_prog) and last["day"] + tau > horizon and N > 0:
            last["trunc"] = 1                  # interval cut by the horizon; label undefined
    rng1 = json.dumps(sim._rng.bit_generator.state, sort_keys=True, default=str)
    return dict(prem=prem, t_prog=t_prog, burden=burden, N_end=N, N_max=Ntrace_max,
                eps=eps, rng0=rng0, rng1=rng1)


def _warm(*_):
    """Pool initialiser: compile the simulator kernels."""
    from threshold_transfer.model.simulator import warmup_jit
    warmup_jit()


def job_logged(a):
    """One logged run, a = (role, tau, thr, s). Returns (run row, decision rows, rng0, rng1)."""
    role, tau, thr, s = a
    r = run_logged(BLOCKS[role] + s, tau, thr)
    run = (tau, thr, s, r["prem"], -1.0 if np.isnan(r["t_prog"]) else r["t_prog"],
           r["burden"], r["N_end"], r["N_max"])
    eps = [(tau, thr, s, x["e"], x["day"], x["N"], x["S"], x["R"], x["Esc"], x["on"],
            x["C"], x["B"], x["B_R"], x["K"], x["Rg"], x["Rmax"], x["d_wall"],
            x["y"], x["trunc"]) for x in r["eps"]]
    return run, eps, r["rng0"], r["rng1"]


def job_ref(a):
    """One unlogged run through construction.run, for the gate's equality check."""
    role, tau, thr, s = a
    r = P.run(BLOCKS[role] + s, "AT", tau=tau, thr=float(thr), N0=N0, horizon=T_END)
    return (r["prem"], -1.0 if np.isnan(r["t_prog"]) else r["t_prog"],
            r["burden"], float(r["N_end"]), r["N_max"])


def job_replica(a):
    """One run with logging off, for the gate's RNG-state check."""
    role, tau, thr, s = a
    r = run_logged(BLOCKS[role] + s, tau, thr, log=False)
    return r["rng0"], r["rng1"]


def save_shard(path, outs):
    """Write one shard file atomically (temporary file, then rename)."""
    runs = np.array([o[0] for o in outs], dtype=RUN_DT)
    eps = np.array([e for o in outs for e in o[1]], dtype=EP_DT)
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, runs=runs, eps=eps)
    os.replace(tmp, path)


def gate():
    """Check the logging code on the gate block, 10 seeds per cell: (1) runs equal construction.run; (2) the
    generator state matches an unlogged rerun; (3) labels match the outcomes; (4) no escape cells and the tumour
    stays at least 3 sites from the edge. Also times the runs. Writes gate.npz and gate.json.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    args = [("gate", t, thr, s) for (t, thr) in CELLS for s in range(NSEED["gate"])]
    with Pool(NW, initializer=_warm) as pool:
        pool.map(_warm, range(NW))                         # compile before timing
        t0 = time.time(); ref = pool.map(job_ref, args, chunksize=2); w_ref = time.time() - t0
        t0 = time.time(); lg = pool.map(job_logged, args, chunksize=2); w_log = time.time() - t0
        rep = pool.map(job_replica, args, chunksize=2)
    runs = np.array([o[0] for o in lg], dtype=RUN_DT)
    eps = np.array([e for o in lg for e in o[1]], dtype=EP_DT)
    np.savez(OUT / "gate.npz", runs=runs, eps=eps)
    # item 1: exact equality with the unlogged reference
    fields = ("prem", "t_prog", "burden", "N_end", "N_max")
    mism = [(a, f) for a, rf, rr in zip(args, ref, runs)
            for f, v in zip(fields, rf) if float(rr[f]) != float(v)]
    # item 2: generator state, logged against unlogged replica; and the state must have moved
    rng_eq = sum(o[3] == r[1] for o, r in zip(lg, rep))
    rng_moved = sum(o[3] != o[2] for o in lg)
    # item 3: labels
    lab_bad = 0
    for rr in runs:
        m = (eps["tau"] == rr["tau"]) & (eps["thr"] == rr["thr"]) & (eps["s"] == rr["s"])
        ysum = int(eps["y"][m & (eps["on"] == 0)].sum())
        lab_bad += int(ysum != rr["prem"]) + int(eps["y"][m & (eps["on"] == 1)].sum() != 0)
    # item 4: no escape cells; distance to the lattice edge
    esc_bad = int((eps["Esc"] != 0).sum()); wall_min = int(eps["d_wall"].min())
    ratio = w_log / w_ref
    proj = 72.6 * ratio                       # projected minutes for the full run (72.6 min uninstrumented)
    res = dict(n=len(args), item1_mismatches=len(mism), item1_examples=[str(x) for x in mism[:5]],
               item2_rng_equal=rng_eq, item2_rng_advanced=rng_moved,
               item3_label_violations=lab_bad, item4_escape_nonzero=esc_bad,
               item4_min_d_wall=wall_min, wall_ref_s=w_ref, wall_log_s=w_log,
               ms_per_run_ref=1000 * w_ref / len(args), ms_per_run_log=1000 * w_log / len(args),
               overhead_ratio=ratio, projected_wall_min=proj,
               n_epochs=len(eps), n_trunc=int(eps["trunc"].sum()), n_events=int(runs["prem"].sum()))
    item2 = ("PASS" if rng_eq == len(args) else "FAIL") if rng_moved == len(args) else "UNINFORMATIVE"
    res["item1"] = "PASS" if not mism else "FAIL"
    res["item2"] = item2
    res["item3"] = "PASS" if lab_bad == 0 else "FAIL"
    res["item4"] = "PASS" if (esc_bad == 0 and wall_min >= 3) else "FAIL"
    res["GATE"] = "PASS" if (res["item1"] == "PASS" and item2 != "FAIL"
                             and res["item3"] == "PASS" and res["item4"] == "PASS") else "FAIL"
    res["LAUNCH"] = bool(res["GATE"] == "PASS" and proj <= 150.0)
    json.dump(res, open(OUT / "gate.json", "w"), indent=1)
    print(json.dumps(res, indent=1), flush=True)


def units():
    """All output files: (role, tau, thr, first seed), training set first."""
    return [(role, t, thr, c) for role in ("train", "test") for (t, thr) in CELLS
            for c in range(0, NSEED[role], SHARD)]


def shard_path(role, t, thr, c):
    """Path of one shard file."""
    return OUT / role / f"{t}_{thr}__{c:03d}.npz"


def supervise(budget_s):
    """Run missing files until all exist or budget_s seconds pass. Needs a passed gate (gate.json)."""
    g = json.load(open(OUT / "gate.json"))
    if not g["LAUNCH"]:
        print("gate.json says LAUNCH = False; refusing to run", flush=True); sys.exit(3)
    for role in ("train", "test"): (OUT / role).mkdir(parents=True, exist_ok=True)
    U = units(); todo = [u for u in U if not shard_path(*u).exists()]
    print(f"units {len(U)}  todo {len(todo)}  budget {budget_s:.0f}s", flush=True)
    t_start = time.time(); done = 0
    if todo:
        with Pool(NW, initializer=_warm) as pool:
            for (role, t, thr, c) in todo:
                if time.time() - t_start > budget_s: break
                outs = pool.map(job_logged, [(role, t, thr, s) for s in range(c, c + SHARD)],
                                chunksize=1)
                save_shard(shard_path(role, t, thr, c), outs); done += 1
    wall = time.time() - t_start
    left = sum(not shard_path(*u).exists() for u in U)
    with open(OUT / "wall_log.jsonl", "a") as f:
        f.write(json.dumps(dict(t=time.strftime("%Y-%m-%d %H:%M:%S"), wall_s=wall,
                                shards=done, left=left)) + "\n")
    print(f"invocation wall {wall:.1f}s  shards {done}  remaining {left}", flush=True)
    print("COMPLETE" if left == 0 else "INCOMPLETE", flush=True)
