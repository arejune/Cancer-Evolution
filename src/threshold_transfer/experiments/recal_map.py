"""Recalibration, step 1: drug-free transit times of the agent-based model.

For each starting size 2200, 2300, ..., 3600 cells, 200 drug-free runs (seed 1,800,000) record the day the
tumour reaches 1.2 N0 and its daily size. Each tumour starts as a compact disk trimmed to the starting size,
with 1% resistant cells; at 3001 cells this is the main run's initial tumour.
Output: data/results/p_t_recal_map.npz.
"""
import os
for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
import json, time
import numpy as np
from multiprocessing import Pool
from threshold_transfer.paths import RESULTS, WORKERS
from threshold_transfer.threshold import construction as P, grid as G
from threshold_transfer.model.ode import L, N0, PROG_MUL
from threshold_transfer.model.simulator import TumorSimulator, SimulatorParams

BLOCK = 1_800_000
PROG = PROG_MUL * N0
R_MAP = 200
MAP_N = list(range(2200, 3601, 100))
M_GRID = (0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15)   # multipliers of the recalibrated sweep
MAP_HOR = 400
DISK = {28: 2453, 29: 2629, 30: 2821, 31: 3001, 32: 3209, 33: 3409, 34: 3625}   # sites in a disk of each radius


def build_at(seed, n_target):
    """Compact disk of exactly n_target cells: the smallest disk that holds them, outermost sites removed
    (ties broken at random), then 1% made resistant. Uses the placement stream seed + 313,000,000.
    """
    sim = TumorSimulator(SimulatorParams(grid_size=(L, L), mu_a=0.0, mu_escape=0.0), seed=int(seed))
    radius = min(r for r, n in sorted(DISK.items()) if n >= n_target)
    sim.initialize_tumor(cluster_radius=radius, n_Ra_init=0, n_escape_init=0)
    g = sim.grid
    idx = np.argwhere(g == P.SENSITIVE)
    d = np.hypot(idx[:, 0] - L // 2, idx[:, 1] - L // 2)
    prng = np.random.default_rng(int(seed) + 313_000_000)
    order = np.lexsort((prng.random(len(d)), -d))     # outermost first
    n_drop = len(idx) - n_target
    for k in order[:n_drop]:
        r, c = idx[k]; g[r, c] = P.EMPTY
    keep = idx[order[n_drop:]]
    pick = prng.choice(len(keep), size=int(round(P.RFRAC * len(keep))), replace=False)
    for k in pick:
        r, c = keep[k]; g[r, c] = P.RA
    return sim


def transit(a):
    """Drug-free run, a = (n_start, seed index). Returns (n_start, s, days to 1.2 N0 or nan, daily sizes)."""
    n_start, s = a
    sim = build_at(BLOCK + s, n_start)
    tr = [float(sim.get_state()["true_total_volume"])]
    for k in range(1, MAP_HOR + 1):
        sim.step(u_A=0.0, u_B=0.0)
        N = float(sim.get_state()["true_total_volume"])
        tr.append(N)
        if N >= PROG:
            return n_start, s, float(k), np.array(tr, dtype=np.float32)
    return n_start, s, np.nan, np.array(tr, dtype=np.float32)


def main(NW=WORKERS):
    t0 = time.time()
    jobs = [(n, s) for n in MAP_N for s in range(R_MAP)]
    with Pool(NW) as pool:
        out = pool.map(transit, jobs, chunksize=8)
    wall = time.time() - t0
    T = {n: np.array([r[2] for r in out if r[0] == n]) for n in MAP_N}
    width = max(len(r[3]) for r in out)
    traces = {n: np.full((R_MAP, width), np.nan, np.float32) for n in MAP_N}
    row = {n: 0 for n in MAP_N}
    for n, s, t, tr in out:
        traces[n][row[n], :len(tr)] = tr; row[n] += 1
    np.savez_compressed(RESULTS / "p_t_recal_map.npz", map_n=np.array(MAP_N), R=R_MAP, block=BLOCK,
                        transit=np.array([T[n] for n in MAP_N]), traces=np.stack([traces[n] for n in MAP_N]))
    json.dump(dict(wall_s=wall, n=len(jobs), workers=NW), open(RESULTS / "p_t_recal_map_cost.json", "w"))
    print("N_start  median days  mean  sd")
    for n in MAP_N:
        f = T[n][np.isfinite(T[n])]
        print(f"{n:7d}  {np.median(f):11.2f}  {f.mean():5.2f}  {f.std():4.2f}")
    print(f"wrote p_t_recal_map.npz ({wall:.1f} s)")
