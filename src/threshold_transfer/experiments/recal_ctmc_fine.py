"""CTMC-nominal under the recalibrated law on a finer grid of multipliers around each crossing of q.

The coarse grid (step 0.05 in m) jumped from p_prem near 0 to near 1 within one step at several intervals, so
each crossing is remeasured in steps of 0.0025 over WIN[tau], 500 replicates from seed 2,000,000.
Round 1 covers tau = 7, 42, 84 and round 2 tau = 14, 28, 56, as in the original runs.
"""
import json, time
import numpy as np
from threshold_transfer.paths import RESULTS
from threshold_transfer.model import ctmc as C
from threshold_transfer.model.ode import N_MAX, P as p
from threshold_transfer.threshold import grid as G
from threshold_transfer.experiments.recal_ctmc import recal_law

BLOCK = 2_000_000
STEP = 0.0025
WIN = {7: (0.9500, 1.0000), 42: (1.0000, 1.0500), 84: (1.0500, 1.1000),
       14: (0.9500, 1.0000), 28: (0.9500, 1.0000), 56: (1.0000, 1.0500)}
ROUNDS = {1: ((7, 42, 84), "p_t_recal_ctmc_fine.json"), 2: ((14, 28, 56), "p_t_recal_ctmc_fine2.json")}


def main(round_=2):
    taus_run, fname = ROUNDS[round_]
    TAUS, NSR, c, pb = recal_law()
    out = {}; t0 = time.time()
    for t in taus_run:
        lo, hi = WIN[t]
        for m in np.round(np.arange(lo, hi + 1e-9, STEP), 4):
            thr = float(m) * NSR[t]
            r = C.batch(p, G.R, BLOCK, policy="AT", tau=t, thr=thr, nmax=N_MAX,
                        horizon_days=float(G.T_END), growth=(c, pb))
            out[f"{t}|{m:.4f}"] = dict(prem=r['prem'].tolist(), thr=thr)
            print(f"tau={t} m={m:.4f}  p={r['prem'].mean():.4f}", flush=True)
    json.dump(dict(block=BLOCK, R=G.R, nmax=N_MAX, step=STEP, c=c, phi=pb,
                   nstar_recal={str(t): float(NSR[t]) for t in TAUS},
                   wall_s=time.time() - t0, cells=out),
              open(RESULTS / fname, "w"))
    print(f"wrote {fname}")
