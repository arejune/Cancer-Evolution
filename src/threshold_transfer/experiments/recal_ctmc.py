"""CTMC-nominal and the ODE under the recalibrated growth law, on the cells of the recalibration sweep.

500 CTMC replicates per cell from seed 1,900,000. Output: data/results/p_t_recal_ctmc.json.
"""
import json, time
import numpy as np
from threshold_transfer.paths import RESULTS
from threshold_transfer.model import ctmc as C
from threshold_transfer.model.ode import N_MAX, P as p, ode_cell
from threshold_transfer.threshold import grid as G
from threshold_transfer.experiments import recal_map as RC

BLOCK = 1_900_000


def recal_law():
    """(taus, N*_recal by tau, refitted c, refitted phi_bulk) from p_t_recal_sweep.npz."""
    z = np.load(RESULTS / "p_t_recal_sweep.npz")
    taus = z['taus'].tolist()
    return taus, dict(zip(taus, z['nstar_recal'])), float(z['c_refit']), float(z['phi_refit'])


def main():
    TAUS, NSR, c, pb = recal_law()
    GRIDR = {t: tuple(m for m in RC.M_GRID if m < RC.PROG / NSR[t]) for t in TAUS}
    out = {}; t0 = time.time()
    for t in TAUS:
        for m in GRIDR[t]:
            thr = m * NSR[t]
            r = C.batch(p, G.R, BLOCK, policy="AT", tau=t, thr=thr, nmax=N_MAX,
                        horizon_days=float(G.T_END), growth=(c, pb))
            o = ode_cell(t, thr, c, pb)
            out[f"{t}|{m:.2f}"] = dict(
                ctmc_prem=r['prem'].tolist(),
                ctmc_tprog=np.where(np.isnan(r['t_prog']), -1.0, r['t_prog']).tolist(),
                ode_prem=int(o['prem']), thr=float(thr))
            print(f"tau={t} m={m:.2f}  ODE prem={o['prem']}  CTMC p={r['prem'].mean():.4f}", flush=True)
    json.dump(dict(block=BLOCK, R=G.R, nmax=N_MAX, c=c, phi=pb,
                   nstar_recal={str(t): float(NSR[t]) for t in TAUS},
                   wall_s=time.time() - t0, cells=out),
              open(RESULTS / "p_t_recal_ctmc.json", "w"))
    print("wrote p_t_recal_ctmc.json")
