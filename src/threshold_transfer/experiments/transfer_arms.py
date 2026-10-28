"""ODE, CTMC-nominal and CTMC-pk arms of the main run, on the same grid of thresholds.

Each cell runs once in the deterministic ODE and with 500 replicates (seed 1,200,000) in each CTMC.
Output: data/results/p_t_main_arms.json.
"""
import json
import numpy as np
from threshold_transfer.paths import RESULTS
from threshold_transfer.threshold import grid as G
from threshold_transfer.model import ctmc as C
from threshold_transfer.model.ode import P as p, nstar, ode_cell

PK = (0.10004, 0.13420)   # (kappa_pk, k_el) per day, fitted to the agent-based model's mean drug level
NS = {t: nstar(float(t)) for t in G.TAUS}
GRID = G.grid(NS)


def _ctmc(policy, pk, **kw):
    r = C.batch(p, G.R, G.BLOCK_ABM, policy=policy, horizon_days=float(G.T_END), pk=pk, **kw)
    return dict(prem=r['prem'].tolist(), t_prog=np.where(np.isnan(r['t_prog']), -1.0, r['t_prog']).tolist())


def main():
    out = {}
    for t in G.TAUS:
        for m in GRID[t]:
            thr = m * NS[t]
            o = ode_cell(t, thr)
            for arm, pk in (("ctmc_nom", None), ("ctmc_pk", PK)):
                out[f"{arm}|{t}|{m:.2f}"] = _ctmc("AT", pk, tau=t, thr=thr)
            out[f"ode|{t}|{m:.2f}"] = o
            print(f"tau={t} m={m:.2f}  ODE prem={o['prem']}  "
                  f"CTMC-nominal p={np.mean(out[f'ctmc_nom|{t}|{m:.2f}']['prem']):.4f}  "
                  f"CTMC-pk p={np.mean(out[f'ctmc_pk|{t}|{m:.2f}']['prem']):.4f}", flush=True)
    for arm, pk in (("ctmc_nom", None), ("ctmc_pk", PK)):
        out[f"{arm}|CT|CT"] = _ctmc("CT", pk)
    json.dump({"nstar": {str(k): v for k, v in NS.items()},
               "grid": {str(k): list(v) for k, v in GRID.items()},
               "pk": list(PK), "cells": out},
              open(RESULTS / "p_t_main_arms.json", "w"))
    print("wrote p_t_main_arms.json")
