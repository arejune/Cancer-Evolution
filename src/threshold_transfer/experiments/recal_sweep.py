"""Recalibration, steps 2-4.

2. Recalibrate the threshold to the agent-based model's drug-free regrowth (p_t_recal_map.npz) in two ways:
   (i) refit the growth constants (c, phi_bulk) to the median trajectories, (ii) invert the median transit
   time. The analysis plan uses (i) only if it fits within 1% and moves N* by more than 2%; otherwise (ii).
3. Find the ODE's safe threshold under the refitted law.
4. Run the size-only policy at m * N*_recal(tau), 500 seeds per cell from seed 1,800,500.
Outputs: p_t_recal_nstar.json and p_t_recal_odeb.json (steps 2-3 only), p_t_recal_sweep.npz (all steps).
"""
import os
for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
import json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import least_squares
from threshold_transfer.paths import RESULTS, WORKERS
from threshold_transfer.threshold import construction as P, grid as G
from threshold_transfer.experiments import recal_map as RC
from threshold_transfer.model.ode import GROWTH, N0, P as p, nstar, phi, regrow_time, ode_boundary

PROG, TAUS = RC.PROG, G.TAUS
C0, PHI0 = GROWTH


def dN(N, c, pb):
    """Untreated growth rate under constants (c, pb), cells per day."""
    return p.r_S * phi(N, c, pb) * N - p.base_death * N


def job(a):
    """One sweep run, a = (tau, m, thr, seed index). Returns (tau, m, s, prem, t_prog or -1, N_max)."""
    tau, m, thr, s = a
    r = P.run(RC.BLOCK + 500 + s, "AT", tau=tau, thr=thr, N0=N0, horizon=G.T_END)
    return (tau, m, s, r['prem'], -1.0 if np.isnan(r['t_prog']) else r['t_prog'], r['N_max'])


def main(no_sweep=False, NW=WORKERS):
    z = np.load(RESULTS / "p_t_recal_map.npz")
    MAPN = z['map_n']; TR = z['transit']; TRACES = z['traces']
    medT = np.array([np.median(v[np.isfinite(v)]) for v in TR])

    # (i) refit (c, phi_bulk) to the daily growth of the median trajectories
    Ns, dNs = [], []
    for i in range(len(MAPN)):
        med = np.nanmedian(TRACES[i], axis=0)
        med = med[np.isfinite(med)]
        for j in range(len(med) - 1):
            if med[j] >= PROG: break
            Ns.append(med[j]); dNs.append(med[j + 1] - med[j])
    Ns = np.array(Ns); dNs = np.array(dNs)
    res = least_squares(lambda th: dN(Ns, th[0], th[1]) - dNs, x0=[C0, PHI0], bounds=([0.0, -1.0], [50.0, 1.0]))
    C1, PHI1 = float(res.x[0]), float(res.x[1])
    rel = []                                   # relative RMSE of the refitted law against each median trajectory
    for i in range(len(MAPN)):
        med = np.nanmedian(TRACES[i], axis=0); med = med[np.isfinite(med)]
        h = 0.2; N = float(med[0]); pred = [N]
        for _ in range(len(med) - 1):
            for _ in range(int(1 / h)):
                k1 = dN(N, C1, PHI1); k2 = dN(N + h/2*k1, C1, PHI1); k3 = dN(N + h/2*k2, C1, PHI1); k4 = dN(N + h*k3, C1, PHI1)
                N = N + h/6*(k1 + 2*k2 + 2*k3 + k4)
            pred.append(N)
        pred = np.array(pred[:len(med)])
        rel.append(np.sqrt(np.mean(((pred - med) / med) ** 2)))
    NS_OLD = {t: nstar(float(t)) for t in TAUS}
    NS_FIT = {t: nstar(float(t), C1, PHI1) for t in TAUS}
    moved = max(abs(NS_FIT[t] / NS_OLD[t] - 1) for t in TAUS)
    clean = (max(rel) < 0.01) and (moved > 0.02)
    print(f"refit (c, phi_bulk) = ({C1:.4f}, {PHI1:.4f}); max relative RMSE {100 * max(rel):.3f}%; route (i) used: {clean}")

    # (ii) invert the median transit time
    order = np.argsort(medT)
    tt, nn = medT[order], MAPN[order].astype(float)
    NS_MAP = {t: float(np.interp(t, tt, nn)) if tt.min() <= t <= tt.max() else np.nan for t in TAUS}

    ROUTE = "i" if clean else "ii"
    NSR = {t: (NS_FIT[t] if ROUTE == "i" else NS_MAP[t]) for t in TAUS}
    ODEB = {t: ode_boundary(t, C1, PHI1)[0] for t in TAUS}
    for t in TAUS:
        print(f"tau={t}  N* {NS_OLD[t]:.1f}  N*_recal {NSR[t]:.1f}  ODE safe threshold {ODEB[t]:.1f}")
    print(f"untreated transit N0 -> 1.2 N0 under the refitted law: {regrow_time(N0, C1, PHI1):.2f} d")

    if no_sweep:
        json.dump({str(t): NSR[t] for t in TAUS}, open(RESULTS / "p_t_recal_nstar.json", "w"))
        json.dump({str(t): ODEB[t] for t in TAUS}, open(RESULTS / "p_t_recal_odeb.json", "w"))
        print("wrote p_t_recal_nstar.json and p_t_recal_odeb.json"); return

    GRIDR = {t: tuple(m for m in RC.M_GRID if m < PROG / NSR[t]) for t in TAUS}
    jobs = [(t, m, m * NSR[t], s) for t in TAUS for m in GRIDR[t] for s in range(G.R)]
    print(f"sweep: {len(jobs)} runs", flush=True)
    t0 = time.time()
    with Pool(NW) as pool:
        out = pool.map(job, jobs, chunksize=16)
    wall = time.time() - t0
    DT = [('tau', 'i4'), ('m', 'f8'), ('seed', 'i4'), ('prem', 'i4'), ('t_prog', 'f8'), ('nmax', 'f8')]
    np.savez_compressed(RESULTS / "p_t_recal_sweep.npz", rec=np.array(out, dtype=DT), route=ROUTE,
                        nstar_recal=np.array([NSR[t] for t in TAUS]),
                        nstar_fit=np.array([NS_FIT[t] for t in TAUS]),
                        nstar_map=np.array([NS_MAP[t] for t in TAUS]),
                        ode_boundary=np.array([ODEB[t] for t in TAUS]),
                        taus=np.array(TAUS), c_refit=C1, phi_refit=PHI1,
                        c_reg=C0, phi_reg=PHI0, block=RC.BLOCK)
    json.dump(dict(wall_s=wall, n=len(jobs), workers=NW, route=ROUTE), open(RESULTS / "p_t_recal_sweep_cost.json", "w"))
    print(f"wrote p_t_recal_sweep.npz ({wall / 60:.1f} min)")
