"""Printed report of the recalibration: p_prem per cell, the safe threshold with a seed bootstrap CI, and the
derivation, transfer and total terms before and after recalibration.
"""
import numpy as np
from threshold_transfer.paths import RESULTS
from threshold_transfer.threshold import grid as G


def main():
    """Print the recalibration report."""
    Q, B, R = G.Q, 4000, G.R
    z  = np.load(RESULTS/"p_t_recal_sweep.npz", allow_pickle=False)
    rec = z['rec']; TAUS = z['taus'].tolist()
    NSR = dict(zip(TAUS, z['nstar_recal'])); ODEB = dict(zip(TAUS, z['ode_boundary']))
    zz = np.load(RESULTS/"p_t_main_raw.npz"); NS_OLD = dict(zip(zz['taus'].tolist(), zz['nstar']))
    # main-run reference values, recomputed from p_t_main_raw.npz
    GO = G.grid(NS_OLD)
    def isafe(ms, p, q=Q):
        ms = np.asarray(ms, float); p = np.maximum.accumulate(np.asarray(p, float))
        if p[0] > q: return -np.inf
        if p[-1] <= q: return np.inf
        i = int(np.argmax(p > q)); p0, p1 = p[i-1], p[i]
        return ms[i-1] if p1 <= p0 else ms[i-1] + (q-p0)*(ms[i]-ms[i-1])/(p1-p0)
    def old_mat(t):
        ms = GO[t]; M = np.zeros((R, len(ms)))
        for j, m in enumerate(ms):
            s = zz['rec'][(zz['rec']['pol']=="AT") & (zz['rec']['tau']==t) & (np.isclose(zz['rec']['m'],m))]
            M[:, j] = s['prem'][np.argsort(s['seed'])]
        return np.array(ms,float), M
    def new_mat(t):
        ms = sorted(set(rec['m'][rec['tau']==t]))
        M = np.zeros((R, len(ms)))
        for j, m in enumerate(ms):
            s = rec[(rec['tau']==t) & (np.isclose(rec['m'],m))]
            M[:, j] = s['prem'][np.argsort(s['seed'])]
        return np.array(ms,float), M

    print("="*100); print("RECALIBRATION CONTROL — BLOCK 1800000, ROUTE (ii)"); print("="*100)
    print(f"\n{'tau':>4} {'m':>5} {'thr':>8} {'p_prem':>7} {'k/n':>9} {'mean t_prog':>12}")
    for t in TAUS:
        ms, _ = new_mat(t)
        for m in ms:
            s = rec[(rec['tau']==t) & (np.isclose(rec['m'],m))]
            d = s['t_prog'][s['t_prog']>0]
            print(f"{t:4d} {m:5.2f} {m*NSR[t]:8.1f} {s['prem'].mean():7.4f} "
                  f"{int(s['prem'].sum()):4d}/{len(s):<4d} {d.mean() if len(d) else np.nan:12.1f}")
        print()

    rng = np.random.default_rng(20260826); idx = rng.integers(0, R, size=(B, R))
    print("="*100); print("SAFE THRESHOLD ON THE RECALIBRATED GRID, AND THE REGISTERED POINT PREDICTION")
    print("="*100)
    print(f"{'tau':>4} {'N*_recal':>9} {'m_safe':>8} {'N*_safe(new)':>13} {'95% CI':>22} "
          f"{'N*_safe(main)':>14} {'diff':>8} {'disc vs N*_recal':>17}")
    NEW, OLD = {}, {}
    for t in TAUS:
        ms, M = new_mat(t); v = isafe(ms, M.mean(0)); NEW[t] = v*NSR[t]
        bs = np.array([isafe(ms, M[i].mean(0)) for i in idx])*NSR[t]; bs = bs[np.isfinite(bs)]
        mso, Mo = old_mat(t); OLD[t] = isafe(mso, Mo.mean(0))*NS_OLD[t]
        lo, hi = np.percentile(bs,2.5), np.percentile(bs,97.5)
        print(f"{t:4d} {NSR[t]:9.1f} {v:8.4f} {NEW[t]:13.1f} [{lo:9.1f},{hi:9.1f}] {OLD[t]:14.1f} "
              f"{NEW[t]-OLD[t]:+8.1f} {100*(v-1):+16.2f}%")

    print("\n"+"="*100); print("THE SPLIT, RECOMPUTED UNDER THE RECALIBRATED LAW"); print("="*100)
    print(f"{'tau':>4} | {'BEFORE (registered law)':>34} | {'AFTER (recalibrated law)':>34}")
    print(f"{'':>4} | {'deriv':>10} {'transfer':>10} {'total':>10} | {'deriv':>10} {'transfer':>10} {'total':>10}")
    ODE_OLD = {7:3454.0, 14:3299.4, 28:3001.0, 42:2787.7, 56:2831.6, 84:2566.2}
    surv = {}
    for t in TAUS:
        d0 = ODE_OLD[t]/NS_OLD[t]-1; x0 = OLD[t]/ODE_OLD[t]-1; t0 = OLD[t]/NS_OLD[t]-1
        d1 = ODEB[t]/NSR[t]-1;      x1 = NEW[t]/ODEB[t]-1;     t1 = NEW[t]/NSR[t]-1
        surv[t] = (t0, t1, d0, d1, x0, x1)
        print(f"{t:4d} | {100*d0:+9.2f}% {100*x0:+9.2f}% {100*t0:+9.2f}% | "
              f"{100*d1:+9.2f}% {100*x1:+9.2f}% {100*t1:+9.2f}%")
    print("\nSURVIVING FRACTION OF THE ORIGINAL DISCREPANCY, per tau (prereg R4 PARTIAL reading)")
    print(f"{'tau':>4} {'total before':>13} {'total after':>12} {'fraction':>10}   note")
    for t in TAUS:
        t0, t1, d0, d1, x0, x1 = surv[t]
        note = "denominator too small to resolve" if abs(t0) < 0.02 else ""
        fr = f"{t1/t0:9.2f}" if abs(t0) >= 0.02 else f"{'--':>9}"
        print(f"{t:4d} {100*t0:+12.2f}% {100*t1:+11.2f}% {fr}   {note}")
    print("\nTERM-BY-TERM SURVIVAL")
    print(f"{'tau':>4} {'deriv before':>13} {'deriv after':>12} {'kept':>7} | "
          f"{'xfer before':>12} {'xfer after':>11} {'kept':>7}")
    for t in TAUS:
        t0, t1, d0, d1, x0, x1 = surv[t]
        print(f"{t:4d} {100*d0:+12.2f}% {100*d1:+11.2f}% {d1/d0:7.2f} | "
              f"{100*x0:+11.2f}% {100*x1:+10.2f}% {x1/x0:7.2f}")
