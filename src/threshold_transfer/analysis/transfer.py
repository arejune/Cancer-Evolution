"""Printed report of the main run: p_prem per cell, the safe threshold N_safe per interval with a seed
bootstrap CI, p_prem at m = 1.00 in each model, and the division counters.
"""
import json
import numpy as np
from threshold_transfer.paths import RESULTS
from threshold_transfer.threshold import grid as G


def main():
    """Print the main-run report."""
    Q = G.Q; B = 4000
    z = np.load(RESULTS / "p_t_main_raw.npz", allow_pickle=False)
    rec, led = z['rec'], z['led']
    NS = {int(t): float(n) for t, n in zip(z['taus'], z['nstar'])}
    GRID = G.grid(NS)
    arms = json.load(open(RESULTS / "p_t_main_arms.json"))


    def wilson(k, n, zz=1.96):
        """Wilson score interval (Wilson 1927) for k events in n runs at z = 1.96."""
        if n == 0: return (np.nan, np.nan)
        p = k / n; den = 1 + zz*zz/n; c = (p + zz*zz/(2*n))/den
        h = zz*np.sqrt(p*(1-p)/n + zz*zz/(4*n*n))/den
        return (max(0., c-h), min(1., c+h))


    def prem_matrix(tau):
        """(R, number of m) matrix of premature-progression indicators at tau, rows ordered by seed."""
        ms = GRID[tau]; M = np.zeros((G.R, len(ms)))
        for j, m in enumerate(ms):
            s = rec[(rec['pol'] == "AT") & (rec['tau'] == tau) & (np.isclose(rec['m'], m))]
            o = np.argsort(s['seed']); M[:, j] = s['prem'][o]
        return ms, M


    def interp_safe(ms, p):
        """Largest m with p_prem <= Q, linear in p_prem between cells; p_prem is made non-decreasing in m first."""
        p = np.maximum.accumulate(np.asarray(p, float)); ms = np.asarray(ms, float)
        if p[0] > Q: return -np.inf
        if p[-1] <= Q: return np.inf
        i = int(np.argmax(p > Q))
        p0, p1 = p[i-1], p[i]
        if p1 <= p0: return ms[i-1]
        return ms[i-1] + (Q - p0) * (ms[i] - ms[i-1]) / (p1 - p0)


    print("=" * 96)
    print("PER-CELL OUTCOMES  (ABM, block 1200000+, R = 500, t_end = 750 d, q = 0.05)")
    print("=" * 96)
    print(f"{'tau':>4} {'m':>5} {'thr':>8} {'p_prem':>7} {'95% Wilson CI':>18} "
          f"{'cens':>5} {'mean t_prog':>12} {'med':>6}")
    safe = {}
    for t in G.TAUS:
        ms, M = prem_matrix(t)
        for m in ms:
            s = rec[(rec['pol'] == "AT") & (rec['tau'] == t) & (np.isclose(rec['m'], m))]
            k, n = int(s['prem'].sum()), len(s)
            lo, hi = wilson(k, n); tp = s['t_prog']; cens = int((tp < 0).sum())
            d = tp[tp > 0]
            print(f"{t:4d} {m:5.2f} {m*NS[t]:8.1f} {k/n:7.4f} [{lo:6.4f},{hi:6.4f}] "
                  f"{cens:5d} {d.mean() if len(d) else np.nan:12.1f} {np.median(d) if len(d) else np.nan:6.0f}")
        print()
    s = rec[rec['pol'] == "CT"]
    d = s['t_prog'][s['t_prog'] > 0]
    print(f"  CT (continuous)  R={len(s)}  censored {int((s['t_prog']<0).sum())}  "
          f"mean t_prog {d.mean() if len(d) else np.nan:.1f}  median {np.median(d) if len(d) else np.nan:.0f}")

    print("\n" + "=" * 96)
    print("EMPIRICAL SAFE THRESHOLD  m_safe(tau) AND RELATIVE DISCREPANCY FROM THE ANALYTIC N*")
    print("bootstrap: 4000 resamples of the 500 SEEDS (cluster unit = seed, CRN across cells)")
    print("=" * 96)
    print(f"{'tau':>4} {'m_safe':>8} {'95% CI':>20} {'rel. discrepancy':>18} {'95% CI':>20} {'N*_safe':>9}")
    rng = np.random.default_rng(20260822)
    for t in G.TAUS:
        ms, M = prem_matrix(t)
        pt = interp_safe(ms, M.mean(0)); safe[t] = pt
        idx = rng.integers(0, G.R, size=(B, G.R))
        bs = np.array([interp_safe(ms, M[i].mean(0)) for i in idx])
        fin = np.isfinite(bs)
        lo, hi = (np.percentile(bs[fin], 2.5), np.percentile(bs[fin], 97.5)) if fin.sum() > 10 else (np.nan, np.nan)
        rd, rlo, rhi = pt - 1.0, lo - 1.0, hi - 1.0
        tag = "" if fin.all() else f"  [{100*(~fin).mean():.1f}% resamples unbracketed]"
        print(f"{t:4d} {pt:8.4f} [{lo:8.4f},{hi:8.4f}] {100*rd:+16.2f}% "
              f"[{100*rlo:+8.2f}%,{100*rhi:+8.2f}%] {pt*NS[t]:9.1f}{tag}")

    print("\n" + "=" * 96)
    print("DECOMPOSITION AT m = 1.00 (the analytic prediction), p_prem BY ARM")
    print("ODE = deterministic | CTMC-nom = well-mixed stochastic, nominal C_eff")
    print("CTMC-pk = well-mixed stochastic + drug field matched to the ABM | ABM = spatial")
    print("=" * 96)
    print(f"{'tau':>4} {'ODE':>6} {'CTMCnom':>9} {'CTMCpk':>8} {'ABM':>7} | "
          f"{'ODE->CTMC':>10} {'CTMC->ABM':>10} | {'% stoch':>8} {'% spatial':>10}")
    dec = {}
    for t in G.TAUS:
        key = f"|{t}|1.00"
        o = float(arms['cells']['ode' + key]['prem'])
        cn = float(np.mean(arms['cells']['ctmc_nom' + key]['prem']))
        cp = float(np.mean(arms['cells']['ctmc_pk' + key]['prem']))
        s = rec[(rec['pol'] == "AT") & (rec['tau'] == t) & (np.isclose(rec['m'], 1.00))]
        ab = float(s['prem'].mean())
        d1, d2 = cn - o, ab - cn; tot = ab - o
        f1 = 100*d1/tot if tot > 1e-12 else np.nan
        f2 = 100*d2/tot if tot > 1e-12 else np.nan
        dec[t] = (o, cn, cp, ab, d1, d2, f1, f2)
        print(f"{t:4d} {o:6.3f} {cn:9.4f} {cp:8.4f} {ab:7.4f} | {d1:+10.4f} {d2:+10.4f} | "
              f"{f1:7.1f}% {f2:9.1f}%")
    print("\nFINER SPLIT (declared before the run): CTMC-nom -> CTMC-pk isolates the drug")
    print("field's spatial transport; CTMC-pk -> ABM is contact limitation + lattice discreteness.")
    print(f"{'tau':>4} {'stochastic':>11} {'drug transport':>15} {'contact/lattice':>16} {'total':>8}")
    for t in G.TAUS:
        o, cn, cp, ab, *_ = dec[t]
        print(f"{t:4d} {cn-o:+11.4f} {cp-cn:+15.4f} {ab-cp:+16.4f} {ab-o:+8.4f}")

    print("\n" + "=" * 96)
    print("COMPETITION LEDGER AT FULL R  (instrumented cells: m = 1.00 column + CT)")
    print("=" * 96)
    sel = (np.isclose(rec['m'], 1.00)) | (rec['pol'] == "CT")
    print(f"{'cell':>14} {'runs':>5} {'attempts':>12} {'blocked':>12} {'blk%':>6} | "
          f"{'S-div intra':>11} {'S-div inter':>11} | {'R-div intra':>11} {'R-div inter':>11}")
    for t in list(G.TAUS) + ["CT"]:
        if t == "CT": mk = rec['pol'] == "CT"; lbl = "CT"
        else: mk = sel & (rec['tau'] == t) & (rec['pol'] == "AT"); lbl = f"AT tau={t}"
        if mk.sum() == 0: continue
        L = led[mk].sum(0).reshape(2, 6)
        att = L[:, 0].sum(); blk = L[:, 2].sum()
        sS, sR = L[0, 4], L[0, 5]        # S divider: slots held by S / by Ra
        rS, rR = L[1, 4], L[1, 5]        # Ra divider: slots held by S / by Ra
        print(f"{lbl:>14} {int(mk.sum()):5d} {att:12,.0f} {blk:12,.0f} {100*blk/att:5.1f}% | "
              f"{100*sS/max(sS+sR,1):10.2f}% {100*sR/max(sS+sR,1):10.2f}% | "
              f"{100*rR/max(rS+rR,1):10.2f}% {100*rS/max(rS+rR,1):10.2f}%")
        if t == 28:
            print(f"{'':>14} n slots: S-divider {int(sS+sR):,}   Ra-divider {int(rS+rR):,}   "
                  f"wall {int(L[:,3].sum()):,}")
