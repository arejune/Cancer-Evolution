"""Prediction study: logistic models of premature progression at each hold decision.

Rows are decisions after the first at which treatment was held and whose interval is not cut by the horizon;
y = 1 if the tumour reached 1.2 N0 before the next decision. Every model is a ridge logistic regression (C = 1)
on a cubic B-spline in N (5 knots at training quantiles), one indicator per threshold and a feature set
(MODELS). Each comparison (COMPS) reports the relative reduction in test log-loss with a seed bootstrap CI
(4,000 resamples), the AUC difference, and the AUC difference within 20-cell bins of N.

main("fit") fits on the training set and stores each model's hash; main("score") scores on the test set;
main("dryrun") runs both on the gate block as a code check.
"""
import json, hashlib, pickle, warnings, time
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import SplineTransformer
from sklearn.exceptions import ConvergenceWarning
from threshold_transfer.paths import RESULTS

CD = RESULTS / "p_s_cells"
OUTD = RESULTS
TAUS = (7, 14, 28, 56)
BAND = {7: (3060, 3120, 3180, 3240, 3300), 14: (2900, 2950, 3000, 3050, 3100),
        28: (2670, 2740, 2810, 2880, 2950), 56: (2660, 2740, 2820, 2900, 2980)}
NBOOT = 4000
BOOT_SEED = 20260923            # seed of the bootstrap weights
BAR = 0.10
MIN_EV = 50                     # fewer events than this: underpowered
BIN_W, BIN_MIN, BIN_NQ = 20, 10, 3   # within-bin AUC: bin width, events and non-events per bin, bins needed

HIST = ["tsf", "never_on", "g", "log_non", "nadir"]
GEO_IMG = ["B", "K", "Rg", "Rmax"]
GEO_FULL = ["B", "B_R", "K", "Rg", "Rmax"]
MODELS = {"A0": [], "A1": HIST, "A0c": ["C"], "A1c": HIST + ["C"],
          "S0": HIST + ["C", "R"], "S1F": HIST + ["C", "R"] + GEO_FULL,
          "S1I": HIST + ["C", "R"] + GEO_IMG,
          "S0o": HIST + ["C"], "S1Io": HIST + ["C"] + GEO_IMG,
          "CLIN": HIST + GEO_IMG}
COMPS = [("Q-ACT", "A1", "A0", "PRIMARY"),
         ("Q-ACT-carryover", "A1c", "A0c", "secondary"),
         ("Q-SP-FULL (R-IN)", "S1F", "S0", "PRIMARY"),
         ("Q-SP-IMAGE (R-IN)", "S1I", "S0", "secondary"),
         ("Q-SP-IMAGE (R-OUT)", "S1Io", "S0o", "secondary"),
         ("CLIN-IMAGE (no R, no C)", "CLIN", "A1", "descriptive")]


def load(role):
    """Runs and decision rows of one block: 'gate', 'train' or 'test'."""
    if role == "gate":
        z = np.load(CD / "gate.npz"); return z["runs"], z["eps"]
    fs = sorted((CD / role).glob("*.npz"))
    fs = [f for f in fs if not f.name.endswith(".tmp.npz")]
    runs = np.concatenate([np.load(f)["runs"] for f in fs])
    eps = np.concatenate([np.load(f)["eps"] for f in fs])
    assert len(runs) == 20 * 500, f"{role}: {len(runs)} runs, expected 10000"
    return runs, eps


def sanity(runs, eps):
    """Checks on a block: labels match outcomes, no events after treated decisions, no escape cells, edge distance."""
    key = lambda a: (a["tau"].astype(np.int64) * 10**8 + a["thr"].astype(np.int64) * 10**4
                     + a["s"].astype(np.int64))
    kr, ke = key(runs), key(eps)
    hold_y = np.bincount(np.searchsorted(np.sort(kr), ke), weights=eps["y"] * (eps["on"] == 0),
                         minlength=len(kr))
    order = np.argsort(kr)
    lab_ok = bool(np.all(hold_y[np.arange(len(kr))] == runs["prem"][order]))
    on_y = int(eps["y"][eps["on"] == 1].sum())
    return dict(label_identity=lab_ok, y_on_ON_epochs=on_y, escape_nonzero=int((eps["Esc"] != 0).sum()),
                min_d_wall=int(eps["d_wall"].min()), n_trunc=int(eps["trunc"].sum()))


def history(eps):
    """History features for each row of eps: days since treatment was last given (tsf), never treated
    (never_on), log(N / N_prev) (g), log(1 + earlier treated decisions) (log_non) and the lowest N since
    treatment was last given (nadir); plus C, R and the shape columns.
    """
    o = np.lexsort((eps["e"], eps["s"], eps["thr"], eps["tau"]))
    E = eps[o]
    n = len(E)
    tsf = np.zeros(n); never = np.zeros(n); g = np.full(n, np.nan)
    non = np.zeros(n); nadir = np.zeros(n)
    rid = (E["tau"].astype(np.int64) * 10**8 + E["thr"].astype(np.int64) * 10**4 + E["s"])
    starts = np.flatnonzero(np.r_[True, rid[1:] != rid[:-1]])
    ends = np.r_[starts[1:], n]
    for a, b in zip(starts, ends):
        N = E["N"][a:b].astype(float); on = E["on"][a:b]; tau = int(E["tau"][a])
        assert np.array_equal(E["e"][a:b], np.arange(b - a))
        last = -1; cnt = 0
        for i in range(b - a):
            if last >= 0:
                tsf[a+i] = (i - last) * tau; never[a+i] = 0
                nadir[a+i] = N[last+1:i+1].min() if i > last else N[i]
            else:
                tsf[a+i] = i * tau; never[a+i] = 1; nadir[a+i] = N[:i+1].min()
            non[a+i] = cnt
            if i >= 1: g[a+i] = np.log(N[i] / N[i-1])
            if on[i]: last = i; cnt += 1
    inv = np.empty(n, dtype=int); inv[o] = np.arange(n)
    H = dict(tsf=tsf[inv], never_on=never[inv], g=g[inv], log_non=np.log1p(non[inv]),
             nadir=nadir[inv])
    for f in ("C", "R", "B", "B_R", "K", "Rg", "Rmax"):
        H[f] = eps[f].astype(float)
    return H


def hold_set(runs, eps):
    """Rows used for modelling (held, e >= 1, not truncated) with their features, tau, thr, seed, N and y."""
    H = history(eps)
    m = (eps["on"] == 0) & (eps["e"] >= 1) & (eps["trunc"] == 0)
    D = {k: v[m] for k, v in H.items()}
    for k in ("tau", "thr", "s", "N", "y"):
        D[k] = eps[k][m].astype(float if k == "N" else int)
    return D


def design(D, idx, feats, tau, prep=None, n_knots=5):
    """Design matrix: B-spline of N, one-hot threshold cells (first cell dropped), then `feats`; standardised
    with the fit sample's mean and SD. `prep` from a previous call reuses that spline and scaling.
    """
    N = D["N"][idx].reshape(-1, 1)
    if prep is None:
        spl = SplineTransformer(n_knots=n_knots, degree=3, knots="quantile",
                                extrapolation="linear", include_bias=False).fit(N)
    else:
        spl = prep["spl"]
    cells = [np.asarray(D["thr"][idx] == t, float) for t in BAND[tau][1:]]
    X = np.column_stack([spl.transform(N)] + cells + [D[f][idx] for f in feats])
    if prep is None:
        mu, sd = X.mean(0), X.std(0); keep = sd > 0
        prep = dict(spl=spl, mu=mu, sd=sd, keep=keep)
    X = (X[:, prep["keep"]] - prep["mu"][prep["keep"]]) / prep["sd"][prep["keep"]]
    return X, prep


def fit_all(D, n_knots=5):
    """Fit every model in MODELS at every tau on the rows of D. Returns ({(tau, model): (prep, fit)}, converged flags)."""
    F = {}; conv = {}
    for tau in TAUS:
        idx = np.flatnonzero(D["tau"] == tau)
        for mname, feats in MODELS.items():
            X, prep = design(D, idx, feats, tau, n_knots=n_knots)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always", ConvergenceWarning)
                lr = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
                                        max_iter=10000, fit_intercept=True).fit(X, D["y"][idx])
            conv[f"{tau}|{mname}"] = not any(issubclass(x.category, ConvergenceWarning) for x in w)
            F[(tau, mname)] = (prep, lr)
    return F, conv


def model_hash(F):
    """SHA-256 of each fitted model's coefficients and intercept (float64 bytes)."""
    h = {}
    for (tau, m), (prep, lr) in sorted(F.items()):
        b = np.concatenate([lr.coef_.ravel(), lr.intercept_]).astype(np.float64).tobytes()
        h[f"{tau}|{m}"] = hashlib.sha256(b).hexdigest()
    return h


def predict(F, D):
    """Predicted risk p_hat, clipped to [1e-15, 1 - 1e-15], for every (tau, model)."""
    P = {}
    for tau in TAUS:
        idx = np.flatnonzero(D["tau"] == tau)
        for m, feats in MODELS.items():
            prep, lr = F[(tau, m)]
            X, _ = design(D, idx, feats, tau, prep=prep)
            P[(tau, m)] = np.clip(lr.predict_proba(X)[:, 1], 1e-15, 1 - 1e-15)
    return P


def boot_weights(nseed):
    """NBOOT x nseed multinomial weights for the seed-cluster bootstrap, from BOOT_SEED."""
    rng = np.random.default_rng(BOOT_SEED)
    return rng.multinomial(nseed, np.full(nseed, 1.0 / nseed), size=NBOOT).astype(float)


def wauc(p, y, Wseed, si):
    """Weighted AUC with mid-rank ties; one value per row of the seed weights Wseed; si is each row's seed index."""
    o = np.argsort(p, kind="mergesort"); ps = p[o]; ys = y[o]; so = si[o]
    _, gid = np.unique(ps, return_inverse=True); G = gid.max() + 1
    Wseed = np.atleast_2d(Wseed)
    out = np.empty(Wseed.shape[0])
    for b in range(Wseed.shape[0]):
        w = Wseed[b][so]
        pos = np.bincount(gid, w * ys, G); neg = np.bincount(gid, w * (1 - ys), G)
        cn = np.cumsum(neg) - neg
        den = pos.sum() * neg.sum()
        out[b] = (pos * (cn + 0.5 * neg)).sum() / den if den > 0 else np.nan
    return out


def score(Dfit_events, D, P, seeds, label, do_auc=True):
    """Every comparison in COMPS per tau and pooled, scored on D with predictions P; seeds are the scoring
    block's seed indices. Returns a dict with per-tau counts, comparisons and calibration.
    """
    nseed = len(seeds); smap = {s: i for i, s in enumerate(seeds)}
    Wb = boot_weights(nseed)
    si = np.array([smap[s] for s in D["s"]])
    y = D["y"].astype(float)
    res = {"label": label, "per_tau": {}, "comparisons": {}}
    L = {}                     # (tau, model) -> per-seed loss sums
    for tau in TAUS:
        idx = np.flatnonzero(D["tau"] == tau)
        ev_s = int(y[idx].sum()); ev_f = Dfit_events[tau]
        res["per_tau"][tau] = dict(hold_epochs=len(idx), events_score=ev_s, events_fit=ev_f,
                                   underpowered=bool(ev_s < MIN_EV or ev_f < MIN_EV))
        for m in MODELS:
            p = P[(tau, m)]
            ll = -(y[idx] * np.log(p) + (1 - y[idx]) * np.log(1 - p))
            L[(tau, m)] = np.bincount(si[idx], ll, nseed)
    for cname, m1, m0, role in COMPS:
        C = {"role": role, "test_model": m1, "baseline": m0, "per_tau": {}}
        tot1 = np.zeros(nseed); tot0 = np.zeros(nseed)
        for tau in TAUS:
            l1, l0 = L[(tau, m1)], L[(tau, m0)]
            tot1 += l1; tot0 += l0
            C["per_tau"][tau] = _delta_block(l1, l0, Wb, res["per_tau"][tau]["underpowered"])
        C["pooled"] = _delta_block(tot1, tot0, Wb, False)
        C["headline"] = _headline(C)
        res["comparisons"][cname] = C
    if do_auc:
        _auc_and_strata(res, D, P, si, Wb, y)
    _calibration(res, D, P, y)
    return res


def _delta_block(l1, l0, Wb, underpowered):
    """Delta = 1 - LL_test / LL_baseline with its bootstrap CI and reading."""
    d = 1 - l1.sum() / l0.sum()
    bd = 1 - (Wb @ l1) / (Wb @ l0)
    lo, hi = np.percentile(bd, [2.5, 97.5])
    if underpowered: rd = "UNDERPOWERED"
    elif d >= BAR and lo > 0: rd = "GAIN"
    elif hi < BAR: rd = "NULL" + (" (immaterial gain detected)" if lo > 0 else "")
    else: rd = "INCONCLUSIVE"
    return dict(delta=float(d), ci=[float(lo), float(hi)], reading=rd,
                ll_test=float(l1.sum()), ll_base=float(l0.sum()))


def _headline(C):
    """Pooled result: it stands if the pooled reading is GAIN or NULL and two intervals share its reading and sign."""
    pr = C["pooled"]["reading"].split(" ")[0]; ps = np.sign(C["pooled"]["delta"])
    sup = [t for t, v in C["per_tau"].items()
           if v["reading"].split(" ")[0] == pr and np.sign(v["delta"]) == ps]
    ok = pr in ("GAIN", "NULL") and len(sup) >= 2
    return dict(pooled_reading=pr, supporting_tau=sup,
                stands=bool(ok), text=(f"pooled {pr} stands (tau {sup})" if ok else
                                       f"no pooled headline (pooled {pr}; same-reading same-sign tau = {sup})"))


def _auc_and_strata(res, D, P, si, Wb, y):
    """Adds AUC, Delta AUC with CI, and the within-bin Delta AUC at matched size to each comparison."""
    _CACHE = {}
    one = np.ones((1, Wb.shape[1]))
    for cname, m1, m0, _ in COMPS:
        C = res["comparisons"][cname]
        for tau in TAUS:
            idx = np.flatnonzero(D["tau"] == tau); s_ = si[idx]; y_ = y[idx]
            p1, p0 = P[(tau, m1)], P[(tau, m0)]
            a1 = wauc(p1, y_, one, s_)[0]; a0 = wauc(p0, y_, one, s_)[0]
            for m in (m1, m0):
                if (tau, m) not in _CACHE: _CACHE[(tau, m)] = wauc(P[(tau, m)], y_, Wb, s_)
            b = _CACHE[(tau, m1)] - _CACHE[(tau, m0)]
            C["per_tau"][tau].update(auc_test=float(a1), auc_base=float(a0), dauc=float(a1 - a0),
                                     dauc_ci=[float(x) for x in np.nanpercentile(b, [2.5, 97.5])])
            # AUC difference within 20-cell bins of N, weighted by events
            bins = np.floor(D["N"][idx] / BIN_W).astype(int)
            q = [bb for bb in np.unique(bins)
                 if y_[bins == bb].sum() >= BIN_MIN and (1 - y_[bins == bb]).sum() >= BIN_MIN]
            st = dict(n_qualifying_bins=len(q), bins=[int(bb * BIN_W) for bb in q])
            if len(q) < BIN_NQ:
                st["reading"] = "UNDERPOWERED"
            else:
                num = 0.0; den = 0.0; bn = np.zeros(Wb.shape[0]); bdn = np.zeros(Wb.shape[0])
                for bb in q:
                    j = bins == bb; yy = y_[j]; sj = s_[j]
                    d_pt = wauc(p1[j], yy, one, sj)[0] - wauc(p0[j], yy, one, sj)[0]
                    num += yy.sum() * d_pt; den += yy.sum()
                    dB = wauc(p1[j], yy, Wb, sj) - wauc(p0[j], yy, Wb, sj)
                    evB = Wb @ np.bincount(sj, yy, Wb.shape[1])
                    ok = ~np.isnan(dB)
                    bn += np.where(ok, dB, 0.0) * evB; bdn += evB * ok
                est = num / den; bb_ = bn / bdn
                lo, hi = np.nanpercentile(bb_, [2.5, 97.5])
                st.update(dauc_within=float(est), ci=[float(lo), float(hi)],
                          reading="separates" if lo > 0 else "CI includes 0")
            C["per_tau"][tau]["stratified"] = st


def _calibration(res, D, P, y):
    """Observed / expected events, calibration slope and decile reliability per (tau, model)."""
    cal = {}
    for tau in TAUS:
        idx = np.flatnonzero(D["tau"] == tau)
        for m in MODELS:
            p = P[(tau, m)]; lg = np.log(p / (1 - p)).reshape(-1, 1)
            try:
                s = LogisticRegression(penalty=None, max_iter=10000).fit(lg, y[idx]).coef_[0, 0]
            except Exception:
                s = np.nan
            dec = np.quantile(p, np.linspace(0, 1, 11)); bi = np.clip(np.searchsorted(dec, p, "right") - 1, 0, 9)
            rel = [[float(p[bi == k].mean()), float(y[idx][bi == k].mean())] if (bi == k).any() else None
                   for k in range(10)]
            cal[f"{tau}|{m}"] = dict(O_over_E=float(y[idx].sum() / p.sum()), slope=float(s), reliability=rel)
    res["calibration"] = cal


def events_by_tau(D):
    """Number of events y = 1 per tau."""
    return {tau: int(D["y"][D["tau"] == tau].sum()) for tau in TAUS}


def cell_table(runs):
    """p_prem and run count per (tau, thr) cell."""
    return {f"{int(t)}|{int(th)}": dict(p_prem=float(runs["prem"][(runs["tau"] == t) & (runs["thr"] == th)].mean()),
                                         n=int(((runs["tau"] == t) & (runs["thr"] == th)).sum()))
            for t in TAUS for th in BAND[t]}


def report(res):
    """Print a scoring result."""
    print(f"\n==== {res['label']} ====")
    for tau, v in res["per_tau"].items():
        print(f"  tau={tau:2d} hold epochs {v['hold_epochs']:6d}  events score/fit "
              f"{v['events_score']}/{v['events_fit']}{'  UNDERPOWERED' if v['underpowered'] else ''}")
    for cname, C in res["comparisons"].items():
        print(f"  [{C['role']}] {cname}: {C['test_model']} vs {C['baseline']}")
        for tau, v in list(C["per_tau"].items()) + [("pool", C["pooled"])]:
            s = f"    tau={tau!s:>4}  Δ={100*v['delta']:+7.2f}%  CI[{100*v['ci'][0]:+7.2f},{100*v['ci'][1]:+7.2f}]  {v['reading']}"
            if "dauc" in v:
                st = v["stratified"]
                s += f"  | ΔAUC {v['dauc']:+.4f}  strat: {st['n_qualifying_bins']} bins"
                if "dauc_within" in st:
                    s += f" ΔAUC_w {st['dauc_within']:+.4f} [{st['ci'][0]:+.4f},{st['ci'][1]:+.4f}]"
                else:
                    s += f" {st['reading']}"
            print(s)
        print(f"    HEADLINE: {C['headline']['text']}")


def _json(o):
    """Convert numpy scalars and tuple keys for JSON output."""
    if isinstance(o, dict): return {str(k): _json(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_json(v) for v in o]
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, np.bool_): return bool(o)
    return o


def main(mode):
    """Run one mode: 'dryrun', 'fit' or 'score'."""
    if mode == "dryrun":
        runs, eps = load("gate"); print("sanity", sanity(runs, eps))
        D = hold_set(runs, eps)
        fit = {k: v[D["s"] < 5] for k, v in D.items()}; sc = {k: v[D["s"] >= 5] for k, v in D.items()}
        F, conv = fit_all(fit); P = predict(F, sc)
        res = score(events_by_tau(fit), sc, P, list(range(5, 10)), "DRY RUN (gate block, NOT a result)")
        report(res); print("converged:", all(conv.values()))
    elif mode == "fit":
        assert not (OUTD / "p_s_models.pkl").exists(), "fit already done; refusing to refit"
        runs, eps = load("train"); san = sanity(runs, eps); print("TRAIN sanity", san)
        D = hold_set(runs, eps)
        F5, c5 = fit_all(D, 5); F8, c8 = fit_all(D, 8)
        h5, h8 = model_hash(F5), model_hash(F8)
        pickle.dump(dict(F5=F5, F8=F8, ev=events_by_tau(D)), open(OUTD / "p_s_models.pkl", "wb"))
        json.dump(_json(dict(t=time.strftime("%Y-%m-%d %H:%M:%S"), sanity=san, conv5=c5, conv8=c8,
                             hash5=h5, hash8=h8, cells_train=cell_table(runs), events_train=events_by_tau(D),
                             hold_train=int(len(D["y"])))),
                  open(OUTD / "p_s_fit.json", "w"), indent=1)
        print("fit written; non-converged:", [k for k, v in {**c5, **c8}.items() if not v])
    elif mode == "score":
        M = pickle.load(open(OUTD / "p_s_models.pkl", "rb"))
        fj = json.load(open(OUTD / "p_s_fit.json"))
        assert model_hash(M["F5"]) == fj["hash5"] and model_hash(M["F8"]) == fj["hash8"]
        rtr, etr = load("train"); Dtr = hold_set(rtr, etr)
        rte, ete = load("test"); san = sanity(rte, ete); print("TEST sanity", san)
        Dte = hold_set(rte, ete)
        seeds = list(range(500))
        out = dict(test_sanity=san, cells_test=cell_table(rte))
        out["primary"] = score(M["ev"], Dte, predict(M["F5"], Dte), seeds, "PRIMARY: fit TRAIN, score TEST, 5 knots")
        report(out["primary"])
        out["knots8"] = score(M["ev"], Dte, predict(M["F8"], Dte), seeds, "SENSITIVITY: 8 knots", do_auc=False)
        report(out["knots8"])
        Fx, cx = fit_all(Dte, 5)
        out["crossfit"] = score(events_by_tau(Dte), Dtr, predict(Fx, Dtr), seeds,
                                "CROSS-FIT: fit TEST, score TRAIN", do_auc=False)
        out["crossfit_converged"] = cx
        report(out["crossfit"])
        json.dump(_json(out), open(OUTD / "p_s_results.json", "w"), indent=1)
        print("\nwritten p_s_results.json")
