"""Growth-rate prediction model: size spline + threshold indicators + dN, the change in N since the previous
decision.

main("fit") fits it on the training set and stores its hash; main("score") scores it on the test set and
reports the share of the history model's log-loss gain over the size model that it recovers.
"""
import json, pickle, time, warnings
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning
from threshold_transfer.analysis import statesuff as A

OUTD = A.OUTD


def with_dN(D):
    """Add dN = N - N_prev, recovering N_prev from the growth ratio g (checked to be a whole number)."""
    Nprev = D["N"] * np.exp(-D["g"])
    r = np.round(Nprev)
    assert np.max(np.abs(Nprev - r)) < 1e-6
    D["dN"] = D["N"] - r
    return D


def fit_vel(D):
    """Fit the growth-rate model at each interval. Returns ({(tau, 'VEL'): (prep, fit)}, converged flags)."""
    F = {}; conv = {}
    for tau in A.TAUS:
        idx = np.flatnonzero(D["tau"] == tau)
        X, prep = A.design(D, idx, ["dN"], tau, n_knots=5)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", ConvergenceWarning)
            lr = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
                                    max_iter=10000, fit_intercept=True).fit(X, D["y"][idx])
        conv[tau] = not any(issubclass(x.category, ConvergenceWarning) for x in w)
        F[(tau, "VEL")] = (prep, lr)
    return F, conv


def main(mode):
    """mode 'fit' (training set) or 'score' (test set)."""
    if mode == "fit":
        assert not (OUTD / "p_s_vel_model.pkl").exists(), "R-VEL already fitted; refusing to refit"
        runs, eps = A.load("train"); D = with_dN(A.hold_set(runs, eps))
        F, conv = fit_vel(D); h = A.model_hash(F)
        pickle.dump(F, open(OUTD / "p_s_vel_model.pkl", "wb"))
        json.dump(dict(t=time.strftime("%Y-%m-%d %H:%M:%S"), converged={str(k): v for k, v in conv.items()},
                       hash=h, note="fitted before any Step 2 block (2400000/2500000/2600000) was read"),
                  open(OUTD / "p_s_vel_fit.json", "w"), indent=1)
        print("R-VEL fit written", h, "converged", conv)
    elif mode == "score":
        FV = pickle.load(open(OUTD / "p_s_vel_model.pkl", "rb"))
        assert A.model_hash(FV) == json.load(open(OUTD / "p_s_vel_fit.json"))["hash"]
        M = pickle.load(open(OUTD / "p_s_models.pkl", "rb"))
        assert A.model_hash(M["F5"]) == json.load(open(OUTD / "p_s_fit.json"))["hash5"]
        rte, ete = A.load("test"); D = with_dN(A.hold_set(rte, ete))
        P = A.predict(M["F5"], D); A.MODELS["VEL"] = ["dN"]   # VEL scored alongside the frozen models
        for tau in A.TAUS:
            idx = np.flatnonzero(D["tau"] == tau)
            prep, lr = FV[(tau, "VEL")]
            X, _ = A.design(D, idx, ["dN"], tau, prep=prep)
            P[(tau, "VEL")] = np.clip(lr.predict_proba(X)[:, 1], 1e-15, 1 - 1e-15)
        A.COMPS[:] = [("R-VEL vs Q-ACT baseline", "VEL", "A0", "amendment-1"),
                      ("R-HIST (A1) vs Q-ACT baseline", "A1", "A0", "reference")]
        res = A.score(A.events_by_tau(D), D, P, list(range(500)), "R-VEL prediction-level (TEST, frozen fits)")
        # recovery = (LL_A0 - LL_VEL) / (LL_A0 - LL_A1), with a seed bootstrap
        y = D["y"].astype(float); si = D["s"]; Wb = A.boot_weights(500); rec = {}
        tot = {m: np.zeros(500) for m in ("A0", "A1", "VEL")}
        for tau in list(A.TAUS) + ["pool"]:
            if tau == "pool":
                L = tot
            else:
                idx = np.flatnonzero(D["tau"] == tau); L = {}
                for m in ("A0", "A1", "VEL"):
                    p = P[(tau, m)]
                    L[m] = np.bincount(si[idx], -(y[idx]*np.log(p) + (1-y[idx])*np.log(1-p)), 500)
                    tot[m] += L[m]
            pt = (L["A0"].sum() - L["VEL"].sum()) / (L["A0"].sum() - L["A1"].sum())
            b = (Wb @ L["A0"] - Wb @ L["VEL"]) / (Wb @ L["A0"] - Wb @ L["A1"])
            n = int(len(y)) if tau == "pool" else int(len(idx))
            rec[str(tau)] = dict(recovery=float(pt), ci=[float(x) for x in np.percentile(b, [2.5, 97.5])],
                                 ll_A0=float(L["A0"].sum() / n), ll_A1=float(L["A1"].sum() / n),
                                 ll_VEL=float(L["VEL"].sum() / n))
        res["loglik_recovery_of_A1_gain"] = rec
        A.report(res)
        print("\nrecovery of R-HIST's log-loss gain by R-VEL:")
        for k, v in rec.items():
            print(f"  tau={k:>4}  LL A0/VEL/A1 = {v['ll_A0']:.5f} / {v['ll_VEL']:.5f} / {v['ll_A1']:.5f}   "
                  f"recovery {100*v['recovery']:.1f}% [{100*v['ci'][0]:.1f}, {100*v['ci'][1]:.1f}]")
        json.dump(A._json(res), open(OUTD / "p_s_vel_prediction.json", "w"), indent=1)
        print("written p_s_vel_prediction.json")
