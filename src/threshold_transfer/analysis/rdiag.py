"""Resistant-count check: refits the size and history models with R added and scores on the test set the
history gain without R, the history gain with R in both models, and R alone over size.
Outputs: p_s_rdiag_models.pkl, p_s_rdiag_fit.json, p_s_rdiag.json.
"""
import json, pickle, time
from threshold_transfer.analysis import statesuff as A

OUT = A.OUTD


def main():
    """Fit, check and score the models. Replaces statesuff.MODELS and COMPS for this run."""
    assert not (OUT / "p_s_rdiag_fit.json").exists(), "already fitted; refusing to refit"
    A.MODELS.clear(); A.MODELS.update(A0=[], A1=A.HIST, A0R=["R"], A1R=A.HIST + ["R"])
    rtr, etr = A.load("train"); Dtr = A.hold_set(rtr, etr)
    F, conv = A.fit_all(Dtr, 5); h = A.model_hash(F)
    for tau in A.TAUS:   # the refitted size and history models must equal the frozen ones
        for m in ("A0", "A1"):
            assert h[f"{tau}|{m}"] == json.load(open(OUT / "p_s_fit.json"))["hash5"][f"{tau}|{m}"]
    pickle.dump(F, open(OUT / "p_s_rdiag_models.pkl", "wb"))
    json.dump(dict(t=time.strftime("%Y-%m-%d %H:%M:%S"), hash=h, converged=conv,
                   note="A0R/A1R frozen before any Step 2 block was read; A0/A1 reproduce Step 1 hashes"),
              open(OUT / "p_s_rdiag_fit.json", "w"), indent=1)
    rte, ete = A.load("test"); Dte = A.hold_set(rte, ete)
    P = A.predict(F, Dte)
    A.COMPS[:] = [("history gain, R NOT in baseline (A1 vs A0)", "A1", "A0", "diag"),
                  ("history gain, R in both (A1R vs A0R)", "A1R", "A0R", "diag"),
                  ("R alone over N+threshold (A0R vs A0)", "A0R", "A0", "diag")]
    res = A.score(A.events_by_tau(Dtr), Dte, P, list(range(500)), "R-CONDITIONING DIAGNOSTIC (TEST, frozen fits)")
    A.report(res)
    json.dump(A._json(res), open(OUT / "p_s_rdiag.json", "w"), indent=1)
    print("written p_s_rdiag.json")
