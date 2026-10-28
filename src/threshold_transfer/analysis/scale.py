"""Absolute scale of the prediction results: base rate, log-loss of a base-rate model, and log-loss and Brier
score of each compared model on the test set. Output: data/results/p_s_absolute_scale.json.
"""
import json, pickle
from pathlib import Path
import numpy as np
from threshold_transfer.analysis import statesuff as A


def main(out_path=None):
    """Compute the absolute-scale table and write it to out_path (default data/results/p_s_absolute_scale.json)."""
    out_path = Path(out_path) if out_path else A.OUTD / "p_s_absolute_scale.json"
    M = pickle.load(open(A.OUTD / "p_s_models.pkl", "rb")); fj = json.load(open(A.OUTD / "p_s_fit.json"))
    assert A.model_hash(M['F5']) == fj['hash5']
    R = json.load(open(A.OUTD / "p_s_results.json"))['primary']
    rte, ete = A.load('test'); D = A.hold_set(rte, ete); P = A.predict(M['F5'], D); y = D['y'].astype(float)
    out = {}
    def ll(p, yy): return float(np.mean(-(yy*np.log(p)+(1-yy)*np.log(1-p))))
    def br(p, yy): return float(np.mean((p-yy)**2))
    for cname, m1, m0, role in A.COMPS:
        C = R['comparisons'][cname]; rows = {}
        for tau in list(A.TAUS) + ['pool']:
            if tau == 'pool':
                p0 = np.concatenate([P[(t, m0)] for t in A.TAUS]); p1 = np.concatenate([P[(t, m1)] for t in A.TAUS])
                yy = np.concatenate([y[D['tau'] == t] for t in A.TAUS]); v = C['pooled']
            else:
                yy = y[D['tau'] == tau]; p0 = P[(tau, m0)]; p1 = P[(tau, m1)]; v = C['per_tau'][str(tau)]
            b = yy.mean(); lnull = float(-(b*np.log(b)+(1-b)*np.log(1-b)))
            r = dict(base_rate=float(b), n=int(len(yy)), events=int(yy.sum()), ll_null=lnull, ll_base=ll(p0, yy), ll_test=ll(p1, yy),
                     brier_null=float(b*(1-b)), brier_base=br(p0, yy), brier_test=br(p1, yy), delta_ll=v['delta'], delta_ll_ci=v['ci'])
            assert abs(1 - r['ll_test']/r['ll_base'] - v['delta']) < 1e-9
            st = v.get('stratified', {})
            if 'dauc_within' in st: r['dauc_within'] = st['dauc_within']; r['dauc_within_ci'] = st['ci']
            rows[str(tau)] = r
        out[cname] = rows
    json.dump(out, open(out_path, 'w'), indent=1)
    print("written", out_path)
