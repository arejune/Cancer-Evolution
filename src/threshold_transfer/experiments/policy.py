"""Policy test: hold rules driven by the frozen prediction models against the size-only rule.

Arms: SIZE holds if N < thr. The rules CLIN, HIST, FULL and VEL hold at a decision after the first if
N <= N_top(tau) and the model's predicted risk is at most p*; at the first decision they copy SIZE. The run
loop is construction.run plus state reads and model evaluations, which draw no random numbers.

Seed blocks: gate 2,600,000 (10 seeds), tuning 2,400,000 (500), scoring 2,500,000 (1,000). Tuning is a
bisection (8 evaluations of 500 seeds) for the most permissive setting with p_prem <= q; it resumes from the
stored files. Scoring compares each rule with SIZE on the same seeds (on-drug fraction, RMST over 750 days,
p_prem) with a seed bootstrap. Output: data/results/p_s_policy/ and p_s_policy_results.json.
"""
import os
for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "1"
import sys, time, json, pickle, hashlib
import numpy as np
from multiprocessing import Pool
from threshold_transfer.paths import RESULTS, WORKERS
from threshold_transfer.threshold import construction as P
from threshold_transfer.experiments import statesuff_log as LG
from threshold_transfer.analysis import statesuff as A

TAUS = (7, 14, 28, 56)
N0, T_END = LG.N0, LG.T_END
N_TOP = {t: LG.BAND[t][-1] for t in TAUS}                 # cap: the highest training threshold
CELL = {7: 3120, 14: 2950, 28: 2740, 56: 2740}            # threshold indicator fixed in the models
NSAFE_MAIN = {7: 3114.6, 14: 2946.1, 28: 2719.9, 56: 2718.3}   # main-run N_safe; sets the SIZE search range only
SIZE_BR = {t: (float(round(NSAFE_MAIN[t] - 200)), float(N_TOP[t])) for t in TAUS}
RULE_BR = (-5.0, 0.0)
EXT_STEP = {"SIZE": 200.0, "rule": 2.0}
ARMS = ("SIZE", "CLIN", "HIST", "FULL", "VEL")
RULES = ARMS[1:]
MODEL_OF = {"CLIN": "CLIN", "HIST": "A1", "FULL": "S1F", "VEL": "VEL"}
FEATS = {**A.MODELS, "VEL": ["dN"]}                         # features of each model
NEED_GEO = {"CLIN", "FULL"}
BLOCKS = {"gate": 2_600_000, "tune": 2_400_000, "score": 2_500_000}
NSEED = {"gate": 10, "tune": 500, "score": 1000}
SHARD = 25
Q = 0.05
NW = WORKERS
OUT = RESULTS / "p_s_policy"
_M = {}


def canon(arm, x):
    """Setting rounded as stored: 4 decimals for SIZE thresholds, 6 for log10 p*."""
    return round(float(x), 4 if arm == "SIZE" else 6)


def sstr(arm, x):
    """Setting as written in shard file names."""
    return f"{canon(arm, x):.4f}" if arm == "SIZE" else f"{canon(arm, x):.6f}"


def load_models(check=True):
    """Frozen models per (tau, rule), each checked against its stored hash."""
    F5 = pickle.load(open(RESULTS / "p_s_models.pkl", "rb"))["F5"]
    FV = pickle.load(open(RESULTS / "p_s_vel_model.pkl", "rb"))
    if check:
        h5 = json.load(open(RESULTS / "p_s_fit.json"))["hash5"]
        hv = json.load(open(RESULTS / "p_s_vel_fit.json"))["hash"]
        assert A.model_hash(F5) == h5 and A.model_hash(FV) == hv, "model hash mismatch"
    M = {}
    for t in TAUS:
        for arm in RULES:
            m = MODEL_OF[arm]
            M[(t, arm)] = FV[(t, "VEL")] if arm == "VEL" else F5[(t, m)]
    return M


def _init():
    """Pool initialiser: compile the simulator kernels and load the frozen models."""
    from threshold_transfer.model.simulator import warmup_jit
    warmup_jit()
    _M.update(load_models(check=True))


def history_feats(Ns, ons, e, tau):
    """History features at decision e from the sizes Ns and treatment flags ons (same as statesuff.history)."""
    last = -1; cnt = 0
    for i in range(e):
        if ons[i]: last = i; cnt += 1
    N = np.asarray(Ns, float)
    if last >= 0:
        tsf = (e - last) * tau; never = 0.0; nadir = N[last+1:e+1].min() if e > last else N[e]
    else:
        tsf = e * tau; never = 1.0; nadir = N[:e+1].min()
    return dict(tsf=float(tsf), never_on=never, g=float(np.log(N[e] / N[e-1])),
                log_non=float(np.log1p(cnt)), nadir=float(nadir), dN=float(N[e] - N[e-1]))


def predict_one(arm, tau, feats):
    """Predicted risk of one rule's model at one decision."""
    prep, lr = _M[(tau, arm)]
    D = {k: np.array([v], float) for k, v in feats.items()}
    D["thr"] = np.array([CELL[tau]])
    X, _ = A.design(D, np.array([0]), FEATS[MODEL_OF[arm]], tau, prep=prep)
    return float(np.clip(lr.predict_proba(X)[0, 1], 1e-15, 1 - 1e-15))


def run_arm(seed, arm, tau, setting, thr_size, log_feats=False, decisions=None):
    """One run of an arm. setting is the SIZE threshold or log10 p*; thr_size is the tuned SIZE threshold,
    used by the rules at the first decision; decisions, if given, replays recorded treatment flags.
    Returns outcomes (prem, t_prog, burden, drug days, run length T, decisions forced by the cap), the
    decision rows and the simulator's generator state before and after.
    """
    sim, _, _ = P.build(seed)
    rng0 = json.dumps(sim._rng.bit_generator.state, sort_keys=True, default=str)
    prog = P.PROG_MUL * N0
    burden = 0.0; on = True; prem = 0; t_prog = np.nan
    Ntrace_max = 0.0; held = False
    Ns, ons, eps = [], [], []
    caps = 0; drug = 0; k_last = 0; N = N0
    pstar = None if arm == "SIZE" else 10.0 ** setting
    for k in range(1, T_END + 1):
        if (k - 1) % tau == 0:
            e = (k - 1) // tau
            V = float(sim.get_state()["true_total_volume"])
            p = np.nan; feats = None; cap = 0
            Ns.append(V)
            if decisions is not None:
                on = bool(decisions[e])
            elif arm == "SIZE":
                on = (V >= setting)
            elif e == 0:
                on = (V >= thr_size)
            elif V > N_TOP[tau]:
                on = True; cap = 1; caps += 1
            else:
                feats = history_feats(Ns, ons + [None], e, tau)
                feats["N"] = V
                if arm in NEED_GEO or log_feats:
                    r = LG.epoch_reads(sim)
                    for f in ("C", "R", "B", "B_R", "K", "Rg", "Rmax"): feats[f] = float(r[f])
                p = predict_one(arm, tau, feats)
                on = not (p <= pstar)
            held = not on
            ons.append(on)
            eps.append(dict(e=e, N=V, on=int(on), p=p, cap=cap, y=0, trunc=0, feats=feats))
        sim.step(u_A=P.DOSE if on else 0.0, u_B=0.0)
        if on: drug += 1
        st = sim.get_state(); N = float(st["true_total_volume"])
        burden += N; k_last = k
        if N > Ntrace_max: Ntrace_max = N
        if N >= prog:
            if np.isnan(t_prog): t_prog = float(k)
            if held: prem = 1
            break
        if N <= 0: break
    if prem == 1:
        eps[-1]["y"] = 1
    elif np.isnan(t_prog) and eps[-1]["e"] * tau + tau > T_END and N > 0:
        eps[-1]["trunc"] = 1
    rng1 = json.dumps(sim._rng.bit_generator.state, sort_keys=True, default=str)
    return dict(prem=prem, t_prog=t_prog, burden=burden, N_end=N, N_max=Ntrace_max, drug=drug,
                T=k_last, caps=caps, eps=eps, rng0=rng0, rng1=rng1)


RUN_DT = [('arm', 'U4'), ('tau', 'i2'), ('setting', 'f8'), ('s', 'i4'), ('prem', 'i1'), ('t_prog', 'f8'),
          ('burden', 'f8'), ('N_end', 'f8'), ('N_max', 'f8'), ('drug', 'i4'), ('T', 'i4'), ('caps', 'i4')]
EP_DT = [('arm', 'U4'), ('tau', 'i2'), ('s', 'i4'), ('e', 'i2'), ('N', 'f8'), ('on', 'i1'), ('p', 'f8'),
         ('cap', 'i1'), ('y', 'i1'), ('trunc', 'i1')]


def shard_job(a):
    """25 seeds of one arm and setting; a = (role, arm, tau, setting, thr_size, first seed, keep_eps)."""
    role, arm, tau, setting, thr_size, c, keep_eps = a
    runs, eps = [], []
    for s in range(c, c + SHARD):
        r = run_arm(BLOCKS[role] + s, arm, tau, setting, thr_size)
        runs.append((arm, tau, setting, s, r["prem"], -1.0 if np.isnan(r["t_prog"]) else r["t_prog"],
                     r["burden"], r["N_end"], r["N_max"], r["drug"], r["T"], r["caps"]))
        if keep_eps:
            eps += [(arm, tau, s, x["e"], x["N"], x["on"], x["p"], x["cap"], x["y"], x["trunc"]) for x in r["eps"]]
    return a, np.array(runs, dtype=RUN_DT), np.array(eps, dtype=EP_DT)


def shard_path(role, arm, tau, setting, c):
    """Path of one shard file."""
    return OUT / role / f"{arm}_{tau}_{sstr(arm, setting)}__{c:03d}.npz"


def save(path, runs, eps):
    """Write one shard file atomically (temporary file, then rename)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name[:-4] + ".tmp.npz")
    np.savez(tmp, runs=runs, eps=eps); os.replace(tmp, path)


def log_wall(stage, wall, n):
    """Append the wall time of one invocation to wall_log.jsonl."""
    with open(OUT / "wall_log.jsonl", "a") as f:
        f.write(json.dumps(dict(t=time.strftime("%Y-%m-%d %H:%M:%S"), stage=stage, wall_s=wall, shards=n)) + "\n")


def _g_ref(a):
    """Gate reference run of SIZE through construction.run."""
    tau, thr, s = a
    r = P.run(BLOCKS["gate"] + s, "AT", tau=tau, thr=float(thr), N0=N0, horizon=T_END)
    return (r["prem"], -1.0 if np.isnan(r["t_prog"]) else r["t_prog"], r["burden"], float(r["N_end"]), r["N_max"])


def _g_run(a):
    """Gate run of one arm, with features logged and the wall time recorded."""
    arm, tau, setting, s = a
    t0 = time.time()
    r = run_arm(BLOCKS["gate"] + s, arm, tau, setting, CELL[tau], log_feats=(arm != "SIZE"))
    r["wall"] = time.time() - t0
    return a, r


def _g_replay(a):
    """Gate replay of a run with its recorded treatment decisions."""
    (arm, tau, setting, s), dec = a
    return run_arm(BLOCKS["gate"] + s, arm, tau, setting, CELL[tau], decisions=dec)


def gate():
    """Check the rule code on the gate block, 10 seeds per arm and interval: (1) SIZE equals construction.run;
    (2) replaying the decisions reproduces each run and its generator state; (3) features and predicted risks
    equal statesuff.history and the models; (4) labels match outcomes. Also projects run times. Writes gate.json.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    size_args = [(t, CELL[t], s) for t in TAUS for s in range(NSEED["gate"])]
    run_args = [("SIZE", t, float(CELL[t]), s) for t in TAUS for s in range(NSEED["gate"])] + \
               [(arm, t, -2.0, s) for arm in RULES for t in TAUS for s in range(NSEED["gate"])]
    with Pool(NW, initializer=_init) as pool:
        ref = pool.map(_g_ref, size_args, chunksize=1)
        t0 = time.time(); runs = pool.map(_g_run, run_args, chunksize=1); w_all = time.time() - t0
        rep = pool.map(_g_replay, [(a, [x["on"] for x in r["eps"]]) for a, r in runs], chunksize=1)
    res = {}
    F = ("prem", "t_prog", "burden", "N_end", "N_max")
    size_runs = [r for a, r in runs if a[0] == "SIZE"]
    mism = sum(float(v) != float(-1.0 if (f == "t_prog" and np.isnan(r[f])) else r[f])
               for rf, r in zip(ref, size_runs) for f, v in zip(F, rf))
    res["item1_size_vs_p_t_pilot_mismatches"] = int(mism)
    rng_eq = sum(r["rng1"] == rp["rng1"] and all(
        (np.isnan(r[f]) and np.isnan(rp[f])) or r[f] == rp[f] for f in F + ("drug", "T"))
        for (a, r), rp in zip(runs, rep))
    moved = sum(r["rng1"] != r["rng0"] for a, r in runs)
    res.update(item2_replay_identical=int(rng_eq), item2_rng_advanced=int(moved), n_runs=len(runs))
    # item 3: features equal statesuff.history, and p_hat recomputed equals the logged value
    _M.update(load_models(check=True))
    bad_feat = 0; bad_p = 0; nchk = 0
    for a, r in runs:
        arm, tau, setting, s = a
        if arm == "SIZE": continue
        eps = r["eps"]; n = len(eps)
        E = np.zeros(n, dtype=[('tau', 'i2'), ('thr', 'i2'), ('s', 'i2'), ('e', 'i2'), ('N', 'i4'), ('on', 'i1'),
                               ('C', 'f8'), ('R', 'f8'), ('B', 'f8'), ('B_R', 'f8'), ('K', 'f8'), ('Rg', 'f8'), ('Rmax', 'f8')])
        E['tau'] = tau; E['thr'] = CELL[tau]; E['s'] = s
        E['e'] = [x["e"] for x in eps]; E['N'] = [int(x["N"]) for x in eps]; E['on'] = [x["on"] for x in eps]
        H = A.history(E)
        for i, x in enumerate(eps):
            if x["feats"] is None: continue
            nchk += 1
            for f in ("tsf", "never_on", "g", "log_non", "nadir"):
                bad_feat += int(H[f][i] != x["feats"][f])
            bad_feat += int(x["feats"]["dN"] != eps[i]["N"] - eps[i-1]["N"])
            fe = {f: float(H[f][i]) for f in ("tsf", "never_on", "g", "log_non", "nadir")}
            fe.update(N=x["N"], dN=x["N"] - eps[i-1]["N"],
                      **{f: x["feats"][f] for f in ("C", "R", "B", "B_R", "K", "Rg", "Rmax")})
            bad_p += int(predict_one(arm, tau, fe) != x["p"])
    res.update(item3_epochs_checked=nchk, item3_feature_mismatches=bad_feat, item3_p_mismatches=bad_p)
    lab = 0
    for a, r in runs:
        eps = r["eps"]
        lab += int(sum(x["y"] for x in eps if not x["on"]) != r["prem"]) + int(sum(x["y"] for x in eps if x["on"]) != 0)
    res["item4_label_violations"] = lab
    # cost pilot: single-run wall time inside workers, microseconds per simulated day, by arm
    per = {}
    for arm in ARMS:
        rr = [r for a, r in runs if a[0] == arm]
        per[arm] = 1e6 * sum(r["wall"] for r in rr) / sum(r["T"] for r in rr)
    res["us_per_day_single_worker"] = per
    res["gate_wall_s"] = w_all
    days = sum(r["T"] for a, r in runs)
    res["us_per_day_throughput_7w"] = 1e6 * w_all / days
    ratio = {arm: per[arm] / per["SIZE"] for arm in ARMS}
    res["arm_cost_ratio_vs_SIZE"] = ratio
    cen = 354.0                               # central run length in days, for the projections
    tune_runs = {"SIZE": 4 * 8 * 500, **{a: 4 * 8 * 500 for a in RULES}}
    base_us = res["us_per_day_throughput_7w"] / (sum(ratio[a] for a in ARMS) / len(ARMS))
    res["tune_projection_min_central_pilot"] = sum(n * cen * base_us * ratio[a] for a, n in tune_runs.items()) / 1e6 / 60
    res["score_projection_min_central_pilot"] = sum(1000 * 4 * cen * base_us * ratio[a] for a in ARMS) / 1e6 / 60
    res["score_projection_min_worst_pilot"] = sum(1000 * 4 * 750 * base_us * ratio[a] for a in ARMS) / 1e6 / 60
    ok = (mism == 0 and rng_eq == len(runs) and moved == len(runs) and bad_feat == 0 and bad_p == 0 and lab == 0)
    res["GATE"] = "PASS" if ok else "FAIL"
    res["TUNE_LAUNCH"] = bool(ok and res["tune_projection_min_central_pilot"] <= 360.0)
    json.dump(res, open(OUT / "gate.json", "w"), indent=1)
    print(json.dumps(res, indent=1), flush=True)


def p_prem_of(role, arm, tau, setting):
    """Tuning p_prem of one arm and setting, or None if its files are not all there."""
    fs = [shard_path(role, arm, tau, setting, c) for c in range(0, NSEED[role], SHARD)]
    if not all(f.exists() for f in fs): return None
    return float(np.concatenate([np.load(f)["runs"]["prem"] for f in fs]).mean())


def search(arm, tau):
    """Next setting to evaluate, ('need', setting), or the result, ('done', info).

    Starts at the lower end of the range (stepping lower up to twice if it is unsafe), checks the upper end,
    then bisects 6 times, keeping the lower end safe.
    """
    kind = "SIZE" if arm == "SIZE" else "rule"
    lo, hi = SIZE_BR[tau] if arm == "SIZE" else RULE_BR
    ev = []
    get = lambda x: p_prem_of("tune", arm, tau, canon(arm, x))
    ext = 0; extended = False
    while True:
        p = get(lo)
        if p is None: return ("need", canon(arm, lo))
        ev.append((canon(arm, lo), p))
        if p <= Q: break
        if ext == 2: return ("done", dict(status="SAFETY-INFEASIBLE", selected=None, evals=ev))
        ext += 1; extended = True; hi = lo; lo = canon(arm, lo - EXT_STEP[kind])
    if not extended:
        p = get(hi)
        if p is None: return ("need", canon(arm, hi))
        ev.append((canon(arm, hi), p))
        if p <= Q: return ("done", dict(status="CAP-LIMITED", selected=canon(arm, hi), p_prem=p, evals=ev))
    for _ in range(6):
        mid = canon(arm, (lo + hi) / 2)
        p = get(mid)
        if p is None: return ("need", mid)
        ev.append((mid, p))
        if p <= Q: lo = mid
        else: hi = mid
    return ("done", dict(status="OK", selected=canon(arm, lo), p_prem=get(lo), evals=ev))


def wilson(k, n, z=1.959964):
    """Wilson score interval (Wilson 1927) for k events in n runs."""
    ph = k / n; d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d; h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return [float(c - h), float(c + h)]


def tune(budget):
    """Run tuning files until every search is done or budget seconds pass. Needs a passed gate."""
    g = json.load(open(OUT / "gate.json"))
    if not g["TUNE_LAUNCH"]: print("gate.json: TUNE_LAUNCH false; refusing"); sys.exit(3)
    t0 = time.time(); done_n = 0
    pool = Pool(NW, initializer=_init)
    try:
        while True:
            state = {}; jobs = []
            for tau in TAUS:
                st = search("SIZE", tau); state[("SIZE", tau)] = st
                if st[0] == "need":
                    jobs += [("tune", "SIZE", tau, st[1], None, c, False) for c in range(0, 500, SHARD)
                             if not shard_path("tune", "SIZE", tau, st[1], c).exists()]
                    continue
                thr_size = st[1]["selected"]
                for arm in RULES:
                    s2 = search(arm, tau); state[(arm, tau)] = s2
                    if s2[0] == "need":
                        jobs += [("tune", arm, tau, s2[1], thr_size, c, False) for c in range(0, 500, SHARD)
                                 if not shard_path("tune", arm, tau, s2[1], c).exists()]
            if not jobs:
                if all(v[0] == "done" for v in state.values()) and len(state) == 20:
                    finish_tune(state); break
                print("no jobs but searches incomplete?", flush=True); break
            over = False
            for a, runs, eps in pool.imap_unordered(shard_job, jobs):
                save(shard_path(*a[:4], a[5]), runs, eps); done_n += 1
                if time.time() - t0 > budget: over = True; break
            if over: break
    finally:
        pool.terminate(); pool.join()
    w = time.time() - t0; log_wall("tune", w, done_n)
    print(f"tune invocation wall {w:.1f}s  shards {done_n}", flush=True)


def finish_tune(state):
    """Write tune_selected.json: setting, tuned p_prem and Wilson interval per arm and interval, and UNMATCHED."""
    sel = {}
    for (arm, tau), (_, info) in sorted(state.items()):
        k = None
        if info["selected"] is not None:
            k = int(round(info["p_prem"] * 500))
        sel[f"{arm}|{tau}"] = dict(status=info["status"], selected=info["selected"],
                                   p_prem_tune=info.get("p_prem"), wilson=wilson(k, 500) if k is not None else None,
                                   evals=info["evals"])
    for tau in TAUS:
        a, b = sel[f"CLIN|{tau}"], sel[f"SIZE|{tau}"]
        sel[f"CLIN|{tau}"]["UNMATCHED"] = bool(a["p_prem_tune"] is not None and abs(a["p_prem_tune"] - b["p_prem_tune"]) > 0.010)
    blob = json.dumps(sel, sort_keys=True)
    json.dump(dict(t=time.strftime("%Y-%m-%d %H:%M:%S"), selected=sel,
                   sha256=hashlib.sha256(blob.encode()).hexdigest()), open(OUT / "tune_selected.json", "w"), indent=1)
    print("TUNE COMPLETE", flush=True)
    for k, v in sel.items():
        print(f"  {k:8s} {v['status']:18s} sel={v['selected']}  p_tune={v['p_prem_tune']}  wilson={v['wilson']}"
              + ("  UNMATCHED" if v.get("UNMATCHED") else ""), flush=True)


def selected():
    """Selected settings, checked against their stored SHA-256."""
    d = json.load(open(OUT / "tune_selected.json"))
    assert hashlib.sha256(json.dumps(d["selected"], sort_keys=True).encode()).hexdigest() == d["sha256"]
    return d["selected"]


def score_jobs():
    """All scoring files for the selected settings."""
    sel = selected(); jobs = []
    for tau in TAUS:
        thr_size = sel[f"SIZE|{tau}"]["selected"]
        for arm in ARMS:
            v = sel[f"{arm}|{tau}"]
            if v["selected"] is None: continue
            jobs += [("score", arm, tau, v["selected"], thr_size, c, True) for c in range(0, 1000, SHARD)]
    return jobs


def score(budget):
    """Run scoring files until all exist or budget seconds pass."""
    jobs = [j for j in score_jobs() if not shard_path(*j[:4], j[5]).exists()]
    print(f"score todo shards {len(jobs)}", flush=True)
    t0 = time.time(); n = 0
    if jobs:
        pool = Pool(NW, initializer=_init)
        try:
            for a, runs, eps in pool.imap_unordered(shard_job, jobs):
                save(shard_path(*a[:4], a[5]), runs, eps); n += 1
                if time.time() - t0 > budget: break
        finally:
            pool.terminate(); pool.join()
    w = time.time() - t0; log_wall("score", w, n)
    left = sum(not shard_path(*j[:4], j[5]).exists() for j in score_jobs())
    print(f"score invocation wall {w:.1f}s  shards {n}  remaining {left}", flush=True)
    print("COMPLETE" if left == 0 else "INCOMPLETE", flush=True)


def _ci(b):
    """95% percentile interval of a bootstrap sample."""
    return [float(x) for x in np.percentile(b, [2.5, 97.5])]


def reading(d, safe_fail, infeasible):
    """Reading of a rule against SIZE, by the rules of the analysis plan."""
    if infeasible: return "NOT ASSESSABLE (SAFETY-INFEASIBLE)"
    if safe_fail: return "NOT ASSESSABLE (SAFETY FAILURE)"
    dp = d["dp"]
    if dp["ci"][1] > 0.02:
        return "NOT ASSESSABLE (" + ("bound fails on width" if dp["mean"] <= 0.02 else "risk not matched") + ")"
    eb = d["df"]["ci"][1] < 0; ew = d["df"]["ci"][0] > 0
    tb = d["drmst"]["ci"][0] > 0; tw = d["drmst"]["ci"][1] < 0
    if (eb and tw) or (tb and ew): return "TRADE-OFF"
    if eb and tb: return "GAIN (exposure + time)"
    if eb: return "GAIN (exposure)"
    if tb: return "GAIN (time)"
    return "NO GAIN" + (" (worse)" if (ew or tw) else "")


def analyse():
    """Paired differences of each rule against SIZE, per interval and pooled (equal-weight mean), readings,
    the growth-rate rule's recovery ratio rho, and observed/expected events on each rule's own holds.
    Writes p_s_policy_results.json.
    """
    sel = selected()
    R = {}
    for j in score_jobs():
        z = np.load(shard_path(*j[:4], j[5]))
        R.setdefault((j[1], j[2]), []).append(z)
    runs = {k: np.sort(np.concatenate([z["runs"] for z in v]), order="s") for k, v in R.items()}
    eps = {k: np.concatenate([z["eps"] for z in v]) for k, v in R.items()}
    Wb = A.boot_weights(1000)
    out = dict(selected=sel, per=dict(), pooled=dict(), diagnostics=dict())
    vec = {}
    for (arm, tau), r in runs.items():
        assert np.array_equal(r["s"], np.arange(1000))
        T = r["T"].astype(float); f = r["drug"] / T
        rm = np.where(r["t_prog"] > 0, np.minimum(r["t_prog"], T_END), T_END)
        vec[(arm, tau)] = dict(prem=r["prem"].astype(float), f=f, D=r["drug"].astype(float), rmst=rm)
        k = int(r["prem"].sum()); w = wilson(k, 1000)
        ph = k / 1000
        out["per"][f"{arm}|{tau}"] = dict(
            p_prem=ph, wilson=w, SAFETY_FAILURE=bool(w[0] > Q), NOMINAL_EXCEEDANCE=bool(ph > Q and w[1] > Q),
            on_drug_fraction=float(f.mean()), drug_days=float(r["drug"].mean()), rmst750=float(rm.mean()),
            censored=int((r["t_prog"] < 0).sum()), cap_forced_mean=float(r["caps"].mean()))
    def paired(a, b_, key):
        d = a[key] - b_[key]; bs = (Wb @ d) / 1000
        rel = (Wb @ d) / (Wb @ b_[key])
        return dict(mean=float(d.mean()), ci=_ci(bs), rel=float(d.mean() / b_[key].mean()), rel_ci=_ci(rel))
    for arm in RULES:
        pool_vec = {k: [] for k in ("prem", "f", "D", "rmst")}
        pool_sz = {k: [] for k in ("prem", "f", "D", "rmst")}
        for tau in TAUS:
            key = f"{arm}|{tau}"
            if (arm, tau) not in vec:
                out["per"][key] = dict(reading="NOT ASSESSABLE (SAFETY-INFEASIBLE)"); continue
            a, s_ = vec[(arm, tau)], vec[("SIZE", tau)]
            d = dict(dp=paired(a, s_, "prem"), df=paired(a, s_, "f"), dD=paired(a, s_, "D"), drmst=paired(a, s_, "rmst"))
            d["SAFETY_SHIFT"] = bool(d["dp"]["ci"][0] > 0)
            d["reading"] = reading(d, out["per"][key]["SAFETY_FAILURE"], sel[key]["selected"] is None)
            out["per"][key].update(d)
            for k in pool_vec: pool_vec[k].append(a[k]); pool_sz[k].append(s_[k])
        if len(pool_vec["f"]) == 4:
            pa = {k: np.mean(v, axis=0) for k, v in pool_vec.items()}
            ps = {k: np.mean(v, axis=0) for k, v in pool_sz.items()}
            d = dict(dp=paired(pa, ps, "prem"), df=paired(pa, ps, "f"), dD=paired(pa, ps, "D"), drmst=paired(pa, ps, "rmst"))
            d["reading"] = reading(d, False, False)
            sup = [t for t in TAUS if out["per"][f"{arm}|{t}"].get("reading", "").split(" (")[0] == d["reading"].split(" (")[0]
                   and np.sign(out["per"][f"{arm}|{t}"]["df"]["mean"]) == np.sign(d["df"]["mean"])]
            d["supporting_tau"] = sup
            d["headline_stands"] = bool(len(sup) >= 2 and d["reading"].startswith(("GAIN", "NO GAIN", "TRADE-OFF")))
            out["pooled"][arm] = d
    # share of the history rule's exposure gain recovered by the growth-rate rule
    rec = {}
    for tau in list(TAUS) + ["pool"]:
        if tau == "pool":
            if "VEL" not in out["pooled"] or "HIST" not in out["pooled"]: continue
            h, v = out["pooled"]["HIST"], out["pooled"]["VEL"]
            dv = np.mean([vec[("VEL", t)]["f"] - vec[("SIZE", t)]["f"] for t in TAUS], axis=0)
            dh = np.mean([vec[("HIST", t)]["f"] - vec[("SIZE", t)]["f"] for t in TAUS], axis=0)
        else:
            if ("VEL", tau) not in vec or ("HIST", tau) not in vec: continue
            h, v = out["per"][f"HIST|{tau}"], out["per"][f"VEL|{tau}"]
            dv = vec[("VEL", tau)]["f"] - vec[("SIZE", tau)]["f"]; dh = vec[("HIST", tau)]["f"] - vec[("SIZE", tau)]["f"]
        rho = float(dv.mean() / dh.mean()); b = (Wb @ dv) / (Wb @ dh)
        assessable = h["reading"].startswith("GAIN (exposure") and not v["reading"].startswith("NOT ASSESSABLE")
        rec[str(tau)] = dict(rho=rho, ci=_ci(b), assessable=bool(assessable),
                             reading=("size + growth rate carries it (rho >= 0.80)" if rho >= 0.8 else "richer history matters")
                             if assessable else "NOT ASSESSABLE")
    out["vel_recovery"] = rec
    # calibration of each rule's model on the holds it chose
    M = load_models(check=True)
    for arm in RULES:
        for tau in TAUS:
            if (arm, tau) not in eps: continue
            E = eps[(arm, tau)]; m = (E["on"] == 0) & ~np.isnan(E["p"]) & (E["trunc"] == 0)
            t = M[(tau, arm)][0]["spl"].bsplines_[0].t; lo, hi = t[3], t[-4]
            out["diagnostics"][f"{arm}|{tau}"] = dict(
                hold_epochs=int(m.sum()), events=int(E["y"][m].sum()), expected=float(E["p"][m].sum()),
                O_over_E=float(E["y"][m].sum() / max(E["p"][m].sum(), 1e-12)),
                frac_N_outside_knots=float(((E["N"][m] < lo) | (E["N"][m] > hi)).mean()) if m.any() else None)
    json.dump(A._json(out), open(RESULTS / "p_s_policy_results.json", "w"), indent=1)
    print("arm|tau  p_prem [Wilson]  f  RMST  | Δf [CI]  ΔRMST [CI]  Δp [CI]  reading")
    for tau in TAUS:
        s_ = out["per"][f"SIZE|{tau}"]
        print(f"SIZE|{tau:<3d} {s_['p_prem']:.3f} [{s_['wilson'][0]:.3f},{s_['wilson'][1]:.3f}]  f={s_['on_drug_fraction']:.3f}  RMST={s_['rmst750']:.1f}")
        for arm in RULES:
            v = out["per"].get(f"{arm}|{tau}")
            if "df" not in v: print(f"  {arm}|{tau}: {v.get('reading')}"); continue
            print(f"  {arm:4s} {v['p_prem']:.3f} [{v['wilson'][0]:.3f},{v['wilson'][1]:.3f}] f={v['on_drug_fraction']:.3f} "
                  f"RMST={v['rmst750']:.1f} | Δf {v['df']['mean']:+.4f} [{v['df']['ci'][0]:+.4f},{v['df']['ci'][1]:+.4f}] "
                  f"({100*v['df']['rel']:+.1f}%) ΔRMST {v['drmst']['mean']:+.1f} [{v['drmst']['ci'][0]:+.1f},{v['drmst']['ci'][1]:+.1f}] "
                  f"Δp {v['dp']['mean']:+.3f} [{v['dp']['ci'][0]:+.3f},{v['dp']['ci'][1]:+.3f}]  {v['reading']}")
    for arm, d in out["pooled"].items():
        print(f"POOLED {arm}: Δf {d['df']['mean']:+.4f} [{d['df']['ci'][0]:+.4f},{d['df']['ci'][1]:+.4f}] rel {100*d['df']['rel']:+.1f}% "
              f"[{100*d['df']['rel_ci'][0]:+.1f},{100*d['df']['rel_ci'][1]:+.1f}]  ΔRMST {d['drmst']['mean']:+.1f}  Δp {d['dp']['mean']:+.3f} "
              f"[{d['dp']['ci'][0]:+.3f},{d['dp']['ci'][1]:+.3f}]  {d['reading']}  support {d['supporting_tau']}  stands={d['headline_stands']}")
    print("VEL recovery:", {k: (round(v["rho"], 3), v["reading"]) for k, v in rec.items()})
    print("written p_s_policy_results.json")
