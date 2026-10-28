"""Number ledger: every number in the paper, its source and how it is obtained (`python run.py ledger`).

kind "recomputed": computed here from data/results or read from a constant in the code that ran.
kind "filing": a value that exists only as text in a dated record under config/ (for example a bootstrap
interval whose random state was not saved); the text is checked to be present there.
Where a recomputed value was also recorded, the two are compared.
"""
import json
import numpy as np
from threshold_transfer.paths import ROOT, SRC, RESULTS
RES = RESULTS
Q = 0.05
L = {}                      # id -> entry


def isafe(x, p):
    """N_safe: largest x with p_prem <= Q, linear in p_prem between grid points (p made non-decreasing first)."""
    p = np.maximum.accumulate(np.asarray(p, float)); x = np.asarray(x, float)
    if p[0] > Q: return -np.inf
    if p[-1] <= Q: return np.inf
    i = int(np.argmax(p > Q)); p0, p1 = p[i-1], p[i]
    return x[i-1] if p1 <= p0 else x[i-1] + (Q - p0) * (x[i] - x[i-1]) / (p1 - p0)


def fmt(v, f):
    """Render a value with a format spec; strings pass through."""
    if isinstance(v, str): return v
    return format(v, f)


def R(id_, value, f, artifact, method, filed=None):
    """Add a recomputed entry; filed = (record file, exact text, recorded value, tolerance) or None."""
    L[id_] = dict(kind="recomputed", value=float(value) if not isinstance(value, str) else value, fmt=f,
                  text=fmt(value, f), artifact=artifact, method=method, filed=filed)


def F(id_, text, file, substring, note=""):
    """Add an entry quoted from a dated record: `substring` must appear in `file`."""
    L[id_] = dict(kind="filing", value=None, fmt=None, text=text, artifact=file,
                  method=f"verbatim in {file}" + (f" ({note})" if note else ""), substring=substring)


TAUS6 = (7, 14, 28, 42, 56, 84); TAUS4 = (7, 14, 28, 56)
RM = "config/filings/RESULTS_2026-08-24_threshold_transfer_main.md"
RL = "config/filings/research_log_excerpts.md"

# configuration
from threshold_transfer.threshold import grid as G, construction as P
R("cfg.N0", G.N0, ".0f", "src/threshold_transfer/model/ode.py", "N0 constant")
R("cfg.L", P.L, "d", "src/threshold_transfer/model/ode.py", "lattice edge L")
R("cfg.nmax", P.N_MAX, ".0f", "src/threshold_transfer/model/ode.py", "N_MAX = L*L, threaded explicitly")
R("cfg.prog", P.PROG_MUL * G.N0, ".1f", "src/threshold_transfer/model/ode.py", "PROG_MUL * N0")
R("cfg.rfrac", 100 * P.RFRAC, ".0f", "src/threshold_transfer/threshold/construction.py", "RFRAC, percent")
R("cfg.ric", P.R_IC, "d", "src/threshold_transfer/threshold/construction.py", "R_IC disk radius")
R("cfg.tend", G.T_END, "d", "src/threshold_transfer/threshold/grid.py", "T_END horizon, days")

# main run
z = np.load(RES / "p_t_main_raw.npz"); rec = z["rec"]; NS = dict(zip(z["taus"].tolist(), z["nstar"]))
R("main.nruns", len(rec), "d", "data/results/p_t_main_raw.npz", "len(rec)")
R("main.ncells", len(set(zip(rec["pol"], rec["tau"], rec["m"]))), "d", "data/results/p_t_main_raw.npz", "distinct (pol,tau,m)")
nstar_json = {int(k): v for k, v in json.load(open(RES / "p_t_nstar.json")).items()}
main = {}
for t in TAUS6:
    s = rec[(rec["pol"] == "AT") & (rec["tau"] == t)]; ms = sorted(set(s["m"]))
    main[t] = isafe([m * NS[t] for m in ms], [s[s["m"] == m]["prem"].mean() for m in ms])
    R(f"nstar.{t}", nstar_json[t], ".1f", "data/results/p_t_nstar.json", "guard-exact N*(tau), n_max = L^2")
    R(f"main.nsafe.{t}", main[t], ".1f", "data/results/p_t_main_raw.npz",
      "interp_safe (linear in p between bracketing m-cells, running max), thr = m*N*",
      filed=(RM, f"| {t} | {main[t]:.1f} |", round(main[t], 1), 0.05))
    step = (ms[1] - ms[0]) * NS[t]
    R(f"main.step.{t}", step, ".1f", "data/results/p_t_main_raw.npz", "0.05 m-step in cells")
CI = {7: "[3110.5, 3119.9]", 14: "[2944.8, 2947.6]", 28: "[2702.3, 2748.0]", 42: "[2662.4, 2693.2]",
      56: "[2666.4, 2739.8]", 84: "[2499.5, 2630.8]"}
for t in TAUS6:
    F(f"main.nsafe_ci.{t}", CI[t], RM, f"| {t} | {main[t]:.1f} | {CI[t]} |", "R3.5 seed-clustered bootstrap")
ODE = {7: 3454.0, 14: 3299.4, 28: 3001.0, 42: 2787.7, 56: 2831.6, 84: 2566.2}   # ODE safe thresholds as filed; recomputed by tests/regression.py --full
DER = {7: "+0.84%", 14: "+1.36%", 28: "+2.30%", 42: "+5.78%", 56: "+20.00%", 84: "+37.02%"}
TRF = {7: "−9.83%", 14: "−10.71%", 28: "−9.37%", 42: "−4.00%", 56: "−4.00%", 84: "−0.15%"}
for t in TAUS6:
    F(f"ode.bnd.{t}", f"{ODE[t]:.1f}", RM, f"| {t} | {nstar_json[t]:.1f} | {ODE[t]:.1f} |", "R0, bisection of the deterministic ODE arm")
    d = 100 * (ODE[t] / nstar_json[t] - 1); tr = 100 * (main[t] / ODE[t] - 1); tot = 100 * (main[t] / nstar_json[t] - 1)
    R(f"lead.der.{t}", d, "+.2f", "p_t_nstar.json + R0 ODE boundary", "ODE boundary / N* − 1",
      filed=(RM, f"| {t} | {DER[t]} |", float(DER[t].replace('%', '')), 0.011))
    R(f"lead.trf.{t}", tr, "+.2f", "p_t_main_raw.npz + R0 ODE boundary", "N_safe(ABM) / ODE boundary − 1",
      filed=(RM, f"**{TRF[t]}**", float(TRF[t].replace('%', '').replace('−', '-')), 0.011))
    R(f"lead.tot.{t}", tot, "+.2f", "p_t_main_raw.npz + p_t_nstar.json", "N_safe(ABM) / N* − 1")
F("lead.transit_abm", "60.3", RM, "60.3 d in the ABM", "R1.2 untreated transit N0 -> 1.2 N0")
F("lead.transit_ode", "24.98", RM, "24.98 d in the ODE", "R1.2")
F("lead.ratio", "2.41", RM, "2.41x slower regrowth", "R1.2")
F("lead.decline_emp", "17.73", RM, "declines **17.73%**", "R3.4")
F("lead.decline_an", "45.32", RM, "analytic **45.32%**", "R3.4")
F("lead.ctmcnom7", "0.200", RM, "CTMC-nom = **0.200**", "R3.4, tau = 7, m = 1.00")
F("lead.abm7m1", "0.976", RM, "0.976 / 0.974 / 0.480", "R1.1")
F("lead.rev", "21 of 64", RM, "below CTMC-nom in 21 of 64 cells", "R3.1")

# recalibration
zr = np.load(RES / "p_t_recal_sweep.npz"); rr = zr["rec"]; NSR = dict(zip(zr["taus"].tolist(), zr["nstar_recal"]))
R("recal.nruns", len(rr), "d", "data/results/p_t_recal_sweep.npz", "len(rec)")
recal = {}
FILED_RECAL = {7: 3156.5, 14: 2994.1, 28: 2723.5, 42: 2645.8, 56: 2706.9, 84: 2581.5}
FILED_MED = {7: -11.87, 14: -14.90, 28: -19.10, 42: -17.60, 56: -11.09, 84: -5.94}
NQ05 = {7: 3553.7, 14: 3499.4, 28: 3333.3, 42: 3166.5, 56: 3000.0, 84: 2677.8}
FILED_Q05 = {7: -11.18, 14: -14.44, 28: -18.30, 42: -16.44, 56: -9.77, 84: -3.60}
for t in TAUS6:
    s = rr[rr["tau"] == t]; ms = sorted(set(s["m"]))
    recal[t] = isafe([m * NSR[t] for m in ms], [s[s["m"] == m]["prem"].mean() for m in ms])
    R(f"recal.nstar.{t}", NSR[t], ".1f", "data/results/p_t_recal_sweep.npz", "N*_recal (median transit-map inversion)")
    R(f"recal.nsafe.{t}", recal[t], ".1f", "data/results/p_t_recal_sweep.npz", "interp_safe on the block-1800000 grid",
      filed=(RL, f"{FILED_RECAL[t]}", FILED_RECAL[t], 0.05))
    R(f"recal.dmed.{t}", 100 * (recal[t] / NSR[t] - 1), "+.2f", "p_t_recal_sweep.npz", "N_safe(recal block) / N*_recal − 1",
      filed=(RL, f"{abs(FILED_MED[t]):.2f}", FILED_MED[t], 0.011))
    F(f"recal.nq05.{t}", f"{NQ05[t]:.1f}", RL, f"{NQ05[t]:.1f}", "5th-percentile transit inversion, RECAL row")
    R(f"recal.dq05.{t}", 100 * (recal[t] / NQ05[t] - 1), "+.2f", "p_t_recal_sweep.npz + filed N_q05", "N_safe / N_q05 − 1",
      filed=(RL, f"{abs(FILED_Q05[t]):.2f}", FILED_Q05[t], 0.011))
zm = np.load(RES / "p_t_recal_map.npz")
for n, row in zip(zm["map_n"], zm["transit"]):
    if n in (2200, 2600, 3000, 3400):
        R(f"map.med.{n}", np.median(row), ".1f", "data/results/p_t_recal_map.npz", f"median drug-free transit from N_start = {n}")
R("map.R", zm["R"], "d", "data/results/p_t_recal_map.npz", "seeds per N_start")
# refined CTMC and the split
odeb = {int(k): v for k, v in json.load(open(RES / "p_t_recal_odeb.json")).items()}
cells = {**json.load(open(RES / "p_t_recal_ctmc_fine.json"))["cells"], **json.load(open(RES / "p_t_recal_ctmc_fine2.json"))["cells"]}
FILED_CT = {7: 3513.2, 14: 3453.2, 28: 3321.6, 42: 3233.9, 56: 3078.7, 84: 2998.6}
FILED_SH = {7: 89.4, 14: 94.2, 28: 85.0, 42: 78.4, 56: 45.9, 84: 99.4}
ctmc = {}
for t in TAUS6:
    ks = sorted([k for k in cells if k.startswith(f"{t}|")], key=lambda k: float(k.split("|")[1]))
    ctmc[t] = isafe([cells[k]["thr"] for k in ks], [np.mean(cells[k]["prem"]) for k in ks])
    R(f"odeR.{t}", odeb[t], ".1f", "data/results/p_t_recal_odeb.json", "ODE boundary under the recalibrated law")
    R(f"ctmc.{t}", ctmc[t], ".1f", "data/results/p_t_recal_ctmc_fine{,2}.json", "interp_safe on the 0.0025-m CTMC grid",
      filed=(RL, f"{FILED_CT[t]}", FILED_CT[t], 0.05))
    R(f"share.{t}", 100 * (ctmc[t] - recal[t]) / (odeb[t] - recal[t]), ".1f", "odeb + ctmc_fine + recal_sweep",
      "(CTMC − ABM)/(ODE − ABM), recalibrated law", filed=(RL, f"{FILED_SH[t]}", FILED_SH[t], 0.051))
BANDS = {7: "[10.9,11.4]", 14: "[6.1,6.9]", 28: "[12.9,15.8]", 56: "[47.1,58.3]"}
for t in TAUS4:
    F(f"band.{t}", BANDS[t], RL, BANDS[t], "variance-share jitter band, CHECK A")
F("recal.cv", "4.49", RL, "median **CV 4.49%**", "transit spread")

# N_safe on the finer grids (correction of 2026-09-24)
from threshold_transfer.analysis import statesuff as A
r1 = np.concatenate([A.load("train")[0], A.load("test")[0]])
s1 = {}
for t in TAUS4:
    th = A.BAND[t]; s1[t] = isafe(th, [r1["prem"][(r1["tau"] == t) & (r1["thr"] == x)].mean() for x in th])
    R(f"corr.s1.{t}", s1[t], ".1f", "data/results/p_s_cells/{train,test} (1,000 seeds/cell)", "interp_safe on the Step 1 grid")
    R(f"corr.s1step.{t}", A.BAND[t][1] - A.BAND[t][0], "d", "src/threshold_transfer/experiments/statesuff_log.py BAND", "Step 1 grid step, cells")
T = json.load(open(RES / "p_s_policy" / "tune_selected.json"))["selected"]
for t in TAUS4:
    ev = T[f"SIZE|{t}"]["evals"]; lo = T[f"SIZE|{t}"]["selected"]; hi = min(x for x, p in ev if p > Q and x > lo)
    R(f"corr.tlo.{t}", lo, ".1f", "data/results/p_s_policy/tune_selected.json", "TUNE SIZE bisection: largest evaluated thr with p_prem <= 0.05")
    R(f"corr.thi.{t}", hi, ".1f", "data/results/p_s_policy/tune_selected.json", "smallest evaluated thr above it with p_prem > 0.05")
    R(f"corr.trf.{t}", 100 * (s1[t] / ODE[t] - 1), "+.2f", "Step 1 grid N_safe + R0 ODE boundary", "corrected transfer term")
    R(f"corr.share.{t}", 100 * (ctmc[t] - s1[t]) / (odeb[t] - s1[t]), ".1f", "Step 1 grid N_safe + odeb + ctmc", "spatial share with Step 1 ABM value")
    R(f"corr.dq05.{t}", 100 * (s1[t] / NQ05[t] - 1), "+.2f", "Step 1 grid N_safe + filed N_q05", "q05 discrepancy with Step 1 ABM value")
    R(f"corr.dnsafe.{t}", s1[t] - main[t], "+.1f", "Step 1 grid vs main run", "shift in N_safe, cells")

# prediction study
S = json.load(open(RES / "p_s_results.json")); P1 = S["primary"]
AB = json.load(open(RES / "p_s_absolute_scale.json"))
RD = json.load(open(RES / "p_s_rdiag.json"))
VP = json.load(open(RES / "p_s_vel_prediction.json"))
FJ = json.load(open(RES / "p_s_fit.json"))
R("s1.nruns", 2 * 20 * 500, "d", "data/results/p_s_cells", "20 cells x 500 seeds x 2 blocks (asserted at load)")
for t in TAUS4:
    v = P1["per_tau"][str(t)]
    R(f"s1.hold.{t}", v["hold_epochs"], "d", "data/results/p_s_results.json", "TEST hold epochs")
    R(f"s1.ev.{t}", v["events_score"], "d", "data/results/p_s_results.json", "TEST events")
    R(f"s1.evtr.{t}", FJ["events_train"][str(t)], "d", "data/results/p_s_fit.json", "TRAIN events")
def comp(prefix, C, taus=TAUS4, strat=True):
    for t in list(taus) + ["pool"]:
        v = C["pooled"] if t == "pool" else C["per_tau"][str(t)]
        R(f"{prefix}.d.{t}", 100 * v["delta"], "+.1f", "data/results/p_s_results.json", "relative held-out log-loss reduction")
        R(f"{prefix}.lo.{t}", 100 * v["ci"][0], "+.1f", "data/results/p_s_results.json", "95% seed-cluster CI lower")
        R(f"{prefix}.hi.{t}", 100 * v["ci"][1], "+.1f", "data/results/p_s_results.json", "95% seed-cluster CI upper")
        if strat and t != "pool" and "dauc_within" in v.get("stratified", {}):
            st = v["stratified"]
            R(f"{prefix}.w.{t}", st["dauc_within"], "+.3f", "data/results/p_s_results.json", "within-bin ΔAUC (event-weighted)")
            R(f"{prefix}.wlo.{t}", st["ci"][0], "+.3f", "data/results/p_s_results.json", "CI lower")
            R(f"{prefix}.whi.{t}", st["ci"][1], "+.3f", "data/results/p_s_results.json", "CI upper")
            R(f"{prefix}.nb.{t}", st["n_qualifying_bins"], "d", "data/results/p_s_results.json", "qualifying bins")
comp("qact", P1["comparisons"]["Q-ACT"])
comp("qactc", P1["comparisons"]["Q-ACT-carryover"])
comp("qsp", P1["comparisons"]["Q-SP-FULL (R-IN)"])
comp("qspi", P1["comparisons"]["Q-SP-IMAGE (R-IN)"])
comp("qspo", P1["comparisons"]["Q-SP-IMAGE (R-OUT)"])
comp("clin", P1["comparisons"]["CLIN-IMAGE (no R, no C)"])
for key, lab in (("knots8", "k8"), ("crossfit", "xf")):
    for cn, pre in (("Q-ACT", "qact"), ("Q-SP-FULL (R-IN)", "qsp")):
        v = S[key]["comparisons"][cn]["pooled"]
        R(f"{lab}.{pre}", 100 * v["delta"], "+.1f", f"data/results/p_s_results.json [{key}]", "pooled relative log-loss reduction")
for t in list(TAUS4) + ["pool"]:
    a = AB["Q-ACT"][str(t)]; b = AB["Q-SP-FULL (R-IN)"][str(t)]
    R(f"abs.base.{t}", 100 * a["base_rate"], ".3f", "data/results/p_s_absolute_scale.json", "event base rate among TEST hold epochs, %")
    R(f"abs.null.{t}", a["ll_null"], ".4f", "data/results/p_s_absolute_scale.json", "base-rate log-loss, nats/epoch")
    R(f"abs.a0.{t}", a["ll_base"], ".4f", "data/results/p_s_absolute_scale.json", "A0 held-out log-loss")
    R(f"abs.a1.{t}", a["ll_test"], ".4f", "data/results/p_s_absolute_scale.json", "A1 held-out log-loss")
    R(f"abs.s0.{t}", b["ll_base"], ".5f", "data/results/p_s_absolute_scale.json", "S0 held-out log-loss")
    R(f"abs.s1.{t}", b["ll_test"], ".5f", "data/results/p_s_absolute_scale.json", "S1-FULL held-out log-loss")
R("abs.br.a0", AB["Q-ACT"]["pool"]["brier_base"], ".5f", "data/results/p_s_absolute_scale.json", "pooled Brier A0")
R("abs.br.a1", AB["Q-ACT"]["pool"]["brier_test"], ".5f", "data/results/p_s_absolute_scale.json", "pooled Brier A1")
R("abs.below", 100 * (1 - AB["Q-ACT"]["pool"]["ll_base"] / AB["Q-ACT"]["pool"]["ll_null"]), ".0f", "data/results/p_s_absolute_scale.json", "A0 log-loss below null, %")
for t in list(TAUS4) + ["pool"]:
    for cn, pre in (("history gain, R in both (A1R vs A0R)", "rd1"), ("R alone over N+threshold (A0R vs A0)", "rd0")):
        C = RD["comparisons"][cn]; v = C["pooled"] if t == "pool" else C["per_tau"][str(t)]
        R(f"{pre}.d.{t}", 100 * v["delta"], "+.1f", "data/results/p_s_rdiag.json", "relative log-loss reduction")
        if t == "pool":
            R(f"{pre}.lo", 100 * v["ci"][0], "+.1f", "data/results/p_s_rdiag.json", "CI lower")
            R(f"{pre}.hi", 100 * v["ci"][1], "+.1f", "data/results/p_s_rdiag.json", "CI upper")
R("rd1.w.28", RD["comparisons"]["history gain, R in both (A1R vs A0R)"]["per_tau"]["28"]["stratified"]["dauc_within"], "+.3f",
  "data/results/p_s_rdiag.json", "within-bin ΔAUC, tau = 28, R in both")
for t in list(TAUS4) + ["pool"]:
    v = VP["loglik_recovery_of_A1_gain"][str(t)]
    R(f"vel.rec.{t}", 100 * v["recovery"], ".1f", "data/results/p_s_vel_prediction.json", "share of A1's log-loss gain recovered by VEL")
R("vel.rec.lo", 100 * VP["loglik_recovery_of_A1_gain"]["pool"]["ci"][0], ".1f", "data/results/p_s_vel_prediction.json", "CI lower")
R("vel.rec.hi", 100 * VP["loglik_recovery_of_A1_gain"]["pool"]["ci"][1], ".1f", "data/results/p_s_vel_prediction.json", "CI upper")

# policy test
PR = json.load(open(RES / "p_s_policy_results.json")); GJ = json.load(open(RES / "p_s_policy" / "gate.json"))
for arm in ("SIZE", "CLIN", "HIST", "FULL", "VEL"):
    for t in TAUS4:
        v = T[f"{arm}|{t}"]
        R(f"tune.sel.{arm}.{t}", v["selected"], ".1f" if arm == "SIZE" else ".3f", "data/results/p_s_policy/tune_selected.json",
          "selected thr (SIZE) or log10 p* (rules)")
        R(f"tune.p.{arm}.{t}", v["p_prem_tune"], ".3f", "data/results/p_s_policy/tune_selected.json", "TUNE p_prem at the selected setting")
        p = PR["per"][f"{arm}|{t}"]
        R(f"sc.p.{arm}.{t}", p["p_prem"], ".3f", "data/results/p_s_policy_results.json", "SCORE p_prem (1,000 seeds)")
        R(f"sc.f.{arm}.{t}", p["on_drug_fraction"], ".3f", "data/results/p_s_policy_results.json", "mean on-drug fraction")
        R(f"sc.rm.{arm}.{t}", p["rmst750"], ".1f", "data/results/p_s_policy_results.json", "RMST to model-defined progression, 750 d")
        if arm != "SIZE":
            for k, key, f_ in (("df", "df", "+.4f"), ("dr", "drmst", "+.1f"), ("dp", "dp", "+.3f")):
                R(f"sc.{k}.{arm}.{t}", p[key]["mean"], f_, "data/results/p_s_policy_results.json", f"paired Δ vs SIZE ({key})")
                R(f"sc.{k}lo.{arm}.{t}", p[key]["ci"][0], f_, "data/results/p_s_policy_results.json", "CI lower")
                R(f"sc.{k}hi.{arm}.{t}", p[key]["ci"][1], f_, "data/results/p_s_policy_results.json", "CI upper")
            R(f"sc.rel.{arm}.{t}", 100 * p["df"]["rel"], "+.1f", "data/results/p_s_policy_results.json", "Δf relative to SIZE, %")
            R(f"sc.read.{arm}.{t}", p["reading"], "s", "data/results/p_s_policy_results.json", "pre-registered reading")
    if arm != "SIZE":
        d = PR["pooled"][arm]
        R(f"pool.rel.{arm}", 100 * d["df"]["rel"], "+.1f", "data/results/p_s_policy_results.json", "pooled Δf relative, %")
        R(f"pool.rello.{arm}", 100 * d["df"]["rel_ci"][0], "+.1f", "data/results/p_s_policy_results.json", "CI lower")
        R(f"pool.relhi.{arm}", 100 * d["df"]["rel_ci"][1], "+.1f", "data/results/p_s_policy_results.json", "CI upper")
        R(f"pool.df.{arm}", d["df"]["mean"], "+.4f", "data/results/p_s_policy_results.json", "pooled Δf absolute")
        R(f"pool.dr.{arm}", d["drmst"]["mean"], "+.1f", "data/results/p_s_policy_results.json", "pooled ΔRMST, d")
        R(f"pool.drlo.{arm}", d["drmst"]["ci"][0], "+.1f", "data/results/p_s_policy_results.json", "CI lower")
        R(f"pool.drhi.{arm}", d["drmst"]["ci"][1], "+.1f", "data/results/p_s_policy_results.json", "CI upper")
        R(f"pool.dp.{arm}", d["dp"]["mean"], "+.3f", "data/results/p_s_policy_results.json", "pooled Δp_prem")
        R(f"pool.dplo.{arm}", d["dp"]["ci"][0], "+.3f", "data/results/p_s_policy_results.json", "CI lower")
        R(f"pool.dphi.{arm}", d["dp"]["ci"][1], "+.3f", "data/results/p_s_policy_results.json", "CI upper")
        R(f"pool.read.{arm}", d["reading"], "s", "data/results/p_s_policy_results.json", "pooled reading")
        R(f"pool.stands.{arm}", "stands" if d["headline_stands"] else "does not stand", "s", "data/results/p_s_policy_results.json",
          f"headline rule; supporting tau = {d['supporting_tau']}")
        R(f"pool.sup.{arm}", ", ".join(str(x) for x in d["supporting_tau"]) or "none", "s", "data/results/p_s_policy_results.json", "supporting tau")
for t in TAUS4:
    R(f"oe.CLIN.{t}", PR["diagnostics"][f"CLIN|{t}"]["O_over_E"], ".2f", "data/results/p_s_policy_results.json", "self-calibration O/E")
for k, v in PR["vel_recovery"].items():
    R(f"rho.{k}", v["rho"], ".2f", "data/results/p_s_policy_results.json", "Δf(VEL)/Δf(HIST)")
    R(f"rho.read.{k}", v["reading"], "s", "data/results/p_s_policy_results.json", "recovery reading")
R("vel7.wlo", PR["per"]["VEL|7"]["wilson"][0], ".3f", "data/results/p_s_policy_results.json", "Wilson lower, VEL tau=7")
R("gate.tune_proj", GJ["tune_projection_min_central_pilot"], ".1f", "data/results/p_s_policy/gate.json", "pilot TUNE projection, min")
R("gate.score_proj", GJ["score_projection_min_central_pilot"], ".1f", "data/results/p_s_policy/gate.json", "pilot SCORE projection, min")
R("gate.nchk", GJ["item3_epochs_checked"], "d", "data/results/p_s_policy/gate.json", "epochs feature-checked")
W = [json.loads(l) for l in open(RES / "p_s_policy" / "wall_log.jsonl")]
R("wall.tune", sum(x["wall_s"] for x in W if x["stage"] == "tune") / 60, ".1f", "data/results/p_s_policy/wall_log.jsonl", "summed invocation wall, min")
R("wall.score", sum(x["wall_s"] for x in W if x["stage"] == "score") / 60, ".1f", "data/results/p_s_policy/wall_log.jsonl", "summed invocation wall, min")
W1 = [json.loads(l) for l in open(RES / "p_s_cells" / "wall_log.jsonl")]
R("wall.step1", sum(x["wall_s"] for x in W1) / 60, ".1f", "data/results/p_s_cells/wall_log.jsonl", "summed invocation wall, min")
R("s2.nruns_score", 5 * 4 * 1000, "d", "data/results/p_s_policy/score (asserted at analyse)", "5 arms x 4 tau x 1000")
R("s2.nruns_tune", sum(len(v["evals"]) for v in T.values()) * 500, "d", "data/results/p_s_policy/tune_selected.json", "settings evaluated x 500")


# values quoted from dated records
OLD = "config/filings/assets_kit_2026-09-23_excerpts.md"
F("cred.400", "400", OLD, "all 400 held-out seeds", "[old §4.1.1]")
F("cred.res", "91–125", RL, "resolution 91–125% of the term", "CHECK A")
F("cred.basis", "179.9", RL, "179.9 ms/run** for a 750-d ABM AT cell", "RECAL row cost basis")
F("cred.kl", "0.377", OLD, "KL 0.377", "[old §4.5]")
F("cred.chi", "17,558", OLD, "17,558", "[old §4.5]")
F("cred.pp", "26.2", OLD, "+26.2 pp", "[old §4.5]")
R("cred.gate173", json.load(open(RES / "p_s_cells" / "gate.json"))["ms_per_run_ref"], ".1f", "data/results/p_s_cells/gate.json", "Step 1 gate, uninstrumented ms/run under load at the time")
F("lit.en20", "20", "config/filings/lit_check_2026-09-23.json", "dose -50% if volume fell >=20%", "Enriquez-Navas 2016 AT-1, PMC4962860")
F("lit.en50", "50", "config/filings/lit_check_2026-09-23.json", "dose -50% if volume fell >=20%", "Enriquez-Navas 2016 AT-1, PMC4962860")

# design constants
# Grid, seed ranges, model parameters and planned settings, read from the source code with ast (nothing is
# executed) or checked word for word in the analysis plan.
import ast
def const(f, name):
    """Value of the module-level constant `name` in source file f (relative to src/threshold_transfer), read with ast."""
    for node in ast.parse(open(SRC / f).read()).body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            tg = node.targets if isinstance(node, ast.Assign) else [node.target]
            if any(getattr(t, "id", None) == name for t in tg):
                return ast.literal_eval(node.value)
    raise KeyError(f"{name} not in {f}")
def sim_default(name):
    """Default of field `name` in simulator.SimulatorParams, read with ast."""
    for node in ast.walk(ast.parse(open(SRC / "model" / "simulator.py").read())):
        if isinstance(node, ast.ClassDef) and node.name == "SimulatorParams":
            for s in node.body:
                if isinstance(s, ast.AnnAssign) and s.target.id == name: return ast.literal_eval(s.value)
    raise KeyError(name)
for t in TAUS6:
    R(f"design.tau.{t}", t, "d", "src/threshold_transfer/threshold/grid.py", "TAUS")
R("design.q", G.Q, ".2f", "src/threshold_transfer/threshold/grid.py", "Q, arbitrary model-reliability level")
R("design.R", G.R, "d", "src/threshold_transfer/threshold/grid.py", "R, seeds per main-run cell")
R("design.progmul", P.PROG_MUL, ".1f", "src/threshold_transfer/model/ode.py", "PROG_MUL, progression limit = PROG_MUL * N0")
R("design.mstep", G.CORE[1] - G.CORE[0], ".2f", "src/threshold_transfer/threshold/grid.py", "multiplier step of the core grid")
R("design.mlo", G.CORE[0], ".2f", "src/threshold_transfer/threshold/grid.py", "lowest multiplier")
for key, f, nm in (("main", "threshold/grid.py", "BLOCK_ABM"), ("recal", "experiments/recal_map.py", "BLOCK"),
                   ("rctmc", "experiments/recal_ctmc.py", "BLOCK"), ("fine", "experiments/recal_ctmc_fine.py", "BLOCK")):
    R(f"design.block.{key}", const(f, nm), "d", f"src/threshold_transfer/{f}", nm)
for role, b in const("experiments/statesuff_log.py", "BLOCKS").items():
    R(f"design.block.s1{role}", b, "d", "src/threshold_transfer/experiments/statesuff_log.py", f"BLOCKS['{role}']")
for role, b in const("experiments/policy.py", "BLOCKS").items():
    R(f"design.block.s2{role}", b, "d", "src/threshold_transfer/experiments/policy.py", f"BLOCKS['{role}']")
F("design.block.reserved", "1700000", "config/preregistrations/PREREG_2026-09-23_state_sufficiency_logging.md", "1700000 stays reserved and untouched", "S1, reserved block")
R("design.ctmcstep", const("experiments/recal_ctmc_fine.py", "STEP"), ".4f", "src/threshold_transfer/experiments/recal_ctmc_fine.py", "STEP, fine CTMC multiplier step")
from threshold_transfer.model.ode import ODE2Params
_op = ODE2Params()
R("par.rS", sim_default("r_S"), ".3f", "src/threshold_transfer/model/simulator.py", "SimulatorParams.r_S, division rate S, /day")
R("par.rR", sim_default("r_Ra"), ".3f", "src/threshold_transfer/model/simulator.py", "SimulatorParams.r_Ra, division rate resistant, /day")
R("par.d", sim_default("base_death"), ".3f", "src/threshold_transfer/model/simulator.py", "SimulatorParams.base_death, /day")
R("par.kill", sim_default("d_A_S"), ".2f", "src/threshold_transfer/model/simulator.py", "SimulatorParams.d_A_S, kill coefficient of drug on S")
R("par.D", sim_default("D_A"), ".2f", "src/threshold_transfer/model/simulator.py", "SimulatorParams.D_A, drug diffusion coefficient")
R("par.kap", sim_default("kappa_A"), ".2f", "src/threshold_transfer/model/simulator.py", "SimulatorParams.kappa_A, systemic clearance/supply")
R("par.kapcell", sim_default("kappa_cell_A"), ".2f", "src/threshold_transfer/model/simulator.py", "SimulatorParams.kappa_cell_A, uptake by occupied sites")
R("par.dt", sim_default("dt"), ".0f", "src/threshold_transfer/model/simulator.py", "SimulatorParams.dt, days")
R("par.c", _op.growth_params[0], ".3f", "src/threshold_transfer/model/ode.py", "growth_params[0], perimeter coefficient c")
R("par.phib", _op.growth_params[1], ".4f", "src/threshold_transfer/model/ode.py", "growth_params[1], bulk-turnover floor phi_bulk")
R("par.ceff", _op.C_eff, ".4f", "src/threshold_transfer/model/ode.py", "C_eff = kappa/(kappa+kappa_cell)")
PM = "config/preregistrations/PREREG_2026-08-22_threshold_transfer_main.md"
F("pk.kpk", "0.10004", PM, "kappa_pk = 0.10004 /day", "A4.1 PK fit")
F("pk.kel", "0.13420", PM, "k_el = 0.13420 /day", "A4.1 PK fit")
F("pk.tc", "7.45", PM, "time constant 7.45 d", "A4.1 PK fit")
PL = "config/preregistrations/PREREG_2026-09-23_state_sufficiency_logging.md"; PP = "config/preregistrations/PREREG_2026-09-23_state_sufficiency_policy.md"
R("s1d.nthr", len(const("experiments/statesuff_log.py", "BAND")[7]), "d", "src/threshold_transfer/experiments/statesuff_log.py", "thresholds per tau in BAND")
F("s1d.R", "500", PL, "R = 500 per cell per block", "S2")
for t in TAUS4:
    R(f"s1d.ntop.{t}", const("experiments/statesuff_log.py", "BAND")[t][-1], "d", "src/threshold_transfer/experiments/statesuff_log.py", "top BAND threshold = N_top (Step 2 cap)")
F("s1d.knots", "5", PL, "n_knots=5", "Amendment C")
F("s1d.knots8", "8", PL, "n_knots=8", "Amendment C sensitivity refit")
F("s1d.ridge", "1.0", PL, "C=1.0", "Amendment C")
R("s1d.boot", const("analysis/statesuff.py", "NBOOT"), ",d", "src/threshold_transfer/analysis/statesuff.py", "NBOOT, cluster-bootstrap resamples")
F("s1d.gain", "10", PL, "**GAIN:** Δ ≥ 0.10", "S8, 10% relative log-loss bar")
F("s1d.binw", "20", PL, "binned in N at width 20 cells", "S9")
F("s1d.binmin", "10", PL, "**≥ 10 events and ≥ 10", "Amendment B")
F("s1d.minbins", "3", PL, "fewer than 3 bins", "Amendment B")
F("design.ci", "95", PP, "percentile 95% CIs", "P6, CI level")
F("s2d.Rtune", "500", PP, "(R = 500 per setting)", "P5")
F("s2d.Rscore", "1,000", PP, "(R = 1000)", "P6")
F("s2d.neval", "8", PP, "8 evaluations per arm", "P5")
F("s2d.ni", "0.02", PP, "≤ +0.02", "sign-off #7, arbitrary margin")
F("s2d.rho", "0.80", PP, "ρ ≥ 0.80", "Amendment 1 recovery bar")
F("s2d.match", "0.010", PP, "more than 0.010", "P5 matching check")
R("pool.dfpp.CLIN", 100 * PR["pooled"]["CLIN"]["df"]["mean"], "+.1f", "data/results/p_s_policy_results.json",
  "pooled Δ on-drug fraction, percentage points")


DESIGN_PREFIXES = ("design.", "par.", "pk.", "s1d.", "s2d.")     # design constants


def check():
    """Counts of entries by kind, filed values reconciled exactly or within rounding, and the list of failures
    (missing filed text, or a recomputed value outside the tolerance of its filed value).
    """
    out = dict(recomputed=0, filing=0, reconciled=0, reconciled_within_rounding=0, recomputed_new=0,
               design_from_code=0, design_from_filing=0, failing=[])
    for k, e in L.items():
        if k.startswith(DESIGN_PREFIXES):
            out["design_from_code" if e["kind"] == "recomputed" else "design_from_filing"] += 1
            if e["kind"] == "filing":
                t = open(ROOT / e["artifact"]).read()
                if e["substring"] not in t: out["failing"].append((k, "substring not found in filing"))
            continue
        if e["kind"] == "filing":
            out["filing"] += 1
            t = open(ROOT / e["artifact"]).read()
            if e["substring"] not in t: out["failing"].append((k, "substring not found in filing"))
        else:
            out["recomputed"] += 1
            if e["filed"]:
                f, sub, fv, tol = e["filed"]
                t = open(ROOT / f).read()
                if sub not in t: out["failing"].append((k, f"filed text '{sub}' not found in {f}"))
                elif abs(e["value"] - fv) <= 1e-9 + 0.5 * 10 ** -(len(str(fv).split('.')[-1]) if '.' in str(fv) else 0):
                    out["reconciled"] += 1
                elif abs(e["value"] - fv) <= tol:
                    out["reconciled_within_rounding"] += 1
                else:
                    out["failing"].append((k, f"recomputed {e['value']:.4f} vs filed {fv}"))
            else:
                out["recomputed_new"] += 1
    return out


def run():
    """Write data/results/ledger.json and print the counts."""
    c = check()
    json.dump(dict(entries=L, counts=c), open(RES / "ledger.json", "w"), indent=1, default=str)
    print(json.dumps(c, indent=1))
    return c
