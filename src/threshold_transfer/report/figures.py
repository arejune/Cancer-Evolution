"""Figures of the paper (`python run.py figures`), from data/results and the number ledger; no simulation.

Writes figures/F1-F5 as PDF and 300-dpi PNG, and figures/PROVENANCE.md with the MD5 of every output and source.
"""
import json, hashlib, platform
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from threshold_transfer.paths import ROOT, SRC, FIGURES
from threshold_transfer.report import ledger as WL

RES = WL.RES; OUT = FIGURES
Q = WL.Q; L = WL.L
TAUS6, TAUS4 = WL.TAUS6, WL.TAUS4

# colours: four categorical colours in fixed order, neutral greys, and a blue ramp for tau (Fig. 3)
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
TAU_RAMP = dict(zip(TAUS4, ("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")))
# figures are drawn at the single-column text width of the Springer Nature template (31 pc, 131 mm) and placed at
# their natural size, so every label prints at 8 pt, the smallest size Springer artwork guidelines allow.
WIDTH = 5.1
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 1.2,
    "lines.markersize": 3.5, "pdf.fonttype": 42, "savefig.dpi": 300, "figure.dpi": 100,
    "svg.hashsalt": "0", "axes.titlesize": 8, "mathtext.fontset": "dejavusans",
})
PROV = {}
# paper symbols (mathtext) and constants read from the ledger
PPREM, NSTAR, NSAFE = r"$p_{\mathrm{prem}}$", r"$N^{*}$", r"$N_{\mathrm{safe}}$"
QTXT, GAIN, NI = L["design.q"]["text"], float(L["s1d.gain"]["text"]), float(L["s2d.ni"]["text"])
MSTEP, TEND = L["design.mstep"]["text"], L["cfg.tend"]["text"]


def md5(p):
    """MD5 of a file."""
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def save(fig, name, sources, ledger_ids, spec):
    """Write figure `name` as PDF and 300-dpi PNG with fixed metadata, and record its sources and ledger ids."""
    missing = [i for i in ledger_ids if i not in L]; assert not missing, missing
    for ext, kw in (("pdf", dict(metadata={"CreationDate": None, "ModDate": None, "Producer": None})),
                    ("png", dict(dpi=300, metadata={"Software": None}))):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight", pad_inches=0.02, **kw)
    plt.close(fig)
    PROV[name] = dict(sources=sources, ledger=sorted(ledger_ids), spec=spec)


def tau_tag(ax, t):
    """Label a panel with its monitoring interval."""
    ax.text(0.04, 0.96, f"τ = {t} d", transform=ax.transAxes, ha="left", va="top", color=INK2)


def panel(ax, letter, x=-0.02):
    """Bold lower-case panel letter above the top-left corner of the axes, as in Springer Nature figures."""
    ax.text(x, 1.03, letter, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, fontweight="bold", color=INK)


def qline(ax):
    """Horizontal line at the tolerance q."""
    ax.axhline(Q, color=MUTED, lw=0.7, ls="-", zorder=1)




def main():
    """Draw F1-F5 and write PROVENANCE.md."""
    # F1: p_prem against the threshold, four models
    OUT.mkdir(parents=True, exist_ok=True)
    rec = WL.rec; NS = WL.NS
    arms = json.load(open(RES / "p_t_main_arms.json"))
    for t in TAUS6:
        assert abs(arms["nstar"][str(t)] - NS[t]) < 1e-6, t
    fig, axs = plt.subplots(3, 2, figsize=(WIDTH, 5.4), sharey=True)
    for ax, t, letter in zip(axs.flat, TAUS6, "abcdef"):
        s = rec[(rec["pol"] == "AT") & (rec["tau"] == t)]; ms = sorted(set(s["m"]))
        ax.plot([m * NS[t] for m in ms], [s[s["m"] == m]["prem"].mean() for m in ms], color=C1, marker="o", zorder=4)
        for arm, col, mk, ls in (("ctmc_nom", C2, "s", "-"), ("ctmc_pk", C3, "^", "-"), ("ode", INK2, None, (0, (2, 1.5)))):
            ks = sorted((k for k in arms["cells"] if k.startswith(f"{arm}|{t}|")), key=lambda k: float(k.split("|")[2]))
            ax.plot([float(k.split("|")[2]) * NS[t] for k in ks], [np.mean(arms["cells"][k]["prem"]) for k in ks],
                    color=col, marker=mk, ls=ls, lw=1.0 if arm == "ode" else 1.2, zorder=3)
        qline(ax)
        ax.axvline(NS[t], color=INK, lw=0.7, ls=(0, (6, 2, 1, 2)), zorder=2)
        ax.axvline(WL.main[t], color=C1, lw=0.8, ls=":", zorder=2)
        tau_tag(ax, t); panel(ax, letter)
        ax.set_ylim(-0.03, 1.03)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    for ax in axs[-1]: ax.set_xlabel("Treatment-hold threshold (cells)")
    for ax in axs[:, 0]: ax.set_ylabel(PPREM)
    fig.legend(handles=[Line2D([], [], color=C1, marker="o", label="ABM"),
                        Line2D([], [], color=C2, marker="s", label="CTMC-nominal"),
                        Line2D([], [], color=C3, marker="^", label="CTMC-pk"),
                        Line2D([], [], color=INK2, ls=(0, (2, 1.5)), label="ODE"),
                        Line2D([], [], color=INK, lw=0.7, ls=(0, (6, 2, 1, 2)), label=f"{NSTAR} (derived threshold)"),
                        Line2D([], [], color=C1, lw=0.8, ls=":", label=f"ABM {NSAFE} (safe threshold)"),
                        Line2D([], [], color=MUTED, lw=0.7, label=f"tolerance q = {QTXT}")],
               loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, handlelength=2.0, columnspacing=1.0)
    fig.tight_layout(h_pad=1.2)
    save(fig, "F1", ["data/results/p_t_main_raw.npz", "data/results/p_t_main_arms.json", "data/results/p_t_nstar.json"],
         [f"nstar.{t}" for t in TAUS6] + [f"main.nsafe.{t}" for t in TAUS6],
         "p_prem vs threshold per τ: ABM, CTMC-nominal, CTMC-pk, ODE; N* and N_safe marked; q = 0.05 line")

    # F2: safe thresholds of the ODE, CTMC and ABM under the recalibrated law
    fig, ax = plt.subplots(figsize=(WIDTH, 3.3))
    x = np.arange(len(TAUS6)); ids = []
    for i, t in enumerate(TAUS6):
        ys = [WL.odeb[t], WL.ctmc[t], WL.recal[t]]
        ax.plot([i, i], [min(ys), max(ys)], color=AXIS, lw=0.8, zorder=1)
        for y, col, mk in zip(ys, (INK2, C2, C1), ("D", "s", "o")):
            ax.plot(i, y, marker=mk, color=col, ms=5, ls="none", zorder=3)
        ids += [f"odeR.{t}", f"ctmc.{t}", f"recal.nsafe.{t}"]
        if t in TAUS4:
            ax.plot(i + 0.18, WL.s1[t], marker="o", mfc="white", mec=C1, mew=1.0, ms=5, ls="none", zorder=4)
            ax.annotate("", xy=(i + 0.18, WL.s1[t]), xytext=(i + 0.02, WL.recal[t]),
                        arrowprops=dict(arrowstyle="-", color=C1, lw=0.5, ls=":"))
            ax.text(i, max(ys) + 60, f"{L[f'share.{t}']['text']}%\n({L[f'corr.share.{t}']['text']}%)",
                    ha="center", va="bottom", color=INK2, linespacing=1.1)
            ids += [f"corr.s1.{t}", f"share.{t}", f"corr.share.{t}"]
        else:
            ax.text(i, max(ys) + 60, "no share", ha="center", va="bottom", color=MUTED)
    ax.set_xticks(x, [str(t) for t in TAUS6])
    ax.set_xlim(-0.5, len(TAUS6) - 0.3)
    ax.set_xlabel("Monitoring interval τ (days)")
    ax.set_ylabel("Safe threshold (cells)")
    lo = min(min(WL.recal.values()), min(WL.s1.values())); hi = max(max(WL.odeb.values()), max(WL.ctmc.values()))
    ax.set_ylim(lo - 60, hi + 280)
    ax.legend(handles=[Line2D([], [], color=INK2, marker="D", ls="none", label="ODE"),
                       Line2D([], [], color=C2, marker="s", ls="none", label="CTMC"),
                       Line2D([], [], color=C1, marker="o", ls="none", label="ABM, recalibration experiment"),
                       Line2D([], [], color=C1, marker="o", mfc="white", ls="none", label="ABM, finer grid (corrected)")],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, columnspacing=1.0, handletextpad=0.3)
    fig.tight_layout()
    save(fig, "F2", ["data/results/p_t_recal_odeb.json", "data/results/p_t_recal_ctmc_fine.json", "data/results/p_t_recal_ctmc_fine2.json",
                     "data/results/p_t_recal_sweep.npz", "data/results/p_s_cells/train/*.npz", "data/results/p_s_cells/test/*.npz"], ids,
         "Safe thresholds of the ODE, CTMC and ABM per τ under the recalibrated law, spatial share annotated; corrected ABM value overlaid")

    # F3: prediction gains on the test set
    S = json.load(open(RES / "p_s_results.json"))["primary"]["comparisons"]
    RD = json.load(open(RES / "p_s_rdiag.json"))["comparisons"]
    ROWS = [("History over size", S["Q-ACT"], "qact"),
            ("History over size,\ndrug level in both", S["Q-ACT-carryover"], "qactc"),
            ("Geometry over history,\ndrug level, resistant count", S["Q-SP-FULL (R-IN)"], "qsp"),
            ("Image geometry over history,\ndrug level, resistant count", S["Q-SP-IMAGE (R-IN)"], "qspi"),
            ("Image geometry over\nhistory, drug level", S["Q-SP-IMAGE (R-OUT)"], "qspo"),
            ("Image geometry\nover history", S["CLIN-IMAGE (no R, no C)"], "clin"),
            ("History over size,\nresistant count in both", RD["history gain, R in both (A1R vs A0R)"], "rd1"),
            ("Resistant count\nover size", RD["R alone over N+threshold (A0R vs A0)"], "rd0")]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(WIDTH, 4.6), sharey=True, gridspec_kw=dict(width_ratios=(1.15, 1)))
    off = dict(zip(TAUS4, (-0.27, -0.13, 0.01, 0.15))); ids = []
    for r, (lab, C, pre) in enumerate(ROWS):
        for t in TAUS4:
            v = C["per_tau"][str(t)]; y = r + off[t]
            a1.errorbar(100 * v["delta"], y, xerr=[[100 * (v["delta"] - v["ci"][0])], [100 * (v["ci"][1] - v["delta"])]],
                        fmt="o", color=TAU_RAMP[t], ms=3.2, elinewidth=0.9, capsize=0, zorder=3)
            st = v.get("stratified", {})
            if "dauc_within" in st:
                a2.errorbar(st["dauc_within"], y, xerr=[[st["dauc_within"] - st["ci"][0]], [st["ci"][1] - st["dauc_within"]]],
                            fmt="o", color=TAU_RAMP[t], ms=3.2, elinewidth=0.9, capsize=0, zorder=3)
        v = C["pooled"]; y = r + 0.29
        a1.errorbar(100 * v["delta"], y, xerr=[[100 * (v["delta"] - v["ci"][0])], [100 * (v["ci"][1] - v["delta"])]],
                    fmt="D", color=INK, ms=3.2, elinewidth=0.9, capsize=0, zorder=3)
        ids += [f"{pre}.d.{t}" for t in list(TAUS4) + ["pool"]]
    a1.axvline(0, color=AXIS, lw=0.8, zorder=1); a2.axvline(0, color=AXIS, lw=0.8, zorder=1)
    a1.axvline(GAIN, color=INK2, lw=0.6, ls=(0, (4, 3)), zorder=1)
    a1.set_yticks(range(len(ROWS)), [r[0] for r in ROWS]); a1.invert_yaxis()
    a1.tick_params(axis="y", labelsize=8)
    a1.set_xlabel("Reduction in test log-loss (%)")
    a2.set_xlabel("ΔAUC within size bins")
    for a in (a1, a2): a.grid(axis="y", visible=False)
    panel(a1, "a", x=0.0); panel(a2, "b", x=0.0)
    fig.legend(handles=[Line2D([], [], color=TAU_RAMP[t], marker="o", ls="none", label=f"τ = {t} d") for t in TAUS4]
               + [Line2D([], [], color=INK, marker="D", ls="none", label="pooled"),
                  Line2D([], [], color=INK2, lw=0.6, ls=(0, (4, 3)), label=f"gain threshold ({L['s1d.gain']['text']}%)")],
               loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, columnspacing=1.0, handletextpad=0.3)
    fig.tight_layout()
    save(fig, "F3", ["data/results/p_s_results.json", "data/results/p_s_rdiag.json"], ids,
         "Prediction study: reduction in test log-loss [CI] and ΔAUC within size bins [CI] per τ for each comparison")

    # F4: policy test, paired differences against the size-only rule
    PR = json.load(open(RES / "p_s_policy_results.json"))
    RULES = (("CLIN", "history and image (primary)", C1, "o"), ("HIST", "history", C2, "s"), ("FULL", "full state", C3, "^"),
             ("VEL", "growth rate", C4, "D"))
    METRICS = (("df", "Δ on-drug fraction"), ("drmst", "Δ RMST (days)"), ("dp", f"Δ {PPREM}"))
    GROUPS = list(TAUS4) + ["pooled"]
    fig, axs = plt.subplots(3, 1, figsize=(WIDTH, 5.0), sharex=True)
    dx = dict(zip([r[0] for r in RULES], (-0.24, -0.08, 0.08, 0.24))); ids = []; readings = set()
    for ax, (key, xl), letter in zip(axs, METRICS, "abc"):
        for arm, lab, col, mk in RULES:
            for gi, g in enumerate(GROUPS):
                p = PR["pooled"][arm] if g == "pooled" else PR["per"][f"{arm}|{g}"]
                m, (lo, hi) = p[key]["mean"], p[key]["ci"]; rd = p["reading"]; readings.add(rd)
                filled = rd.startswith("GAIN")
                ax.errorbar(gi + dx[arm], m, yerr=[[m - lo], [hi - m]], fmt=mk, color=col, mfc=col if filled else "white",
                            mec=col, mew=1.0, ms=3.6, elinewidth=0.9, capsize=0, zorder=3)
                if rd.startswith("NOT ASSESSABLE"):
                    ax.plot(gi + dx[arm], m, marker="x", color=INK, ms=3.6, mew=0.9, ls="none", zorder=4)
        ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
        ax.set_xticks(range(len(GROUPS)), [str(g) for g in TAUS4] + ["pooled"])
        ax.set_ylabel(xl)
        ax.grid(axis="x", visible=False)
        panel(ax, letter, x=-0.1)
    axs[2].axhline(NI, color=INK2, lw=0.6, ls=(0, (4, 3)), zorder=1)
    axs[2].set_xlabel("Monitoring interval τ (days)")
    for arm in [r[0] for r in RULES]:
        ids += [f"sc.{k}.{arm}.{t}" for k in ("df", "dr", "dp") for t in TAUS4] + [f"pool.{k}.{arm}" for k in ("df", "dr", "dp")]
    handles = [Line2D([], [], color=c, marker=mk, ls="none", label=lab) for _, lab, c, mk in RULES]
    handles += [Line2D([], [], color=INK2, marker="o", ls="none", label="filled: gain"),
                Line2D([], [], color=INK2, marker="o", mfc="white", ls="none", label="open: no gain"),
                Line2D([], [], color=INK, marker="x", ls="none", label="cross: not assessable"),
                Line2D([], [], color=INK2, lw=0.6, ls=(0, (4, 3)), label=f"non-inferiority margin (+{L['s2d.ni']['text']})")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, columnspacing=1.2, handletextpad=0.3)
    fig.align_ylabels(axs)
    fig.tight_layout()
    save(fig, "F4", ["data/results/p_s_policy_results.json"], ids,
         "Policy test: paired Δ on-drug fraction, ΔRMST, Δp_prem [CI] against the size-only rule per τ for four rules; non-inferiority margin; readings "
         f"(marker fill; reading strings present: {sorted(readings)})")

    # F5: grid resolution of N_safe
    T = WL.T; r1 = WL.r1; BAND = WL.A.BAND
    fig, axs = plt.subplots(2, 2, figsize=(WIDTH, 4.2), sharey=True)
    ids = []
    for ax, t, letter in zip(axs.flat, TAUS4, "abcd"):
        s = rec[(rec["pol"] == "AT") & (rec["tau"] == t)]; ms = sorted(set(s["m"]))
        xm = np.array([m * NS[t] for m in ms]); ym = np.array([s[s["m"] == m]["prem"].mean() for m in ms])
        ax.plot(xm, ym, color=C1, marker="o", zorder=3)
        th = BAND[t]; y1 = [r1["prem"][(r1["tau"] == t) & (r1["thr"] == v)].mean() for v in th]
        ax.plot(th, y1, color=C2, marker="s", zorder=4)
        ev = np.array(T[f"SIZE|{t}"]["evals"]); o = np.argsort(ev[:, 0])
        ax.plot(ev[o, 0], ev[o, 1], color=C3, marker="^", ls="none", ms=4, zorder=5)
        # opaque equivalent of C3 at alpha 0.18 on white, so the PDF and EPS versions draw the same band
        ax.axvspan(WL.L[f"corr.tlo.{t}"]["value"], WL.L[f"corr.thi.{t}"]["value"], color="#d6f1e7", lw=0, zorder=0.5)
        ax.axvline(WL.main[t], color=C1, lw=0.8, ls=":", zorder=2)
        ax.axvline(WL.s1[t], color=C2, lw=0.8, ls=":", zorder=2)
        qline(ax)
        step = (ms[1] - ms[0]) * NS[t]
        allx = list(th) + list(ev[:, 0]) + [WL.main[t]]
        ax.set_xlim(min(allx) - 1.2 * step, max(allx) + 1.2 * step)
        ax.set_ylim(-0.03, 1.03); tau_tag(ax, t); panel(ax, letter)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ids += [f"main.nsafe.{t}", f"corr.s1.{t}", f"corr.tlo.{t}", f"corr.thi.{t}"]
    for ax in axs[1]: ax.set_xlabel("Treatment-hold threshold (cells)")
    for ax in axs[:, 0]: ax.set_ylabel(PPREM)
    fig.legend(handles=[Line2D([], [], color=C1, marker="o", label=f"main grid (step {MSTEP} {NSTAR})"),
                        Line2D([], [], color=C2, marker="s", label="finer grid"),
                        Line2D([], [], color=C3, marker="^", ls="none", label="tuning search"),
                        plt.Rectangle((0, 0), 1, 1, color="#d6f1e7", lw=0, label="final tuning interval"),
                        Line2D([], [], color=INK2, lw=0.8, ls=":", label="interpolated crossing"),
                        Line2D([], [], color=MUTED, lw=0.7, label=f"tolerance q = {QTXT}")],
               loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, columnspacing=1.2)
    fig.tight_layout()
    save(fig, "F5", ["data/results/p_t_main_raw.npz", "data/results/p_s_cells/train/*.npz", "data/results/p_s_cells/test/*.npz",
                     "data/results/p_s_policy/tune_selected.json"], ids,
         "N_safe resolution: p_prem vs threshold at τ = 7/14/28/56 on the main grid, the finer grid and the tuning bisection; crossings marked")

    lines = ["# Figure provenance", "",
             "Generated by `python run.py figures` (no simulation; every value read from the files below, "
             "directly or through `src/threshold_transfer/report/ledger.py`). No titles or captions are drawn.", "",
             f"- Script md5: `{md5(SRC / 'report' / 'figures.py')}`; ledger md5: `{md5(SRC / 'report' / 'ledger.py')}`",
             f"- Interpreter: Python {platform.python_version()}, numpy {np.__version__}, matplotlib {matplotlib.__version__}",
             "- Colours: four categorical colours in fixed order and a blue ramp for τ in F3; "
             "every series also carries a distinct marker or line style.", ""]
    for name, p in PROV.items():
        lines += [f"## {name}", "", f"Spec: {p['spec']}", "", "Outputs:"]
        for ext in ("pdf", "png"):
            lines.append(f"- `figures/{name}.{ext}` md5 `{md5(OUT / f'{name}.{ext}')}`")
        lines += ["", "Sources (md5 at generation):"]
        for s in p["sources"]:
            fs = sorted(ROOT.glob(s))
            if len(fs) == 1:
                lines.append(f"- `{s}` `{md5(fs[0])}`")
            else:
                h = hashlib.md5("".join(md5(f) for f in fs).encode()).hexdigest()
                lines.append(f"- `{s}` ({len(fs)} files; md5 of concatenated per-file md5s, sorted by path: `{h}`)")
        lines += ["", "Ledger ids whose values are plotted (or used for markers):", "", ", ".join(f"`{i}`" for i in p["ledger"]), ""]
    (OUT / "PROVENANCE.md").write_text("\n".join(lines))
    print("written", sorted(str(p.name) for p in OUT.iterdir()))
