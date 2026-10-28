"""Tables of the paper, written to tables/<name>.tex from the number ledger (`python run.py tables`).

Every number is a ledger entry ({L:id}); none is typed here. Captions belong to the paper.
"""
import json
from threshold_transfer.paths import RESULTS, TABLES
from threshold_transfer.report import richtext as RT
LED = json.load(open(RESULTS / "ledger.json"))
assert not LED["counts"]["failing"], "ledger has failing entries; rebuild it first"
T6, T4 = (7, 14, 28, 42, 56, 84), (7, 14, 28, 56)
FLOOR = (42, 84)   # intervals where the ODE threshold is set by the decision calendar
ARMS = (("CLIN", "history and image"), ("HIST", "history"), ("FULL", "full state"), ("VEL", "growth rate"))
# Readings and headline outcomes are ledger text in upper case (the analysis code's labels); tables print them in words.
# Every ledger text that occurs must be listed here, so a new reading cannot reach a table unmapped.
READING = {"GAIN (exposure + time)": "gain (exposure, time)", "GAIN (exposure)": "gain (exposure)",
           "GAIN (time)": "gain (time)", "NO GAIN": "no gain", "NO GAIN (worse)": "no gain (worse)",
           "NOT ASSESSABLE (bound fails on width)": "not assessable",
           "NOT ASSESSABLE (SAFETY FAILURE)": "safety failure",
           "stands": "met", "does not stand": "not met"}


def words(i):
    """A reading or headline outcome from the ledger, in plain words."""
    assert i in RT.LEDGER, i
    return READING[RT.LEDGER[i]["text"]]


def v(i):
    """Token for ledger entry i (asserts the entry exists)."""
    assert i in RT.LEDGER, i
    return "{L:" + i + "}"


def ci(lo, hi):
    """Interval token '[lo, hi]' from two ledger ids."""
    return f"[{v(lo)}, {v(hi)}]"
def tau(t):
    """Monitoring-interval label for a table row."""
    return v(f"design.tau.{t}") if t != "pool" else "pooled"


TAU = "τ (d)"
TABLES_DEF = {}


def table(name, header, rows, align=None, wide=False):
    """Register a table: header and rows of rich strings; checks that no plain text contains a digit."""
    for s in header + [c for r in rows for c in r]:
        RT.check_no_digits(s, name)
    assert all(len(r) == len(header) for r in rows), name
    TABLES_DEF[name] = dict(header=header, rows=rows, align=align or "r" * len(header), wide=wide)


# methods
table("tab_config", ["parameter", "symbol", "value", "unit"], [
    ["initial tumour size", "$N_0$", v("cfg.N0"), "cells"], ["initial disk radius", "", v("cfg.ric"), "sites"],
    ["lattice edge", "$L$", v("cfg.L"), "sites"], ["lattice sites", "$n_{\\max}$", v("cfg.nmax"), "sites"],
    ["progression limit", "$N_{\\mathrm{prog}}$", v("cfg.prog"), "cells"],
    ["pre-existing resistant fraction", "", v("cfg.rfrac"), "%"],
    ["simulation horizon", "", v("cfg.tend"), "days"], ["time step", "$\\Delta t$", v("par.dt"), "days"],
    ["division rate, sensitive cells", "$r_S$", v("par.rS"), "per day"],
    ["division rate, resistant cells", "$r_R$", v("par.rR"), "per day"],
    ["background death rate", "$d$", v("par.d"), "per day"],
    ["drug kill coefficient", "$k$", v("par.kill"), "per day per unit concentration"],
    ["drug diffusion coefficient", "$D$", v("par.D"), "lattice units"],
    ["drug exchange rate", "$\\kappa$", v("par.kap"), "per day"],
    ["drug uptake at occupied sites", "$\\kappa_{\\mathrm{cell}}$", v("par.kapcell"), "per day"],
    ["growth law, rim coefficient", "$c$", v("par.c"), ""],
    ["growth law, interior turnover", "$\\phi_{\\mathrm{bulk}}$", v("par.phib"), ""]],
    align="llrl")
table("tab_models", ["model", "information added", "terms in the model"], [
    ["size", "none (baseline)", "spline in $N$; threshold-cell indicator"],
    ["history", "earlier decisions", "size terms; days since treatment was last given; never-treated indicator; "
                                     "$\\log(N/N_{\\mathrm{prev}})$; number of earlier treated decisions; nadir since treatment"],
    ["state", "drug level and resistant count", "history terms; mean drug concentration $C$; resistant count $R$"],
    ["geometry", "shape of the tumour", "state terms; rim cells $B$; resistant rim cells $B_R$; components $K$; "
                                       "radius of gyration $R_g$; maximum radius $R_{\\max}$"],
    ["image geometry", "shape visible on an image", "geometry terms without $B_R$"]],
    align="llp")

# random-number seeds
table("tab_seeds", ["experiment", "first seed", "runs"], [
    ["main experiment, ABM and the ODE and CTMC arms", v("design.block.main"), v("main.nruns")],
    ["recalibration: drug-free transit times and threshold sweep", v("design.block.recal"), v("recal.nruns")],
    ["CTMC under the recalibrated law", v("design.block.rctmc"), "--"],
    ["CTMC under the recalibrated law, finer grids", v("design.block.fine"), "--"],
    ["prediction study, check of the logging code", v("design.block.s1gate"), "--"],
    ["prediction study, training set and test set", v("design.block.s1train") + ", " + v("design.block.s1test"), v("s1.nruns")],
    ["policy test, check of the rule code", v("design.block.s2gate"), "--"],
    ["policy test, tuning", v("design.block.s2tune"), v("s2.nruns_tune")],
    ["policy test, scoring", v("design.block.s2score"), v("s2.nruns_score")],
    ["reserved and never used", v("design.block.reserved"), "--"]], align="lrr")

# threshold transfer
table("tab_lead_main", [TAU, "$N^{*}$ (cells)", "$N_{\\mathrm{safe}}$ (cells)", "CI (cells)", "ODE safe threshold (cells)",
                        "derivation (%)", "transfer (%)", "total (%)", "grid step (cells)"],
      [[tau(t), v(f"nstar.{t}"), v(f"main.nsafe.{t}"), v(f"main.nsafe_ci.{t}"), v(f"ode.bnd.{t}"), v(f"lead.der.{t}"),
        v(f"lead.trf.{t}"), v(f"lead.tot.{t}"), v(f"main.step.{t}")] for t in T6], wide=True)
table("tab_recal", [TAU, "$N^{*}_{\\mathrm{recal}}$ (cells)", "$N_{\\mathrm{safe}}$ (cells)", "vs median (%)",
                    "$N_{q05}$ (cells)", "vs $N_{q05}$ (%)"],
      [[tau(t), v(f"recal.nstar.{t}"), v(f"recal.nsafe.{t}"), v(f"recal.dmed.{t}"), v(f"recal.nq05.{t}"),
        v(f"recal.dq05.{t}")] for t in T6])
table("tab_lead_split", [TAU, "ODE safe threshold (cells)", "CTMC safe threshold (cells)", "ABM $N_{\\mathrm{safe}}$ (cells)",
                         "spatial share (%)", "ODE-to-CTMC share, range (%)", "ODE threshold set by"],
      [[tau(t), v(f"odeR.{t}"), v(f"ctmc.{t}"), v(f"recal.nsafe.{t}"), v(f"share.{t}"),
        v(f"band.{t}") if t in T4 else "--", "decision calendar" if t in FLOOR else "regrowth"] for t in T6], align="rrrrrrl",
      wide=True)
table("tab_correction", [TAU, "main-grid $N_{\\mathrm{safe}}$ (cells)", "finer-grid $N_{\\mathrm{safe}}$ (cells)",
                         "finer-grid step (cells)", "final tuning interval (cells)", "shift (cells)"],
      [[tau(t), v(f"main.nsafe.{t}"), v(f"corr.s1.{t}"), v(f"corr.s1step.{t}"), ci(f"corr.tlo.{t}", f"corr.thi.{t}"),
        v(f"corr.dnsafe.{t}")] for t in T4], wide=True)
table("tab_corrected", [TAU, "transfer, main grid (%)", "transfer, corrected (%)", "spatial share, main grid (%)",
                        "spatial share, corrected (%)", "vs $N_{q05}$, main grid (%)", "vs $N_{q05}$, corrected (%)"],
      [[tau(t), v(f"lead.trf.{t}"), v(f"corr.trf.{t}"), v(f"share.{t}"), v(f"corr.share.{t}"),
        v(f"recal.dq05.{t}"), v(f"corr.dq05.{t}")] for t in T4], wide=True)

# prediction study
def dl(pre, t):
    """Relative log-loss reduction with its CI."""
    return v(f"{pre}.d.{t}") + " " + ci(f"{pre}.lo.{t}", f"{pre}.hi.{t}")
def wb(pre, t):
    """Within-bin ΔAUC with its CI and number of bins."""
    return v(f"{pre}.w.{t}") + " " + ci(f"{pre}.wlo.{t}", f"{pre}.whi.{t}") + " (" + v(f"{pre}.nb.{t}") + ")"


table("tab_qact", [TAU, "hold decisions, test set", "events, test / training", "Δ log-loss (%) [CI]",
                   "ΔAUC within size bins [CI] (bins)", "Δ log-loss, $C$ in both models (%)"],
      [[tau(t), v(f"s1.hold.{t}"), f"{v(f's1.ev.{t}')} / {v(f's1.evtr.{t}')}", dl("qact", t), wb("qact", t), v(f"qactc.d.{t}")] for t in T4]
      + [["pooled", "", "", dl("qact", "pool"), "", v("qactc.d.pool")]], align="rrrlll", wide=True)
table("tab_qsp", [TAU, "geometry, Δ log-loss (%) [CI]", "geometry, ΔAUC within size bins [CI] (bins)",
                  "image geometry, $R$ in both (%)", "image geometry, no $R$ (%)", "image geometry over history (%)"],
      [[tau(t), dl("qsp", t), wb("qsp", t), v(f"qspi.d.{t}"), v(f"qspo.d.{t}"), v(f"clin.d.{t}")] for t in T4]
      + [["pooled", dl("qsp", "pool"), "", v("qspi.d.pool"), v("qspo.d.pool"), v("clin.d.pool")]], align="rllrrr", wide=True)
table("tab_scale", [TAU, "event base rate (%)", "base-rate model (nats)", "size model (nats)", "history model (nats)",
                    "state model (nats)", "geometry model (nats)"],
      [[tau(t), v(f"abs.base.{t}"), v(f"abs.null.{t}"), v(f"abs.a0.{t}"), v(f"abs.a1.{t}"), v(f"abs.s0.{t}"), v(f"abs.s1.{t}")]
       for t in list(T4) + ["pool"]], wide=True)
table("tab_rcond", [TAU, "history over size, $R$ in both (%)", "$R$ over size (%)", "growth-rate model, share of history gain (%)"],
      [[tau(t), v(f"rd1.d.{t}"), v(f"rd0.d.{t}"), v(f"vel.rec.{t}")] for t in T4]
      + [["pooled", v("rd1.d.pool") + " " + ci("rd1.lo", "rd1.hi"), v("rd0.d.pool") + " " + ci("rd0.lo", "rd0.hi"),
          v("vel.rec.pool") + " " + ci("vel.rec.lo", "vel.rec.hi")]], align="rlll", wide=True)

# policy test
table("tab_tune", ["rule", "tuned quantity"] + ["τ = " + v(f"design.tau.{t}") + " d" for t in T4],
      [["size-only", "threshold (cells)"] + [f"{v(f'tune.sel.SIZE.{t}')} ({v(f'tune.p.SIZE.{t}')})" for t in T4]]
      + [[name, "$\\log_{10} p^{*}$"] + [f"{v(f'tune.sel.{a}.{t}')} ({v(f'tune.p.{a}.{t}')})" for t in T4] for a, name in ARMS],
      align="llllll")
table("tab_policy_primary", [TAU, "$p_{\\mathrm{prem}}$", "relative Δ$f$ (%) [CI]", "ΔRMST (d) [CI]",
                             "Δ$p_{\\mathrm{prem}}$ [CI]", "reading"],
      [[tau(t), v(f"sc.p.CLIN.{t}"),
        v(f"sc.rel.CLIN.{t}"),
        v(f"sc.dr.CLIN.{t}") + " " + ci(f"sc.drlo.CLIN.{t}", f"sc.drhi.CLIN.{t}"),
        v(f"sc.dp.CLIN.{t}") + " " + ci(f"sc.dplo.CLIN.{t}", f"sc.dphi.CLIN.{t}"), words(f"sc.read.CLIN.{t}")] for t in T4]
      + [["pooled", "", v("pool.rel.CLIN") + " " + ci("pool.rello.CLIN", "pool.relhi.CLIN"),
          v("pool.dr.CLIN") + " " + ci("pool.drlo.CLIN", "pool.drhi.CLIN"), v("pool.dp.CLIN") + " " + ci("pool.dplo.CLIN", "pool.dphi.CLIN"),
          words("pool.read.CLIN")]], align="rrllll", wide=True)
table("tab_policy_secondary", ["rule", TAU, "$p_{\\mathrm{prem}}$", "relative Δ on-drug fraction (%)", "ΔRMST (d)", "Δ$p_{\\mathrm{prem}}$ [CI]", "reading"],
      [[name, tau(t), v(f"sc.p.{a}.{t}"), v(f"sc.rel.{a}.{t}"), v(f"sc.dr.{a}.{t}"),
        v(f"sc.dp.{a}.{t}") + " " + ci(f"sc.dplo.{a}.{t}", f"sc.dphi.{a}.{t}"), words(f"sc.read.{a}.{t}")]
       for a, name in ARMS[1:] for t in T4], align="lrrrrll", wide=True)
table("tab_policy_pooled", ["rule", "relative Δ on-drug fraction (%) [CI]", "ΔRMST (d) [CI]", "Δ$p_{\\mathrm{prem}}$ [CI]", "reading",
                            "rule for a pooled result"],
      [[name, v(f"pool.rel.{a}") + " " + ci(f"pool.rello.{a}", f"pool.relhi.{a}"),
        v(f"pool.dr.{a}") + " " + ci(f"pool.drlo.{a}", f"pool.drhi.{a}"), v(f"pool.dp.{a}") + " " + ci(f"pool.dplo.{a}", f"pool.dphi.{a}"),
        words(f"pool.read.{a}"), words(f"pool.stands.{a}")] for a, name in ARMS], align="llllll", wide=True)


def colspec(align):
    """LaTeX column specification; 'p' is a ragged-right column that wraps at half the line width (needs the array package)."""
    return align.replace("p", ">{\\raggedright\\arraybackslash}p{0.5\\linewidth}")


def latex(t):
    """LaTeX source of one table (booktabs, scaled to the text width)."""
    header = [RT.to_tex(h) for h in t["header"]]
    if t["wide"]:   # two-line headers: break before a parenthesis or bracket
        header = [f"\\makecell{{{h.replace(' (', '\\\\(').replace(' [', '\\\\{}[')}}}" if (" (" in h or " [" in h) else h for h in header]
    body = [f"\\begin{{tabular}}{{{colspec(t['align'])}}}", "\\toprule", " & ".join(header) + " \\\\", "\\midrule"]
    body += [" & ".join(RT.to_tex(c) for c in r) + " \\\\" for r in t["rows"]] + ["\\bottomrule", "\\end{tabular}"]
    body = (["\\footnotesize\\setlength{\\tabcolsep}{3pt}"] if t["wide"] else []) + ["\\begin{adjustbox}{max width=\\textwidth}"] + body + ["\\end{adjustbox}"]
    return "% generated by run.py tables from data/results/ledger.json; do not edit\n" + "\n".join(body) + "\n"


def main():
    """Write every table to TABLES and the list of ledger ids they use."""
    TABLES.mkdir(parents=True, exist_ok=True)
    used = set()
    for name, t in TABLES_DEF.items():
        (TABLES / f"{name}.tex").write_text(latex(t))
        for s in t["header"] + [c for r in t["rows"] for c in r]:
            used.update(RT.ledger_ids(s))
    json.dump(sorted(used), open(TABLES / "ledger_ids_used.json", "w"), indent=0)
    print(f"tables written: {len(TABLES_DEF)}; ledger ids used: {len(used)}")
