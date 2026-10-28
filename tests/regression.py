"""Regression checks: this code reproduces the stored results and the research code's outputs exactly.

  python tests/regression.py --quick     about 2 minutes; uses the results in data/results
  python tests/regression.py --full      about 6 minutes; adds refits, rescoring and every ODE/CTMC step in full

Each of seven experiments reruns a subset of at least 200 runs or replicates and compares (1) every field with
the stored raw results and (2) a SHA-256 of all outputs, including the simulator's final random-number state,
with the digest of the same subset run by the original research code (tests/expected/subsets.json).
Analysis checks rerun a step on a scratch copy of data/results and compare its output with the stored file.
Exit status 1 if any check fails.
"""
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "data" / "results"
PY = sys.executable
EXPECTED = json.load(open(ROOT / "tests" / "expected" / "subsets.json"))
SIM_EXPERIMENTS = ["main_abm", "main_arms", "recal_map", "recal_sweep", "recal_ctmc", "statesuff_log", "policy"]


def worker(exp):
    r = subprocess.run([PY, str(ROOT / "tests" / "_subset_worker.py"), exp],
                       capture_output=True, text=True, cwd=ROOT)
    if r.returncode: raise RuntimeError(r.stderr[-2000:])
    return json.loads(r.stdout.strip().splitlines()[-1])


def f(x): return float(x)


# comparisons with the stored raw results
def archive_main_abm(out):
    z = np.load(RES / "p_t_main_raw.npz"); rec, led = z["rec"], z["led"]; bad = 0
    for (row, _state) in out:
        pol, tau, m, s, prem, tp, burden, nmax, l = row
        k = np.flatnonzero((rec["pol"] == pol) & (rec["tau"] == tau) & np.isclose(rec["m"], f(m)) & (rec["seed"] == s))
        assert len(k) == 1
        r = rec[k[0]]
        bad += int((r["prem"], r["t_prog"], r["burden"], r["nmax"]) != (prem, f(tp), f(burden), f(nmax)))
        bad += int(not np.array_equal(led[k[0]], np.array([f(x) for x in l])))
    return len(out), bad


def archive_main_arms(out):
    cells = json.load(open(RES / "p_t_main_arms.json"))["cells"]; bad = 0; n = 0
    for key, v in out.items():
        a = cells[key]
        if key.startswith("ode"):
            n += 1; bad += int(a["prem"] != v["prem"] or not (f(v["t_prog"]) == a["t_prog"] or (np.isnan(f(v["t_prog"])) and np.isnan(a["t_prog"]))))
        else:
            n += len(a["prem"]); bad += sum(x != y for x, y in zip(a["prem"], v["prem"]))
            bad += sum(x != f(y) for x, y in zip(a["t_prog"], v["t_prog"]))
    return n, bad


def archive_recal_map(out):
    z = np.load(RES / "p_t_recal_map.npz"); idx = {int(n): i for i, n in enumerate(z["map_n"])}; bad = 0
    for (n_start, s, t, tr), _state in out:
        i = idx[n_start]; a = z["transit"][i, s]
        bad += int(not (a == f(t) or (np.isnan(a) and np.isnan(f(t)))))
        trace = np.array([f(x) for x in tr], dtype=np.float32); ar = z["traces"][i, s, :len(trace)]
        bad += int(not np.array_equal(ar, trace)) + int(not np.all(np.isnan(z["traces"][i, s, len(trace):])))
    return len(out), bad


def archive_recal_sweep(out):
    rec = np.load(RES / "p_t_recal_sweep.npz")["rec"]; bad = 0
    for (tau, m, s, prem, tp, nmax), _state in out:
        k = np.flatnonzero((rec["tau"] == tau) & np.isclose(rec["m"], f(m)) & (rec["seed"] == s)); assert len(k) == 1
        r = rec[k[0]]; bad += int((r["prem"], r["t_prog"], r["nmax"]) != (prem, f(tp), f(nmax)))
    return len(out), bad


def archive_recal_ctmc(out):
    a = json.load(open(RES / "p_t_recal_ctmc.json"))["cells"]["28|0.80"]
    b = json.load(open(RES / "p_t_recal_ctmc_fine.json"))["cells"]["7|0.9750"]
    c = json.load(open(RES / "p_t_recal_ctmc_fine2.json"))["cells"]["56|1.0250"]
    bad = sum(x != y for x, y in zip(a["ctmc_prem"], out[0]["prem"])) + sum(x != f(y) for x, y in zip(a["ctmc_tprog"], out[0]["t_prog"]))
    bad += sum(x != y for x, y in zip(b["prem"], out[1]["prem"])) + sum(x != y for x, y in zip(c["prem"], out[2]["prem"]))
    return 3 * 500, bad


def archive_statesuff_log(out):
    bad = 0; cache = {}
    for (run, eps, _r0, _r1) in out:
        tau, thr, s = run[0], run[1], run[2]
        role = "train" if (tau, thr) == (28, 2810) else "test"
        fn = RES / "p_s_cells" / role / f"{tau}_{thr}__{(s // 25) * 25:03d}.npz"
        if fn not in cache: cache[fn] = np.load(fn)
        z = cache[fn]; r = z["runs"][z["runs"]["s"] == s][0]
        bad += int(tuple(r.tolist()) != (tau, thr, s, run[3], f(run[4]), f(run[5]), f(run[6]), f(run[7])))
        e = z["eps"][z["eps"]["s"] == s]
        mine = [tuple(f(x) if isinstance(x, str) else x for x in row) for row in eps]
        bad += int(len(e) != len(mine) or any(tuple(a.tolist()) != b for a, b in zip(e, mine)))
    return len(out), bad


def archive_policy(out):
    sel = json.load(open(RES / "p_s_policy" / "tune_selected.json"))["selected"]; bad = 0; cache = {}
    for i, r in enumerate(out):
        arm, tau = ("CLIN", 7) if i < 100 else ("SIZE", 28); s = i % 100
        setting = sel[f"{arm}|{tau}"]["selected"]
        sstr = f"{round(float(setting), 4):.4f}" if arm == "SIZE" else f"{round(float(setting), 6):.6f}"
        fn = RES / "p_s_policy" / "score" / f"{arm}_{tau}_{sstr}__{(s // 25) * 25:03d}.npz"
        if fn not in cache: cache[fn] = np.load(fn)
        z = cache[fn]; a = z["runs"][z["runs"]["s"] == s][0]
        tp = -1.0 if r["t_prog"] == "nan" else f(r["t_prog"])
        bad += int((a["prem"], a["t_prog"], a["burden"], a["N_end"], a["N_max"], a["drug"], a["T"], a["caps"]) !=
                   (r["prem"], tp, f(r["burden"]), f(r["N_end"]), f(r["N_max"]), r["drug"], r["T"], r["caps"]))
        e = z["eps"][z["eps"]["s"] == s]
        bad += int(len(e) != len(r["eps"]) or any((x["e"], x["N"], x["on"], x["cap"], x["y"], x["trunc"]) != (y["e"], f(y["N"]), y["on"], y["cap"], y["y"], y["trunc"])
                                                    or not (x["p"] == f(y["p"]) or (np.isnan(x["p"]) and np.isnan(f(y["p"]))))
                                                    for x, y in zip(e, r["eps"])))
    return len(out), bad


ARCHIVE = dict(main_abm=archive_main_abm, main_arms=archive_main_arms, recal_map=archive_recal_map,
               recal_sweep=archive_recal_sweep, recal_ctmc=archive_recal_ctmc, statesuff_log=archive_statesuff_log,
               policy=archive_policy)


# analysis steps rerun on a scratch copy of data/results
def scratch(files):
    d = Path(tempfile.mkdtemp(prefix="tt_regress_"))
    for rel in files:
        src = RES / rel; dst = d / rel; dst.parent.mkdir(parents=True, exist_ok=True)
        (shutil.copytree if src.is_dir() else shutil.copy2)(src, dst)
    return d


def run_script(args, results, extra_env=None):
    env = dict(os.environ, TT_RESULTS=str(results), **(extra_env or {}))
    r = subprocess.run([PY] + args, capture_output=True, text=True, cwd=ROOT, env=env)
    if r.returncode: raise RuntimeError(" ".join(args) + "\n" + r.stderr[-3000:])
    return r.stdout


def same_json(a, b, drop=()):
    A, B = json.load(open(a)), json.load(open(b))
    for k in drop: A.pop(k, None); B.pop(k, None)
    return A == B


def a_scale():
    d = scratch(["p_s_models.pkl", "p_s_fit.json", "p_s_results.json", "p_s_cells/test"])
    run_script(["run.py", "predict", "scale"], d)
    return (d / "p_s_absolute_scale.json").read_bytes() == (RES / "p_s_absolute_scale.json").read_bytes()


def a_policy_analyse():
    d = scratch(["p_s_models.pkl", "p_s_fit.json", "p_s_vel_model.pkl", "p_s_vel_fit.json", "p_s_policy/score", "p_s_policy/tune_selected.json"])
    run_script(["run.py", "policy", "analyse"], d)
    return (d / "p_s_policy_results.json").read_bytes() == (RES / "p_s_policy_results.json").read_bytes()


def a_recal_fit():
    """N*_recal must be identical. The ODE boundary depends on the least-squares refit of (c, phi_bulk), which
    on the test machine differs from the stored refit by about 1e-7 relative (the research code gives the same
    value as this code); so (1) the boundary from the stored refit constants must equal the stored boundary
    exactly, and (2) a fresh refit's boundary must agree within 1e-6 relative."""
    d = scratch(["p_t_recal_map.npz", "p_t_recal_sweep.npz"])
    run_script(["run.py", "recal", "fit"], d)
    ok = (d / "p_t_recal_nstar.json").read_bytes() == (RES / "p_t_recal_nstar.json").read_bytes()
    old, new = json.load(open(RES / "p_t_recal_odeb.json")), json.load(open(d / "p_t_recal_odeb.json"))
    ok &= all(abs(new[k] / old[k] - 1) < 1e-6 for k in old)
    code = ("import json, numpy as np, sys; sys.path.insert(0, 'src'); from threshold_transfer.model.ode import ode_boundary; "
            "from threshold_transfer.paths import RESULTS; z = np.load(RESULTS / 'p_t_recal_sweep.npz'); "
            "c, pb = float(z['c_refit']), float(z['phi_refit']); "
            "print(json.dumps({str(t): ode_boundary(t, c, pb)[0] for t in (7, 14, 28, 42, 56, 84)}))")
    exact = json.loads(run_script(["-c", code], d).strip().splitlines()[-1])
    return ok and exact == old


def a_vel():
    d = scratch(["p_s_models.pkl", "p_s_fit.json", "p_s_cells/train", "p_s_cells/test"])
    run_script(["run.py", "predict", "vel-fit"], d)
    ok = json.load(open(d / "p_s_vel_fit.json"))["hash"] == json.load(open(RES / "p_s_vel_fit.json"))["hash"]
    run_script(["run.py", "predict", "vel-score"], d)
    return ok and (d / "p_s_vel_prediction.json").read_bytes() == (RES / "p_s_vel_prediction.json").read_bytes()


def a_rdiag():
    d = scratch(["p_s_models.pkl", "p_s_fit.json", "p_s_cells/train", "p_s_cells/test"])
    run_script(["run.py", "predict", "rdiag"], d)
    ok = same_json(d / "p_s_rdiag_fit.json", RES / "p_s_rdiag_fit.json", drop=("t",))
    return ok and (d / "p_s_rdiag.json").read_bytes() == (RES / "p_s_rdiag.json").read_bytes()


def a_fit():
    d = scratch(["p_s_cells/train"])
    run_script(["run.py", "predict", "fit"], d)
    return same_json(d / "p_s_fit.json", RES / "p_s_fit.json", drop=("t",))


def a_score():
    d = scratch(["p_s_models.pkl", "p_s_fit.json", "p_s_cells/train", "p_s_cells/test"])
    run_script(["run.py", "predict", "score"], d)
    return (d / "p_s_results.json").read_bytes() == (RES / "p_s_results.json").read_bytes()


def a_reports():
    """Printed reports of the main run and the recalibration control, against the research code's output."""
    out = {}
    for name, args in (("transfer_report", ["run.py", "transfer", "report"]),
                       ("recal_report", ["run.py", "recal", "report"])):
        out[name] = hashlib.sha256(run_script(args, RES).encode()).hexdigest() == EXPECTED["stdout"][name]
    return all(out.values())


def a_ledger_figures_tables():
    d = Path(tempfile.mkdtemp(prefix="tt_out_"))
    env = dict(TT_FIGURES=str(d / "figures"), TT_TABLES=str(d / "tables"))
    led = scratch([]); shutil.rmtree(led); shutil.copytree(RES, led, ignore=shutil.ignore_patterns("ledger.json"))
    run_script(["run.py", "ledger"], led, env)
    a = json.load(open(led / "ledger.json")); ok = not a["counts"]["failing"]
    ok &= hashlib.sha256(json.dumps({k: [e["kind"], e["text"], e["value"]] for k, e in a["entries"].items()},
                                    sort_keys=True).encode()).hexdigest() == EXPECTED["ledger_values_sha256"]
    run_script(["run.py", "figures"], led, env)
    for n in ("F1", "F2", "F3", "F4", "F5"):
        for ext in ("pdf", "png"):
            ok &= (d / "figures" / f"{n}.{ext}").read_bytes() == (ROOT / "figures" / f"{n}.{ext}").read_bytes()
    run_script(["run.py", "tables"], led, env)
    for p in (ROOT / "tables").glob("*"):
        ok &= (d / "tables" / p.name).read_bytes() == p.read_bytes()
    return ok


def a_ctmc_full():
    """Every ODE/CTMC step in full (arms, recalibrated CTMC, fine grid rounds 1 and 2) from the stored
    growth constants: each output's cells must equal the stored file's."""
    d = scratch(["p_t_nstar.json", "p_t_recal_sweep.npz"])
    for args in (["transfer", "arms"], ["recal", "ctmc"], ["recal", "ctmc-fine", "--round", "1"], ["recal", "ctmc-fine", "--round", "2"]):
        run_script(["run.py"] + args, d)
    ok = (d / "p_t_main_arms.json").read_bytes() == (RES / "p_t_main_arms.json").read_bytes()
    for n in ("p_t_recal_ctmc.json", "p_t_recal_ctmc_fine.json", "p_t_recal_ctmc_fine2.json"):
        a, b = json.load(open(d / n)), json.load(open(RES / n))
        ok &= all(a[k] == b[k] for k in b if k != "wall_s")
    return ok


def a_ode_boundaries():
    """The ODE safe thresholds under the original growth law, which the ledger quotes from the dated record."""
    code = ("import json, sys; sys.path.insert(0, 'src'); from threshold_transfer.model.ode import ode_boundary; "
            "print(json.dumps({t: round(ode_boundary(t)[0], 1) for t in (7, 14, 28, 42, 56, 84)}))")
    got = json.loads(run_script(["-c", code], RES).strip().splitlines()[-1])
    filed = {"7": 3454.0, "14": 3299.4, "28": 3001.0, "42": 2787.7, "56": 2831.6, "84": 2566.2}
    return got == filed


def a_transit_map_full():
    """The full drug-free transit map (3,000 runs): every array equal to the stored one."""
    d = scratch([])
    run_script(["run.py", "recal", "map"], d)
    a, b = np.load(d / "p_t_recal_map.npz"), np.load(RES / "p_t_recal_map.npz")
    return sorted(a.files) == sorted(b.files) and all(np.array_equal(a[k], b[k], equal_nan=a[k].dtype.kind == "f") for k in a.files)


QUICK_ANALYSES = [("ledger, figures and tables", a_ledger_figures_tables), ("absolute scale", a_scale),
                  ("policy-test analysis", a_policy_analyse), ("main-run and recalibration reports", a_reports),
                  ("R-VEL fit and score", a_vel)]
FULL_ANALYSES = [("recalibration fit", a_recal_fit), ("R-conditioning fit and score", a_rdiag),
                 ("Step 1 fit (model hashes)", a_fit), ("Step 1 score", a_score),
                 ("ODE and CTMC steps in full", a_ctmc_full), ("ODE safe thresholds, original law", a_ode_boundaries),
                 ("transit map in full", a_transit_map_full)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(); g.add_argument("--quick", action="store_true"); g.add_argument("--full", action="store_true")
    ap.add_argument("--only", help="comma-separated subset of checks (names as printed)")
    a = ap.parse_args()
    only = set(a.only.split(",")) if a.only else None
    fails = 0; t_all = time.time(); report = []
    for exp in SIM_EXPERIMENTS:
        if only and exp not in only: continue
        t0 = time.time(); r = worker(exp)
        n, bad = ARCHIVE[exp](r["out"]); dig = r["sha256"] == EXPECTED["subsets"][exp]
        line = f"{exp:15s} {n:5d} runs  archive mismatches {bad}  research digest {'equal' if dig else 'DIFFERENT'}"
        ok = bad == 0 and dig
        fails += not ok; line += f"  {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)"; print(line, flush=True); report.append(line)
    for name, fn in QUICK_ANALYSES + (FULL_ANALYSES if a.full else []):
        if only and name not in only: continue
        t0 = time.time()
        try: ok = fn()
        except Exception as e: ok = False; print(e)
        fails += not ok; line = f"{name:45s} {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)"; print(line, flush=True); report.append(line)
    print(f"\n{'ALL PASS' if not fails else f'{fails} FAILED'}  total {time.time() - t_all:.0f} s")
    return fails


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
