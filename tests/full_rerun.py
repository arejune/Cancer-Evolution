"""Full reproduction from scratch: every simulation and analysis step, timed, into an empty folder, then a
comparison of every output with the archived results in data/results.

  python tests/full_rerun.py OUT_DIR          about 85 minutes on an 8-core laptop

Arrays in .npz files are compared element by element (file bytes differ because numpy writes the current time
into each archive member); fitted models are compared through the SHA-256 of their coefficients; JSON files are
compared field by field, leaving out wall-clock times and the cost projections derived from them.
Writes OUT_DIR/rerun_report.json with the time of each step and the result of each comparison.
"""
import json, os, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "data" / "results"
STEPS = [("main run, ABM", ["transfer", "abm"]),
         ("main run, ODE and CTMC arms", ["transfer", "arms"]),
         ("recalibration, transit map", ["recal", "map"]),
         ("recalibration, fit and sweep", ["recal", "sweep"]),
         ("recalibration, N* and ODE threshold files", ["recal", "fit"]),
         ("CTMC, recalibrated law", ["recal", "ctmc"]),
         ("CTMC, fine grid round 1", ["recal", "ctmc-fine", "--round", "1"]),
         ("CTMC, fine grid round 2", ["recal", "ctmc-fine", "--round", "2"]),
         ("prediction study, gate", ["predict", "gate"]),
         ("prediction study, training and test sets", ["predict", "run"]),
         ("prediction models, fit", ["predict", "fit"]),
         ("prediction models, score", ["predict", "score"]),
         ("absolute scale", ["predict", "scale"]),
         ("resistant-count check", ["predict", "rdiag"]),
         ("growth-rate model, fit", ["predict", "vel-fit"]),
         ("growth-rate model, score", ["predict", "vel-score"]),
         ("policy test, gate", ["policy", "gate"]),
         ("policy test, tuning", ["policy", "tune"]),
         ("policy test, scoring", ["policy", "score"]),
         ("policy test, analysis", ["policy", "analyse"]),
         ("number ledger", ["ledger"])]
TIMING = {"t", "wall_s", "wall_ref_s", "wall_log_s", "ms_per_run_ref", "ms_per_run_log", "overhead_ratio", "projected_wall_min",
          "us_per_day_single_worker", "gate_wall_s", "us_per_day_throughput_7w", "arm_cost_ratio_vs_SIZE",
          "tune_projection_min_central_pilot", "score_projection_min_central_pilot", "score_projection_min_worst_pilot"}


def strip(x):
    if isinstance(x, dict): return {k: strip(v) for k, v in x.items() if k not in TIMING}
    if isinstance(x, list): return [strip(v) for v in x]
    return x


def same_array(a, b):
    """Element-wise equality, with NaN equal to NaN, also inside structured arrays."""
    if a.dtype != b.dtype or a.shape != b.shape: return False
    if a.dtype.names: return all(same_array(a[n], b[n]) for n in a.dtype.names)
    return np.array_equal(a, b, equal_nan=a.dtype.kind == "f")


def compare(out):
    res = {}
    for p in sorted(RES.rglob("*")):
        if not p.is_file() or p.name in ("ledger.json",) or p.suffix == ".jsonl": continue
        rel = p.relative_to(RES); q = out / rel
        if p.name.endswith("_cost.json"): continue                      # wall-clock records only
        if not q.exists(): res[str(rel)] = "missing"; continue
        if p.suffix == ".npz":
            a, b = np.load(p), np.load(q)
            ok = sorted(a.files) == sorted(b.files) and all(same_array(a[k], b[k]) for k in a.files)
        elif p.suffix == ".json":
            ok = strip(json.load(open(p))) == strip(json.load(open(q)))
        elif p.suffix == ".pkl":
            continue                                                       # compared through the hashes in the fit files
        else:
            ok = p.read_bytes() == q.read_bytes()
        res[str(rel)] = "identical" if ok else "DIFFERENT"
    return res


if __name__ == "__main__":
    out = Path(sys.argv[1]).resolve(); out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TT_RESULTS=str(out), TT_FIGURES=str(out / "figures"), TT_TABLES=str(out / "tables"))
    times = {}
    for name, args in STEPS:
        t0 = time.time()
        r = subprocess.run([sys.executable, str(ROOT / "run.py")] + args, cwd=ROOT, env=env,
                           capture_output=True, text=True)
        times[name] = round(time.time() - t0, 1)
        (out / "logs").mkdir(exist_ok=True); (out / "logs" / ("_".join(args) + ".log")).write_text(r.stdout + r.stderr)
        print(f"{name:45s} {times[name]:8.1f} s  exit {r.returncode}", flush=True)
        if r.returncode and name != "number ledger": sys.exit(f"step failed: {name}")
    cmp_ = compare(out)
    bad = {k: v for k, v in cmp_.items() if v != "identical"}
    led_a, led_b = json.load(open(RES / "ledger.json"))["entries"], json.load(open(out / "ledger.json"))["entries"]
    led_diff = sorted(k for k in led_a if led_a[k]["text"] != led_b[k]["text"])
    rep = dict(step_seconds=times, total_minutes=round(sum(times.values()) / 60, 1), files_compared=len(cmp_),
               files_not_identical=bad, ledger_entries_with_different_text=led_diff)
    json.dump(rep, open(out / "rerun_report.json", "w"), indent=1)
    print(json.dumps(rep, indent=1))
