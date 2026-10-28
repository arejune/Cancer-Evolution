"""Run any step of the paper. Steps that need no simulation take seconds; see the README for run times.

  python run.py ledger                      recompute every number in the paper (data/results/ledger.json)
  python run.py figures                     figures/F1-F5
  python run.py tables                      tables/*.tex (run `ledger` first)

  python run.py transfer abm|arms|report    main run: agent-based model, ODE and CTMC arms, printed report
  python run.py recal map|fit|sweep|ctmc|ctmc-fine|report
                                            recalibration to the agent-based model's regrowth
  python run.py predict gate|run|fit|score|scale|rdiag|vel-fit|vel-score|dryrun
                                            prediction study
  python run.py policy gate|tune|score|analyse
                                            policy test

`transfer abm`, `predict run`, `policy tune` and `policy score` save their output as they go and skip
files that already exist, so they can be stopped and restarted. `predict run`, `policy tune` and
`policy score` take a time budget in seconds (default: no limit).
"""
import argparse, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

STEPS = {"ledger": None, "figures": None, "tables": None,
         "transfer": ["abm", "arms", "report"],
         "recal": ["map", "fit", "sweep", "ctmc", "ctmc-fine", "report"],
         "predict": ["gate", "run", "fit", "score", "scale", "rdiag", "vel-fit", "vel-score", "dryrun"],
         "policy": ["gate", "tune", "score", "analyse"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("experiment", choices=STEPS)
    ap.add_argument("step", nargs="?")
    ap.add_argument("--budget", type=float, default=float("inf"), help="seconds before stopping (run, tune, score)")
    ap.add_argument("--chunk", type=int, default=25, help="seeds per output file (transfer abm)")
    ap.add_argument("--round", type=int, choices=[1, 2], default=2, help="fine CTMC grid round (recal ctmc-fine)")
    ap.add_argument("--out", help="output path (predict scale)")
    a = ap.parse_args()
    if STEPS[a.experiment] and a.step not in STEPS[a.experiment]:
        ap.error(f"{a.experiment} needs one of: {', '.join(STEPS[a.experiment])}")

    if a.experiment == "ledger":
        from threshold_transfer.report import ledger
        sys.exit(1 if ledger.run()["failing"] else 0)
    if a.experiment == "figures":
        from threshold_transfer.report import figures
        return figures.main()
    if a.experiment == "tables":
        from threshold_transfer.report import tables
        return tables.main()

    if a.experiment == "transfer":
        if a.step == "abm":
            from threshold_transfer.experiments import transfer_main as m; m.main(CHUNK=a.chunk)
        elif a.step == "arms":
            from threshold_transfer.experiments import transfer_arms as m; m.main()
        else:
            from threshold_transfer.analysis import transfer as m; m.main()

    elif a.experiment == "recal":
        if a.step == "map":
            from threshold_transfer.experiments import recal_map as m; m.main()
        elif a.step in ("fit", "sweep"):
            from threshold_transfer.experiments import recal_sweep as m; m.main(no_sweep=(a.step == "fit"))
        elif a.step == "ctmc":
            from threshold_transfer.experiments import recal_ctmc as m; m.main()
        elif a.step == "ctmc-fine":
            from threshold_transfer.experiments import recal_ctmc_fine as m; m.main(round_=a.round)
        else:
            from threshold_transfer.analysis import recal as m; m.main()

    elif a.experiment == "predict":
        if a.step == "gate":
            from threshold_transfer.experiments import statesuff_log as m; m.gate()
        elif a.step == "run":
            from threshold_transfer.experiments import statesuff_log as m; m.supervise(a.budget)
        elif a.step in ("fit", "score", "dryrun"):
            from threshold_transfer.analysis import statesuff as m; m.main(a.step)
        elif a.step == "scale":
            from threshold_transfer.analysis import scale as m; m.main(a.out)
        elif a.step == "rdiag":
            from threshold_transfer.analysis import rdiag as m; m.main()
        else:
            from threshold_transfer.analysis import vel as m; m.main(a.step.split("-")[1])

    elif a.experiment == "policy":
        from threshold_transfer.experiments import policy as m
        {"gate": m.gate, "tune": lambda: m.tune(a.budget), "score": lambda: m.score(a.budget),
         "analyse": m.analyse}[a.step]()


if __name__ == "__main__":
    main()
