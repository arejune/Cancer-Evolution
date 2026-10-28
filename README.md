# Threshold transfer in a spatial tumour model

Code and results for the paper *What a size-only threshold for holding cancer treatment loses in a spatial
tumour model* (Arjun Grandhe). Every number, table and figure in the paper can be recomputed from this
repository, which includes the raw results (66 MB).

The paper takes a treatment threshold derived in an ordinary differential equation (ODE) model (Gallagher
et al. 2025) and measures how safe it is in a spatial agent-based model (ABM) of a tumour. It then tests
whether the history of earlier measurements and the shape of the tumour predict progression during a
treatment hold, and whether a hold rule that uses them does better than a rule based on size alone.

## Install

Python 3.13 with the package versions in `requirements.txt`:

    python3.13 -m venv .venv
    source .venv/bin/activate          # Windows: .venv\Scripts\activate
    pip install -r requirements.txt

or with conda: `conda env create -f environment.yml`, then `conda activate threshold-transfer`.

The raw results are in `data/results`. `python data/check_data.py` checks each of the 1,607 raw files
against its MD5 in `data/MD5SUMS`.

## Check that it works (about 2 minutes)

    python tests/regression.py --quick

This reruns part of every experiment (at least 200 runs each) and checks each run against the stored
results and against a digest of the same runs made by the original research code. It then rebuilds the
number ledger, the figures and the tables and checks they match the files here. `--full` adds the slower
checks (about 6 minutes). `pytest tests` runs the same checks plus unit tests.

The stored results were made on macOS (Apple M3) with the conda packages. With the PyPI packages on the
same machine, every simulated run has the same outcome and every number in the paper is the same; predicted
risks differ by up to 1 part in 10^14, model fits at the solver's tolerance, and the figure files in font
rendering. Other operating systems were not tested; the quick check shows what differs.

## Reproduce the paper from the stored results (seconds)

    python run.py ledger       # every number in the paper -> data/results/ledger.json (833 entries)
    python run.py figures      # figures/F1-F5
    python run.py tables       # tables/*.tex

Each ledger entry records the file a number comes from and how it is computed. Numbers that exist only as
text in a dated record (for example a bootstrap interval whose random state was not saved) are checked
word for word against that record in `config/filings`.

## Rerun the simulations (about 85 minutes)

`python tests/full_rerun.py OUT_DIR` runs every step below into an empty folder and compares every output
with the stored results. The steps can also be run one at a time into `data/results`. Times are from an
8-core laptop.

| step | command | first seed | time |
|---|---|---:|---:|
| main run, ABM | `python run.py transfer abm` | 1,200,000 | 11 min |
| main run, ODE and CTMC | `python run.py transfer arms` | 1,200,000 | 31 s |
| recalibration: drug-free transit times | `python run.py recal map` | 1,800,000 | 18 s |
| recalibration: refit and threshold sweep | `python run.py recal sweep` | 1,800,500 | 9 min |
| recalibration: N* and ODE threshold files | `python run.py recal fit` | | 19 s |
| CTMC under the recalibrated law | `python run.py recal ctmc` | 1,900,000 | 13 s |
| CTMC on the finer grid | `python run.py recal ctmc-fine --round 1`, then `--round 2` | 2,000,000 | 27 s |
| prediction study: code check | `python run.py predict gate` | 2,300,000 | 14 s |
| prediction study: training and test runs | `python run.py predict run` | 2,100,000 and 2,200,000 | 10 min |
| prediction models | `python run.py predict fit`, then `score`, `scale`, `rdiag`, `vel-fit`, `vel-score` | | 2 min |
| policy test: code check | `python run.py policy gate` | 2,600,000 | 14 s |
| policy test: tuning | `python run.py policy tune` | 2,400,000 | 39 min |
| policy test: scoring | `python run.py policy score` | 2,500,000 | 10 min |
| policy test: analysis | `python run.py policy analyse` | | 2 s |

The long steps save their output as they go and skip files that exist, so they can be stopped and
restarted (`--budget SECONDS` stops them after a set time). The pieces written by `transfer abm`
(`p_t_cells/`) and `policy tune` (`p_s_policy/tune/`) are not stored here: the first are merged into
`p_t_main_raw.npz` and the second are summarised in `p_s_policy/tune_selected.json`, so a rerun of either
step starts from the beginning. The fitting steps refuse to overwrite fitted
models; delete the `.pkl` files to refit. `python run.py transfer report` and `python run.py recal report`
print the per-cell results behind the first part of the paper.

A rerun from scratch on the test machine reproduced 1,617 of the 1,622 stored output files exactly. The exception
is the least-squares refit of the growth law in `recal sweep`, which now gives constants that differ from
the stored ones by about 1 part in 10^7 (the research code gives the same new values). That small change
gives 5 of the 177 CTMC cells new random draws, which moves the CTMC safe threshold at tau = 28 from 3,321.6
to 3,323.9 cells and the spatial share there from 85.0% to 85.4%. Every other number is reproduced exactly,
and with the stored constants the CTMC steps match exactly (`tests/regression.py --full` checks this).

## Random-number seeds

Each experiment uses its own range of seeds (table above), and run s of an experiment uses the first seed
plus s, so every threshold in a grid starts from the same tumours. A rerun therefore gives the same runs bit
for bit. Bootstrap intervals use the fixed seeds 20260822, 20260826 and 20260923. The analysis plan for each
experiment was written and dated before its simulations were run (`config/preregistrations`); analyses added
later are marked as such in the paper.

## Files

| file | contents |
|---|---|
| `run.py` | command line for every step |
| `src/threshold_transfer/model/simulator.py` | the agent-based model |
| `src/threshold_transfer/model/ode.py` | shared constants, growth law, N*(tau), and the deterministic ODE |
| `src/threshold_transfer/model/ctmc.py` | the well-mixed stochastic model (CTMC) |
| `src/threshold_transfer/threshold/construction.py` | initial tumour and the size-only policy |
| `src/threshold_transfer/threshold/grid.py` | intervals, tolerance, seeds and threshold grid of the main run |
| `src/threshold_transfer/experiments/transfer_main.py` | main run in the ABM |
| `src/threshold_transfer/experiments/transfer_arms.py` | main run in the ODE and CTMC |
| `src/threshold_transfer/experiments/recal_map.py` | drug-free transit times of the ABM |
| `src/threshold_transfer/experiments/recal_sweep.py` | recalibrated threshold, ODE threshold and ABM sweep |
| `src/threshold_transfer/experiments/recal_ctmc.py` | CTMC under the recalibrated law |
| `src/threshold_transfer/experiments/recal_ctmc_fine.py` | the same on a finer grid |
| `src/threshold_transfer/experiments/statesuff_log.py` | prediction study runs, with the state recorded at each decision |
| `src/threshold_transfer/experiments/policy.py` | policy test: rules, tuning, scoring and analysis |
| `src/threshold_transfer/analysis/statesuff.py` | prediction models and their comparison |
| `src/threshold_transfer/analysis/scale.py` | absolute scale of the prediction results |
| `src/threshold_transfer/analysis/rdiag.py` | check of the history gain against the resistant count |
| `src/threshold_transfer/analysis/vel.py` | growth-rate prediction model |
| `src/threshold_transfer/analysis/transfer.py`, `recal.py` | printed reports of the main run and the recalibration |
| `src/threshold_transfer/report/ledger.py` | the number ledger |
| `src/threshold_transfer/report/figures.py` | the figures |
| `src/threshold_transfer/report/tables.py`, `richtext.py` | the tables |
| `src/threshold_transfer/paths.py` | folder locations and the number of worker processes |
| `tests/regression.py` | the quick and full checks |
| `tests/_subset_worker.py` | runs the fixed subsets used by the checks |
| `tests/full_rerun.py` | rerun everything from scratch and compare |
| `tests/test_units.py`, `test_regression.py`, `conftest.py` | the same checks under pytest, plus unit tests |
| `tests/expected/subsets.json` | digests of the research code's outputs |
| `data/results/` | results: raw simulation output, fitted models and summary files |
| `data/check_data.py`, `MD5SUMS` | checksums of the raw results |
| `config/preregistrations/` | the dated analysis plans |
| `config/filings/` | dated result records and excerpts, used by the ledger |
| `figures/`, `tables/` | outputs of `run.py figures` and `run.py tables` |

Names in the code are shorter than in the paper. The prediction study is "statesuff" and the recalibration
"recal". The models A0, A1, S0, S1F and S1I are the size, history, state, geometry and image-geometry models,
and the rules SIZE, CLIN, HIST, FULL and VEL are the size-only, history-and-image, history, full-state and
growth-rate rules. A "gate" is a check of a program against the original code on its own seeds.

| paper | file |
|---|---|
| Figures 1-4, Figure C1 | `figures/F1`-`F4`, `figures/F5` |
| Tables 1-4 | `tab_config`, `tab_models`, `tab_lead_main`, `tab_policy_primary` |
| Table A1 | `tab_seeds` |
| Tables C1-C2 | `tab_correction`, `tab_corrected` |
| Tables D1-D9 | `tab_recal`, `tab_lead_split`, `tab_qact`, `tab_rcond`, `tab_qsp`, `tab_scale`, `tab_tune`, `tab_policy_pooled`, `tab_policy_secondary` |

## Citation

Grandhe A. What a size-only threshold for holding cancer treatment loses in a spatial tumour model.
bioRxiv preprint, 2026. DOI: N/A. Code and data: https://github.com/arejune/Cancer-Evolution. See also `CITATION.cff`.

## Licence

Code: MIT (`LICENSE`). Data, figures, tables and text: CC BY 4.0 (`LICENSE-data`).
