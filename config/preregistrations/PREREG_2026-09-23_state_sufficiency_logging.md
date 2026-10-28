# PRE-REGISTRATION — State sufficiency, Step 1 (logging run)
**Filed 2026-09-23, before any seed of blocks 2100000 / 2200000 / 2300000 was read.
Status: SIGNED OFF 2026-09-23 (all eight S13 items approved) with Amendments A–D below, filed before any seed was read.**

Parent: the state-insufficiency gap named in `research_log.md` §7 RECAL row (2026-08-27), which
filed *"Total burden is not a sufficient state variable for this control problem."* That sentence
rests on a threshold-calibration residual, not on a test of state at matched N. This filing is
the test.

---

## S0. Step 0 — closed, zero compute [VERIFIED]

`ode/results/p_t_main_raw.npz` and all 1,300 shards in `ode/results/p_t_cells/` store per-RUN
scalars only (`pol, tau, m, seed, prem, t_prog, burden, nmax`) plus a whole-run (12,) ledger.
`p_t_pilot.run()` (lines 107–128) reads N at each epoch and discards it. No N, drug state,
composition or geometry exists at epoch level. The replay of block 1200000 was **declined**: the
block is spent and N_safe was tuned on it. No block was consumed by Step 0.

## S1. Seed blocks [VERIFIED clean 2026-09-23]

| block | role | read |
|---|---|---|
| **2100000** + s, s = 0..499 | TRAIN | once |
| **2200000** + s, s = 0..499 | TEST | once |
| **2300000** + s, s = 0..9 | GATE + COST PILOT + analysis-code dry run. Not a registered experiment block; its outcomes feed no estimate. | once |

Verification: a repo-wide scan for every integer in [2,000,000, 2,999,999] (with `_` / `,`
separators) returns only `2000000` (spent). The `22…`/`24…` substring hits are the
`+222_000_000` / `+424_000_000` placement offsets, not blocks. Placement stream for these blocks
is seed + 313,000,000 (inherited from `p_t_pilot.build`), i.e. 315,100,000+ / 315,200,000+ /
315,300,000+, which collides with no simulator seed in use.

CRN: seed index s is the same initial condition in every cell of a block. **The cluster unit
everywhere below is the seed index**, across all cells and all tau.

Step 2 will require its own fresh blocks (see S12). 1700000 stays reserved and untouched.

## S2. Run design

Policy, IC, horizon and progression rule are **exactly** `p_t_pilot.run(policy="AT")` as used in
the main run: N0 = 3001 radius-31 disk, 1% resistant converted by the dedicated placement stream,
no mutation, DOSE = 1.0 of drug A when on, u_B = 0, decision every tau days from day 0
(`on = N >= thr`), progression at N >= 1.2·N0 = 3601.2 checked daily, horizon T_END = 750 d.
n_max = L² = 43,264 is inherited from `p_t_pilot` (closed files untouched).

**Threshold band — absolute N, 5 per tau, 20 cells.** Lower edge ≈ N_safe − 50; upper edge ≈ the
threshold at which main-run p_prem reached ~0.5, so both outcomes occur. Band chosen from
block-1200000 cell summaries (design use of a spent block, disclosed; no outcome of the new
blocks informs it).

| tau | N_safe (main) | thresholds | expected p_prem by cell [UNCERTAIN, linear interp. of main run] | expected events / block |
|---:|---:|---|---|---:|
| 7  | 3114.6 | 3060 3120 3180 3240 3300 | 0.000 0.058 0.152 0.246 0.458 | ~457 |
| 14 | 2946.1 | 2900 2950 3000 3050 3100 | 0.000 0.062 0.215 0.368 0.520 | ~582 |
| 28 | 2719.9 | 2670 2740 2810 2880 2950 | 0.019 0.063 0.153 0.338 0.519 | ~545 |
| 56 | 2718.3 | 2660 2740 2820 2900 2980 | 0.033 0.069 0.138 0.264 0.423 | ~463 |

R = 500 per cell per block. 20 cells × 500 seeds × 2 blocks = **20,000 runs**.
Thresholds are defined at module level in the runner (the macOS-spawn KeyError of the recal run,
research_log §7 RECAL, must not recur).

## S3. Instrumentation — NEW file `ode/p_s_log.py`

The runner imports `p_t_pilot.build`, `DOSE`, `PROG_MUL` and reproduces the `run()` loop
statement-for-statement. The only additions are **reads**: `sim.get_state()`, `sim.grid`,
`sim.C_A`, copied into numpy arrays. No call to `sim._rng`, no call that consumes randomness, no
write to any simulator attribute. `quad2/simulator2.py`, `ode/growth_forms.py`,
`ode/wellmixed_ctmc2.py`, `ode/p_t_pilot.py` are not edited.

**Logged at every decision epoch e (day d_e = e·tau), from the state read that makes the decision:**

| field | definition |
|---|---|
| `N` | `true_total_volume` |
| `S`, `R` | `n_S`, `n_Ra`; `n_Escape` logged and asserted 0 (no mutation) |
| `on` | decision at this epoch (`N >= thr`) |
| `C` | mean of `sim.C_A` over occupied sites (grid != 0) |
| `B` | occupied sites with ≥ 1 empty Moore (8-) neighbour; out-of-lattice counts as empty (convention of `p_t_pilot.geometry`) |
| `B_R` | resistant sites satisfying the `B` criterion |
| `K` | number of 8-connected components of the occupied set (`scipy.ndimage.label`, 3×3 structure) |
| `Rg` | sqrt of mean squared distance of occupied sites from their centroid |
| `Rmax` | max distance of an occupied site from that centroid |
| `d_wall` | min distance from any occupied site to the lattice edge — sanity only, not a feature |

**History features, computed from epoch observations and past decisions only (clinically
available — no daily trace is used):**

| feature | definition |
|---|---|
| `tsf` | days since the most recent ON epoch; if none, days since day 0, with indicator `never_on` |
| `N_prev` | N at epoch e−1; entered as the growth ratio `g = log(N / N_prev)` |
| `n_on` | number of ON epochs before e ("prior firings") |
| `nadir` | min epoch-observed N over epochs strictly after the most recent ON epoch, up to and including e; if never ON, over epochs 0..e |

Per-run scalars are also written (`prem, t_prog, burden, N_end, N_max`), identical in meaning to
the main run, for the gate. Checkpointed in shards of (block, cell, 50 seeds) under a resuming
supervisor; seeds are a pure function of the index, so a restart is not a second read.

## S4. Gate, cost pilot, launch rule

**Gate (block 2300000, seeds 0–9, all 20 cells = 200 runs), VOID if any item fails:**
1. `p_t_pilot.run()` (uninstrumented, the code that produced N_safe) and `p_s_log` run on every
   gate seed/cell: `burden`, `t_prog`, `prem`, `N_end`, `N_max` **exactly equal**, 200/200.
2. RNG state: final `sim._rng.bit_generator.state` of the instrumented loop equals that of a
   non-logging replica of the same loop, seed-for-seed. **Non-vacuity check:** the final state
   must differ from the post-`build` state. If it does not (the njit kernel may not advance the
   Python-visible state), item 2 is reported as UNINFORMATIVE and the gate rests on item 1,
   which remains mandatory.
3. Label identity: for every run, Σ(hold-epoch labels) == `prem` (S5 makes this an identity).
4. `n_Escape == 0` at every epoch; `d_wall ≥ 3` at every epoch (tumour never touches the edge).

Failure of 3 or 4 is a logging bug, not a block read; fix and re-run the gate on the same block
(deterministic, bit-identical, not a second read). Failure of 1 or 2 = VOID; report and stop.

**Cost projection [UNCERTAIN, from the contended basis]:** 20,000 runs, 6.70 M simulated days
(cell-mean days interpolated from the main run, 295–414 d/run), ~473,000 epochs. At the
contended **650 µs/simulated-day** basis: **72.6 min** uninstrumented. At 179.9 ms/run: 60.0 min
(the per-day basis is used as the governing figure because these cells run longer than the
recal sweep's mean). Instrumentation overhead is unmeasured; the precedent (WRITEUP_ASSETS §21)
is a **2.8× overrun** from per-epoch Python recording, so the overhead is measured, not assumed.
Gate + pilot: 400 runs, ~2–3 min.

**Launch rule, fixed now:** the 200 gate runs are timed both ways under current contention.
Projected wall = 72.6 min × (instrumented / uninstrumented ms per run). If gate PASS and projected
wall ≤ **150 min**, launch TRAIN then TEST without a further stop. If > 150 min, stop and report;
no cut to R or to the band is made unilaterally.

## S5. Analysis unit and label

- **Unit: HOLD epochs** (on = False), epoch index e ≥ 1. Epoch 0 is excluded (every run starts at
  N0 = 3001 with no history; it carries no information about the question).
- **Label y = 1** iff N ≥ 3601.2 at any daily state d_e+1 … d_e+tau, i.e. progression before the
  next decision reads. This matches the main-run `prem` logic exactly (progression on the day of
  the next epoch's read is attributed to the current hold).
- **Truncated intervals** (d_e + tau > 750 with no progression): label undefined, excluded,
  count reported per cell. Main run had 0 censoring in these regions, so this is expected ~0.
- No epoch-level pooling for inference: every CI resamples seed indices (S8).

## S6. Questions and feature sets

N enters every model through the same basis (S7). `cell` = one-hot threshold indicator within tau.

| model | features |
|---|---|
| **A0** (Q-ACT baseline) | N-basis, cell |
| **A1** (Q-ACT test) | A0 + tsf, never_on, g, n_on, nadir |
| **S0** (Q-SPATIAL baseline) | A1 + C + **R** |
| **S1-FULL** | S0 + B, B_R, K, Rg, Rmax |
| **S1-IMAGE** | S0 + B, K, Rg, Rmax (no composition-derived geometry) |
| A0c / A1c (secondary, descriptive) | A0 + C / A1 + C — history beyond carryover |

**Q-ACTIONABLE** = A1 vs A0. **Q-SPATIAL** = S1-FULL vs S0 (primary) and S1-IMAGE vs S0
(secondary). Only Q-SPATIAL can support a claim that spatial state matters.

**Deviation from the brief, flagged for sign-off: R is in the Q-SPATIAL baseline.** The brief's
baseline is N + threshold + history + C. Resistant count is non-spatial state; without it, `B_R`
(and to a lesser degree `B`, since R cells sit at the rim) can register a "spatial" gain that is
composition. The brief asks whether geometry adds "beyond everything non-spatial", which R is.
If you strike this, S0 reverts to A1 + C and a positive S1-FULL must be worded "geometry or
composition".

## S7. Estimation, fixed now

- Per tau, fitted on TRAIN hold epochs only. Logistic regression, L2, C = 1.0 on features
  standardised by train mean/sd, `lbfgs`, max_iter 10,000 (scikit-learn, /opt/anaconda3 python).
  The weak ridge exists only to keep fits finite under quasi-separation at low N.
- **N-basis: cubic B-spline in N, 5 knots at train hold-epoch quantiles, linear extrapolation**,
  identical in every model. **Deviation from "N alone", reason:** progression risk is sharply
  nonlinear in N near the limit; a linear-N baseline is misspecified, and any feature correlated
  with N (g, nadir, Rg, B) would then "gain" by repairing the N curve rather than by carrying
  state. Sensitivity (not a gate): refit all models with 8 knots; report whether any reading changes.
- **Threshold as cell fixed effect** (one-hot, 5 levels per tau) in every model, baseline
  included — fully absorbs the threshold confound rather than assuming it is linear.
- `n_on` enters as log1p; `tsf` in days; others as defined.
- **Hygiene:** all models are fitted on TRAIN, serialised with a SHA-256 of the coefficients, and
  the hashes recorded before the TEST shard directory is opened. The analysis script is developed
  and debugged only on block 2300000 gate data.

## S8. Metrics, uncertainty, bar

Scored on TEST hold epochs, per tau and pooled. Pooled = the concatenation of the four per-tau
models' held-out predictions (no separate pooled model).

- **Primary metric:** relative held-out log-loss reduction Δ = 1 − LL(test model)/LL(baseline).
- Also: AUC and ΔAUC; calibration intercept and slope (logistic recalibration of the held-out
  logit); 10-bin reliability table. Descriptive.
- **CI:** 4,000-resample cluster bootstrap over TEST seed indices (a resample carries every epoch
  of that seed index in every cell of every tau), models held fixed at the train fit. Percentile
  95% CI. Limitation stated now: this CI is conditional on the train fit. **Cross-fit check:**
  swap roles (fit on TEST, score on TRAIN) and report Δ; a reading that does not survive the swap
  is reported with that caveat.

**Reading per comparison (per tau and pooled), three states:**
- **GAIN:** Δ ≥ 0.10 **and** CI lower bound > 0.
- **NULL:** CI upper bound < 0.10 — a material gain is excluded. (If also lower bound > 0, append
  "immaterial gain detected".)
- **INCONCLUSIVE:** neither. **Deviation, reason:** the brief's bar defines GAIN only; "not GAIN"
  is not evidence of sufficiency, and the NULL branch below makes a sufficiency claim, so it
  needs material gain to be excluded, not merely unshown.

**Power floor:** a tau with < 50 events in TRAIN or in TEST is reported UNDERPOWERED with no reading.
**Headline rule (replication):** a pooled GAIN or NULL headline stands only if the same reading
holds at ≥ 2 of the 4 tau individually. Otherwise per-tau readings are reported and no pooled
headline is made. No pattern across tau is claimed.
**Primary comparisons:** Q-ACT pooled and Q-SPATIAL-FULL pooled. All else secondary. No
multiplicity correction; primaries are designated instead.

## S9. Stratified check (matched N), secondary

Per tau, TEST hold epochs binned in N at width 20 cells. A bin is eligible if it holds ≥ 5 events
and ≥ 5 non-events. Within each eligible bin compute AUC of the test model and of the baseline;
statistic = event-weighted mean ΔAUC over eligible bins, same seed-cluster bootstrap. Reported
alongside every primary reading. **Pre-filed consequence:** if a comparison reads GAIN but the
stratified ΔAUC CI includes 0, the gain is flagged as possibly between-bin (N-curve) rather than
state at matched N, and is not described as "at matched N".

Filed now: **whole-sample AUC will be near 1 for every model** because most hold epochs sit far
below the progression limit with no risk. Whole-sample AUC is reported but carries no inferential
weight; within-bin ΔAUC is the matched-N discrimination measure.

## S10. Readings — both branches, filed before compute

**Q-ACTIONABLE (A1 vs A0).**
- GAIN ⇒ clinically available history (epoch sizes and past decisions) improves prediction of
  premature progression over size and threshold alone. **At tau = 7 and 14 this is flagged as
  potentially carryover-driven regardless of outcome** (drug half-life ~7.4 d; C ≈ 0.21 at 7 d
  off, 0.08 at 14 d off, PREREG 2026-08-22 A4.1). A0c/A1c is reported beside it: if history's
  gain vanishes once C is in both models, the tau = 7/14 gain is worded as carryover.
- NULL ⇒ at matched size and threshold, history adds no material predictive information; a
  size-only rule is sufficient for prediction in this model.
- If Q-ACT is GAIN and Q-SPATIAL is NULL: **"temporal information improves on size alone"** —
  worded that way, never as a spatial result.

**Q-SPATIAL (S1 vs S0).**
- GAIN ⇒ spatial state carries predictive information beyond size, history, drug carryover and
  composition. FULL and IMAGE-ONLY reported separately; only IMAGE-ONLY GAIN supports a statement
  about imaging-observable geometry. Prediction only — no mechanism is licensed by it.
- NULL ⇒ size plus non-spatial information is sufficient for prediction here. **This corrects the
  §7 RECAL wording "Total burden is not a sufficient state variable for this control problem"**
  and must be filed as a dated correction to that row and to every WRITEUP_ASSETS sentence that
  repeats it.
- INCONCLUSIVE ⇒ neither the correction nor the claim is licensed; the RECAL sentence is flagged
  "untested at matched N" rather than left standing unqualified.

## S11. Carried constraints and flags

- **Coverage limit, for Step 2:** the models see hold epochs only for N below the top band
  threshold (3300 / 3100 / 2950 / 2980). Any rule built from them is defined only on that range;
  extrapolation above it is out of scope.
- Prediction ≠ control: a GAIN here does not show that a rule improves on the size-only rule
  (Step 2's question).
- q = 0.05 remains an arbitrary model-reliability level with no clinical justification.
- Every number in the results filing is traced to `ode/results/p_s_*.npz` with VERIFIED/UNCERTAIN tags.

## S12. Step 2 note — filed now so it cannot be skipped

Any rule built from these models (hold iff p̂ ≤ p*) must be evaluated by **running it** in the ABM,
not by scoring predictions offline: the rule changes which tumours reach which states, so offline
scores on Step 1 data are not its performance. p* tuning therefore also requires running, which
means Step 2 needs **two fresh blocks** (tune, score), neither of which is 2100000/2200000/2300000.
Step 2 is registered only if Step 1 reads GAIN on Q-ACT or Q-SPATIAL, and only after sign-off.

## S13. Deviations from the brief, collected

1. R (resistant count) in the Q-SPATIAL baseline (S6).
2. Three-state reading with a NULL that requires CI upper < 10% (S8).
3. N entered as a 5-knot spline, not linearly (S7).
4. Threshold as a 5-level cell fixed effect, not a linear covariate (S7).
5. Epoch 0 excluded (S5).
6. Block 2300000 consumed for gate / cost pilot / analysis dry run (S1).
7. Pooled headline requires replication at ≥ 2 tau (S8).
8. Secondary A0c/A1c comparison to interpret the tau = 7/14 carryover flag (S6).

## S14. Constraints in force

Closed files stay closed. No paper prose. Each block read once. Nothing runs until sign-off.

---

## AMENDMENTS — filed 2026-09-23 at sign-off, before any seed of 2100000 / 2200000 / 2300000 was read

**Sign-off record.** All eight S13 items approved as filed. Launch rule S4 confirmed. Blocks
confirmed: GATE 2300000, TRAIN 2100000, TEST 2200000.

### Amendment A — two baselines for the geometry test
- **R-IN** (S0 as filed: A1 + C + R) is the baseline for the SPATIAL ATTRIBUTION claim.
  Comparisons: S1-FULL vs S0 (primary), S1-IMAGE vs S0 (secondary).
- **R-OUT** (no composition: S0-ROUT = A1 + C) is the baseline for any ACTIONABLE image-only
  model, since a clinic cannot observe R. Comparison: S1-IMAGE-ROUT = S0-ROUT + B, K, Rg, Rmax
  vs S0-ROUT. Reported beside R-IN, same metrics, CIs and readings.
- **Step 2's image-only rule and its comparator use R-OUT.**
- **Flag raised by me at filing, not a change to the instruction:** C (mean occupied-site drug
  concentration) is also not clinically observable, and R-OUT as instructed still contains it.
  One additional DESCRIPTIVE comparison is registered so the question can be settled at Step 2
  sign-off without re-analysis: **CLIN-IMAGE = A1 + B, K, Rg, Rmax vs A1** (no R, no C). It
  carries no reading of its own and does not alter any primary or secondary reading.

### Amendment B — bin inclusion for the stratified check (replaces the S9 ≥ 5 / ≥ 5 rule)
A 20-cell N bin (edges at multiples of 20) qualifies only if it holds **≥ 10 events and ≥ 10
non-events on TEST** (fixed on the observed TEST sample; bin qualification is not re-decided per
bootstrap resample). The number of qualifying bins is reported per tau. **If fewer than 3 bins
qualify at a tau, that tau's stratified check is reported UNDERPOWERED, not null.**

### Amendment C — fitting method, fixed now, no change after seeing fits
- Estimator: scikit-learn 1.6.1 `LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
  max_iter=10000, fit_intercept=True)`. Objective = ½‖w‖² + C · Σᵢ logloss_i (a **sum** over
  training epochs), intercept unpenalised. **Fixed ridge, λ = 1/C = 1.0** on standardised
  columns. With ≥ 10⁴ training hold epochs per tau the penalty is ≲ 10⁻⁴ of the data term per unit
  coefficient²; it exists only to keep fits finite under quasi-separation.
- Every column (spline basis, cell indicators, all features) standardised by TRAIN mean and sd
  of that tau; constant columns dropped. Cell indicators: one-hot with the lowest threshold as
  reference.
- **N-basis:** `SplineTransformer(n_knots=5, degree=3, knots="quantile",
  extrapolation="linear", include_bias=False)` **fitted on TRAIN hold-epoch N of that tau only**;
  TEST N passes through the fitted transformer. Sensitivity refit: identical with `n_knots=8`.
- Cross-fit check (S8) re-derives knots and standardisation from whichever block is the fitting
  block in that direction.
- If lbfgs reports non-convergence at 10,000 iterations, the fit is reported as such; the
  estimator is not changed.

### Amendment D — sign agreement for the pooled headline
The ≥ 2 tau that carry a pooled headline must reach the same reading as the pooled comparison
**and** have a point Δ of the **same sign** as the pooled Δ. A pooled NULL whose supporting tau
split in sign does not stand as a headline; per-tau readings are reported instead.

### Execution notes (no design content)
- Runner `ode/p_s_log.py`; analysis `ode/p_s_analyse.py`. Shards of 25 seeds per (role, cell),
  written atomically (**supersedes S3's "50 seeds"**, which conflicted with the §13.5
  checkpoint-granularity rule of 25 runs per unit; no design content changes); the supervisor skips completed shards. Because background processes are
  suspended between tool calls (research_log §13.5 note), the supervisor runs in foreground
  invocations with an internal wall budget; realised wall = sum of invocation walls, JIT included.
- `NUMBA_NUM_THREADS = OMP_NUM_THREADS = MKL_NUM_THREADS = 1`, 7 workers, `pkill -9 -f spawn_main`
  check before each launch (WRITEUP_ASSETS §5.C).
