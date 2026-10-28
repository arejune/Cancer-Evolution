# PRE-REGISTRATION — RECALIBRATION CONTROL FOR THE THRESHOLD-TRANSFER RESULT

Filed 2026-08-26, **before any run**. Nothing in this document has been executed.
Seed block **1800000**. Supersedes nothing; it is a control for
`RESULTS_2026-08-24_threshold_transfer_main.md`.

## R0. What prompted it

Four outside reviewers independently proposed the same control: recalibrate the
analytic rule to the ABM's own untreated regrowth, then re-test. This document
registers that control, and registers one fact discovered in Checks 1–3 that
changes its predicted outcome and must be on record before the run.

## R1. The fact that changes the prediction

The ODE arm's **own** safe threshold, obtained by bisection on `thr` in
`p_t_arms.ode_cell` (deterministic RK4, stored artifact, no new compute), is not
at m = 1.00:

| tau | N*_analytic | thr_safe(ODE) | m_safe(ODE) | boundary mechanism |
|----:|------------:|--------------:|------------:|--------------------|
|   7 |      3425.2 |        3454.0 |      1.0084 | mid-course, t_prog = 202 d |
|  14 |      3255.3 |        3299.4 |      1.0135 | mid-course, t_prog = 616 d |
|  28 |      2933.5 |        3001.0 |      1.0230 | never fires (thr > N0) |
|  42 |      2635.3 |        2787.7 |      1.0578 | mid-course, t_prog = 373 d |
|  56 |      2359.7 |        2831.6 |      1.2000 | mid-course, t_prog = 429 d |
|  84 |      1872.9 |        2566.2 |      1.3702 | mid-course, t_prog = 477 d |

`nstar_guard` and `ode_cell` integrate the **identical** growth law
(`surface_turnover_guarded`, c = 3.5340, phi_bulk = 0.2712, n_max = 43264). The
+0.84% → +37.02% gap between them therefore **cannot** be a growth-law
calibration error. It is derivation-method error: the derivation prices the worst
case "sitting exactly at thr⁻ at an epoch, then a full tau of free growth", and a
fixed decision calendar never attains that worst case. The shortfall grows with
epoch spacing, which is why the gap grows monotonically in tau.

Splitting the reported discrepancy at this reference point:

| tau | N*→ODE (derivation) | ODE→ABM (transfer) | N*→ABM (as reported) |
|----:|--------------------:|-------------------:|---------------------:|
|   7 |              +0.84% |             -9.83% |               -9.07% |
|  14 |              +1.36% |            -10.71% |               -9.50% |
|  28 |              +2.30% |             -9.37% |               -7.28% |
|  42 |              +5.78% |             -4.00% |               +1.55% |
|  56 |             +20.00% |             -4.00% |              +15.20% |
|  84 |             +37.02% |             -0.15% |              +36.82% |

**The sign reversal is not present in the transfer term.** ODE→ABM is negative at
every tau and shrinks monotonically toward zero. The reversal in the reported
N*→ABM column is the sum of a growing positive derivation term and a shrinking
negative transfer term.

**Registered consequence:** route (i) below is predicted to be close to a no-op,
because the growth law in the derivation is *already* the ABM-fit law from Phase
0.5 (`ode/measure_growth.py`, `ode/fit_growth.py`; c and phi_bulk in
`ode_model2.py:82` come from that fit). Registering this in advance so that a
null result from route (i) cannot be presented as a surprise.

## R2. The two routes

**(i) Growth-law re-fit.** Fit the growth law directly to the ABM's drug-free
N(t) over the burden range the policy actually visits (2300–3650, the span of the
registered grid), rather than over the Phase 0.5 range. Derive a matched
threshold by the same method: integrate regrowth from N* to 1.2·N0, set transit
time = tau, invert. Route (i) is **clean only if** the re-fit residual over
2300–3650 is below 1% in N(t) and the re-fit parameters change the derived N* by
more than 2% at some tau. If neither holds, it is a no-op and route (ii) governs.

**(ii) Empirical transit map (the substantive control).** Measure the ABM's
transit time from N_start to 1.2·N0 **with no drug**, across a grid of starting
burdens N_start ∈ {2200, 2300, …, 3600} (15 points), R = 200 seeds each,
resistant fraction held at the registered 1%, founders placed by the registered
in-disk conversion. Build T_ABM(N_start) as the **median** transit time (median,
not mean, registered now — the distribution is right-skewed by construction).
The recalibrated threshold is N*_ABM(tau) = T_ABM⁻¹(tau), by monotone
interpolation on the measured map; where tau exceeds T_ABM(2200) the cell is
reported as out of range, not extrapolated.

## R3. The re-test

Re-run the multiplier sweep against the recalibrated threshold at all six tau,
grid m ∈ {0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15} of N*_ABM(tau), truncated by
the same structural rule m < 1.2·N0 / N*_ABM(tau). R = 500, t_end = 750 d,
q = 0.05, CRN on seed block 1800000, same estimand (`prem` = progression while
the policy was below threshold at the last epoch). Safe threshold by the same
linear interpolation with p forced monotone, seed-clustered bootstrap, 4000
resamples.

**Also reported, mandatory:** the ODE arm's bisected boundary under the
recalibrated law, so the derivation/transfer split in R1 can be recomputed after
recalibration. Without it the run cannot distinguish which term moved.

## R4. Pre-registered readings

- **DISCREPANCY LARGELY VANISHES** (|N*_ABM,safe / N*_ABM − 1| < 5% at every tau,
  no sign change): growth-law calibration mismatch. The finding becomes "an
  ODE-derived threshold fails in a spatial model because the growth law does not
  transfer, and recalibrating restores it." Real, useful, narrower. **The
  spatial-state claim is withdrawn.**
- **DISCREPANCY AND SIGN REVERSAL PERSIST**: total burden is not a sufficient
  state variable for this control problem, since a threshold calibrated to the
  ABM's own aggregate regrowth still fails. Stated as: the sufficiency failure is
  demonstrated, the spatial mechanism behind it is not.
- **PARTIAL**: report the surviving fraction of the original discrepancy per tau,
  `(N*_ABM,safe/N*_ABM − 1) / (N*_safe/N* − 1)`, to two significant figures. Do
  not round toward either reading. A tau where the original discrepancy is under
  2% is reported as "denominator too small to resolve", not as a fraction.
- **FOURTH READING, registered because R1 predicts it** — the transfer term
  ODE→ABM survives at roughly its present size (−10% at tau = 7 shrinking toward
  0 at tau = 84) while the derivation term N*→ODE collapses. Then the correct
  statement is: the reported sign reversal was an artifact of the reference
  point, the real transfer failure is one-signed and largest at frequent
  monitoring, and the headline is rewritten to that. This is the outcome R1
  predicts and it is registered as a **success** condition for the control, not a
  failure.

## R5. Seed block

`grep` over the repo: 1000000 (`p_r_robust`), 1100000 (pilot), 1200000 (main),
1300000 (`p_x_bracket`), 1400000 (`p_m_mediator`), 1500000 (§23 structural),
1600000 (§23 scoring) are spent. 1700000 is reserved by
`PREREG_2026-08-22_threshold_transfer_main.md` §A8 for a confirmation read of the
main result. **This control uses 1800000**, which appears nowhere in the repo.

1700000 is deliberately not spent here. A recalibration control is a new
question, not a confirmation of the existing one; spending the reserved block on
it would leave the main result with no independent read, which is the one thing
§A8 was written to prevent.

## R6. Constraints in force

`quad2/simulator2.py`, `ode/growth_forms.py`, `ode/wellmixed_ctmc2.py` are not
edited; n_max is threaded from experiment code in every arm. Block 1800000 is
read ONCE. The stop rule stands: over budget cuts the m-grid to
{0.90, 1.00, 1.10} before it cuts R, and if still over, reports that the design
cannot resolve the recalibrated safe threshold and stops. Neither may be relaxed
after the variance is seen.

## R7. Cost

Transit map: 15 N_start × 200 = 3000 runs, horizon ≤ 120 d, ≈ 1/6 the 750-day
cost, ≈ 60 s at the measured 112.8 ms/run for 750-day runs. Re-test sweep:
6 tau × ≤ 7 m × 500 = ≤ 21,000 runs ≈ 2400 s. ODE bisection: seconds. Total
under 45 min wall at 7 workers, well inside the 8-hour ceiling.

---

## R8. AMENDMENT, 2026-08-26, BEFORE BLOCK 1800000 IS READ — MY R1 PREDICTION IS WITHDRAWN

**R1's registered prediction that route (i) is "close to a no-op" is WITHDRAWN. It was wrong,
and it was contradicted by a number already filed before this pre-registration was written.**

`RESULTS_2026-08-24_threshold_transfer_main.md` §R1.2 records, and I have re-verified directly
from `p_t_main_raw.npz`, that the ABM's untreated transit N0 → 1.2·N0 is **60.33 d**
(sd 2.59, median 60.0, n = 500 in each of four independently never-treated cells: tau = 42
m = 1.35, tau = 84 m = 1.70/1.80/1.90 — all identical, as CRN requires). The ODE's is
**24.977 d**. The growth law overpredicts this ABM's drug-free growth by **2.415x at this
initial condition.**

R1 reasoned that because c = 3.5340 and phi_bulk = 0.2712 were fitted to a spatial ABM in
Phase 0.5, the derivation's growth law is "already the ABM-fit law" and a re-fit could not move
much. That inference does not hold: the Phase 0.5 fit was made on `simulator3` in a different
geometry and burden range, and it does not reproduce **this** IC (3001-cell disk, R_IC = 31,
L = 208). Route (i) is therefore **substantive, not a formality**, and its predicted-null status
is removed.

**What R1's split does and does not survive.** The measurement stands — the ODE arm's boundary
is where R1 says it is, and the two-term split is unaffected, because both terms were computed
inside a single growth law. What is withdrawn is only the *forecast* that recalibration cannot
move the derivation term.

**A consequence that must be registered now, because it is predictable and will otherwise look
like a post-hoc rescue.** Under a recalibrated law with a ~60 d free transit, the never-fires
corner moves: at tau = 28 the first post-baseline epoch (day 29) now falls *before* free
progression (~day 60) rather than after it, so the R0.1 floor caveat that currently applies at
tau ≤ 28 is expected to **retreat to tau ≤ 14 or vanish**. If it does, the derivation term
becomes a measurement over more of the grid, which is a strengthening of the design and must not
be reported as though it were a result about tumours.

**Logged, not erased**, per the project's standing convention. The error is the same class as
the one R4 of the results filing describes: an inference that matched the fact it was built on
(the constants came from an ABM fit) and was never checked against a quantity that could
contradict it (the ABM transit time, already measured and already filed).

---

## R9. AMENDMENT, 2026-08-26, AFTER STAGE 1–2, BEFORE THE STAGE-4 SWEEP IS RUN

**Trigger, stated before the amendment is used.** Stage 1 (transit map, 3000 drug-free runs,
zero censoring) and Stage 2 put the recalibrated threshold **4.57% to 46.54% ABOVE** the
registered `N*`. The registered re-test grid `m ∈ {0.85 … 1.15}` was written against the *old*
reference and, at the new one, sits entirely above the ABM's already-filed safe threshold at
tau = 14, 28 and 42. Predicted `m_safe = N*_safe(main run) / N*_recal` = 0.870 / 0.837 / 0.808 /
0.833 / 0.893 / 0.934 at tau = 7/14/28/42/56/84 — three of six below the grid's floor. Run as
registered, the sweep would return "unbracketed" at half the grid and measure nothing.

**Amendment.** The grid is extended **downward** by `{0.70, 0.75, 0.80}`, per tau, truncated by
the unchanged structural rule `m < 1.2·N0 / N*_recal(tau)`. This is the mirror image of the
2026-08-22 pre-registration's own option (a) per-tau *upward* extension and is made on the same
basis: the registered window does not bracket, and the extension is fixed by a rule rather than
by inspection of any outcome. **No outcome from block 1800000's sweep has been observed** — only
the drug-free transit map, which contains no policy and no estimand.

Resulting grid: 7 / 7 / 8 / 9 / 10 / 10 cells at tau = 7/14/28/42/56/84 = **51 cells × R = 500 =
25,500 runs**.

**The stop rule is amended with it**, because the registered fallback `{0.90, 1.00, 1.10}` no
longer brackets: over budget, the grid cuts to the three multipliers nearest the predicted
bracket at each tau, and if still over, reports that the design cannot resolve the recalibrated
safe threshold and stops.

**Registered point prediction, filed now so the sweep can fail it.** Because a safe threshold in
absolute cell counts is a property of the ABM, the policy and tau — not of what the threshold is
divided by — the sweep should return `m_safe × N*_recal` within bootstrap error of the main
run's `N*_safe` = 3114.6 / 2946.1 / 2719.9 / 2676.2 / 2718.3 / 2562.4. If it does not, on an
independent seed block, something is wrong with one of the two runs and the control is void.

**Route (ii) governs the re-test.** Route (i) re-fit to `(c, phi_bulk) = (6.3217, 0.1029)` and
fails the registered clean-if test on residual: relative RMSE in N(t) is **1.294% mean, 1.784%
max**, above the 1% bar, though it clears the second bar easily (moves N* by 3.20–50.02%). The
two routes corroborate each other to within 1.3–2.4% at every tau, which is reported as
agreement between an independent parametric and non-parametric recalibration, not as a pass.
