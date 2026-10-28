# RESULTS FILING — Threshold-transfer main run
**Block 1200000–1200499, read once, 2026-08-24. Design filed in
`PREREG_2026-08-22_threshold_transfer_main.md`; that filing is not amended here.**

Artifacts: `ode/results/p_t_main_raw.npz` (32,500 runs), `ode/results/p_t_main_arms.json`
(ODE + both CTMC arms), `ode/results/p_t_main_cost.json`.

---

## R0. HEADLINE — TWO FINDINGS, NOT ONE
**Filed 2026-08-26. This section SUPERSEDES the single-finding "the analytic threshold is wrong
in both directions / the discrepancy changes sign with tau" framing that was carried in
conversation and put to four outside reviews. That framing is WITHDRAWN. Zero new compute; the
correction comes from bisecting the stored ODE arm, `p_t_main_arms.json` + `p_t_arms.ode_cell`.**

The reported discrepancy was measured against `N*`, the analytic prediction. `N*` is not where
the ODE arm's own safe boundary lies. Bisecting `thr` in the deterministic arm gives:

| tau | N*_analytic | thr_safe(ODE) | m_safe(ODE) | boundary mechanism |
|----:|---:|---:|---:|---|
| 7 | 3425.2 | 3454.0 | 1.0084 | mid-course, t_prog = 202 d |
| 14 | 3255.3 | 3299.4 | 1.0135 | mid-course, t_prog = 616 d |
| 28 | 2933.5 | 3001.0 | 1.0230 | **never fires** (thr > N0) |
| 42 | 2635.3 | 2787.7 | 1.0578 | mid-course, t_prog = 373 d |
| 56 | 2359.7 | 2831.6 | 1.2000 | mid-course, t_prog = 429 d |
| 84 | 1872.9 | 2566.2 | 1.3702 | mid-course, t_prog = 477 d |

Splitting the reported discrepancy at that reference point gives two separate findings:

| tau | DERIVATION N*→ODE | TRANSFER ODE→ABM | total N*→ABM (as previously reported) |
|----:|---:|---:|---:|
| 7 | +0.84% | **−9.83%** | −9.07% |
| 14 | +1.36% | **−10.71%** | −9.50% |
| 28 | +2.30% | **−9.37%** | −7.28% |
| 42 | +5.78% | **−4.00%** | +1.55% |
| 56 | +20.00% | **−4.00%** | +15.20% |
| 84 | +37.02% | **−0.15%** | +36.82% |

**FINDING 1 — TRANSFER.** Spatial structure makes the threshold **anti-conservative at every
monitoring interval**: −9.83% to −0.15%. One-signed at all six tau, shrinking monotonically as
tau grows. The ABM needs a lower threshold than the ODE it was derived from, most at frequent
monitoring, converging to agreement at tau = 84.

**FINDING 2 — DERIVATION.** The analytic method is **conservative by an amount that grows with
tau**: +0.84% to +37.02%. `nstar_guard` and `ode_cell` integrate the *identical* growth law
(`surface_turnover_guarded`, c = 3.5340, phi_bulk = 0.2712, n_max = 43264), so this is not a
calibration error. The derivation prices the worst case "sitting at thr⁻ at an epoch, then a
full tau of free growth"; a fixed decision calendar never attains it, and the shortfall grows
with epoch spacing.

**The apparent reversal in the total is Finding 2 crossing Finding 1.** A growing positive
derivation term crosses a shrinking negative transfer term between tau = 28 and tau = 42.
Nothing about the spatial model reverses sign. Any sentence attributing the reversal to spatial
structure is withdrawn.

**R0.1 CAVEAT ON FINDING 2, load-bearing and not a footnote.** At tau ≤ 28 the ODE boundary is
set by the never-fires / too-late corner, not by a dynamical safety limit: the crossing cells at
tau = 7, 14 and 28 all have t_prog = 26 d against a free-growth transit of 24.977 d, i.e. the
drug never fires before progression. The derivation term is therefore a **floor at tau ≤ 28 and
a measurement at tau ≥ 42**. That matters because tau ≤ 28 is exactly where the transfer term is
largest (−9.83 / −10.71 / −9.37%), so the split is least well resolved precisely where Finding 1
is strongest. Wherever the two-term split is used, this caveat travels with it.

---

## R1. THE THREE HONESTY CONSTRAINTS — filed WITH the result, not as later caveats

These bound what the headline numbers may be said to mean. They are part of the result.

**R1.1 The ODE -> CTMC -> ABM decomposition is decisive only at tau = 7, 14, 28.**
At those tau the total discrepancy in p_prem at m = 1.00 is 0.976 / 0.974 / 0.480 and the
CTMC -> ABM share is 79.5% / 94.9% / 98.8%. At tau = 42, 56, 84 the TOTAL discrepancy is
0.014 / 0.010 / 0.010 — that is 7, 5 and 5 premature events in 500 runs. The "85.7% / 100% /
100% spatial" entries at those tau are ratios of near-zero quantities and carry no
inferential weight. Any statement of the form "the failure is spatial" is licensed by the
three short-tau cells only, and must not be presented as a property of the whole tau grid.

**R1.2 "2.41x slower regrowth" is a statement about THIS ABM against THIS ODE.**
The untreated transit time N0 -> 1.2*N0 is 60.3 d in the ABM (n = 500 in each of four
independently never-treated cells) against 24.98 d in the ODE from which N*(tau) was derived.
The ratio 2.41 describes a disagreement between two models in this repository at this initial
condition on this lattice. It is NOT a claim about tumour growth, not a calibration constant,
and not transportable to another N0, another lattice size, or another growth-form fit. It is
reported because it is the mechanism that organises the TRANSFER term, not as a measurement of
anything biological. **[CORRECTED 2026-08-26 — this sentence originally read "the mechanism
that organises the sign change." There is no sign change in the transfer term; see R0. The
2.41x gap organises the one-signed transfer term only, and the apparent reversal in the total
comes from the derivation term, which this ratio does not touch.]**

**R1.3 The competition ledger describes ADJACENCY, not effect.**
The ledger counts which cell type occupies the blocking neighbour of a division attempt. It
does not measure whether that blocking changes any outcome. No causal or therapeutic reading
is licensed by it. The specific reason this constraint is filed: the pilot reported "95.9% of
resistant blocking slots are other resistant cells" and that number was treated as a general
property of the model. At full R it is tau-specific and collapses monotonically —
interspecific share of resistant-divider blocked slots runs 61.17% (tau = 7) -> 27.19% (14)
-> 8.87% (28) -> 3.83% (42) -> 2.58% (56) -> 0.34% (84) -> 0.06% (CT), tracking the size of
the surviving sensitive population. A single-cell ledger percentage generalises to nothing.

---

## R2. AMENDED DEGENERACY CRITERION — the empirical one supersedes A2

**A2 as filed** declared a cell structurally degenerate if `m >= m_max(tau) = 1.2*N0/N*(tau)`.
That criterion is **necessary but not sufficient** and is superseded by the empirical one.

> **A cell is structurally degenerate if its mean time to model-defined progression equals
> the untreated value (60.3 d in this configuration). The drug never fires, the run measures
> an untreated tumour, and p_prem = 1.000 by construction rather than by measurement.**

**Six cells meet it in the main grid, all admitted by the A2 form:**

| tau | m | threshold | mean t_prog | p_prem |
|---:|---:|---:|---:|---:|
| 7 | 1.05 | 3596.4 | 62.4 | 1.000 |
| 14 | 1.10 | 3580.8 | 69.8 | 1.000 |
| 42 | 1.35 | 3557.6 | 60.3 | 1.000 |
| 84 | 1.70 | 3183.9 | 60.3 | 1.000 |
| 84 | 1.80 | 3371.1 | 60.3 | 1.000 |
| 84 | 1.90 | 3558.4 | 60.3 | 1.000 |

Two distinct mechanisms produce the same signature: a threshold within ~1% of the progression
limit (tau = 7/14/42), and a threshold above N0 combined with tau exceeding the untreated
transit time, so progression arrives before the second monitoring epoch ever occurs
(tau = 84, all three cells).

**Why the pre-run prediction failed, recorded rather than quietly corrected.** Before the run
I predicted which cells would be degenerate using the **ODE** transit time, 24.98 d. That
predicted 17 cells across tau = 28/42/56/84 and was wrong in both directions. The correct
boundary uses the **ABM** transit time, 60.3 d, which selects only tau = 84 for the
second mechanism — matching the observation exactly. The error is the same 2.41x gap that is
the run's headline mechanism (R1.2): I used the model under test to predict the behaviour of
the model doing the testing.

**Effect on the primary estimand: none.** Every m_safe(tau) lies strictly below the lowest
degenerate multiplier at its tau — 0.909 < 1.05, 0.905 < 1.10, 1.016 < 1.35, 1.368 < 1.70.
No degenerate cell enters any bracketing pair. The six cells are excluded from reporting and
retained in the artifact.

---

## R3. REPORTING CORRECTIONS — filed 2026-08-26, zero compute, from stored counts

**R3.1 CTMC → ABM dominance is scoped to the m = 1.00 column.**
At m = 1.00 the ABM exceeds both CTMC arms at every tau. Across the whole 64-cell grid it does
not: the ABM is **below CTMC-nom in 21 of 64 cells and below CTMC-pk in 19 of 64**. The
reversals are large and systematic in the upper half of the grid — at **tau = 84, m = 1.30 the
ABM is 0.020 against CTMC-nom 0.392**, a −0.372 gap; at tau = 56, m = 1.20, 0.148 against 0.564.
Since the estimand searches across m, "spatial structure adds risk" is a statement about the
m = 1.00 column only. In the upper grid the spatial model is *more* tolerant of an over-set
threshold than the well-mixed one.

**R3.2 q is a robustness check, not a caveat.** Safe thresholds re-derived at q = 0.01 / 0.05 /
0.10 by the same interpolation. The crossing between tau = 28 and 42 survives q = 0.05 and
q = 0.10. It does **not** survive q = 0.01 — but that is a resolution limit, not a
contradiction: at tau = 84 the first non-zero cell is 1/500 = 0.002, so q = 0.01 falls inside a
single-count step and the interpolation returns m = 1.0000 by construction. Reportable form:
**survives at q = 0.05 and q = 0.10; q = 0.01 is below the resolution of R = 500 and is
evidence in neither direction.**

**R3.3 Monotonicity holds empirically.** p_prem is monotone non-decreasing in m at all six tau,
with no violations. Every safe set {m : p_prem ≤ q} is contiguous and a prefix of the grid at
all three q. The `np.maximum.accumulate` in `interp_safe` is a no-op on this data — verified by
recomputing without it, identical to four decimals. The supremum is reportable at every tau.

**R3.4 Wording, adopted.**
- "Flat" → **substantially weaker tau-dependence than predicted**: the empirical safe threshold
  declines **17.73%** (3114.6 → 2562.4 over tau = 7 → 84) against the analytic **45.32%**
  (3425.2 → 1872.9).
- "The analytic threshold is wrong" → **threshold transfer failure**. The derivation is correct
  under its own assumptions and non-transportable.
- **The decomposition is order-dependent.** The arms form a chain, not a factorial: there is no
  spatial-deterministic arm and no PK-matched ODE, so the drug-transport increment is observable
  only inside the CTMC and the stochastic increment only first. Increments are path-dependent
  and must not be read as variance shares.
- **Stochasticity contributes materially at tau = 7**: CTMC-nom = **0.200** at m = 1.00 against
  ABM 0.976. A fifth of the event probability arrives before any spatial structure.

**R3.5 The tau = 28–56 plateau, replacing "a dip at 42".**
Seed-clustered paired bootstrap, 4000 resamples, CRN across cells:

| tau | N*_safe | 95% CI |
|---:|---:|---|
| 7 | 3114.6 | [3110.5, 3119.9] |
| 14 | 2946.1 | [2944.8, 2947.6] |
| 28 | 2719.9 | [2702.3, 2748.0] |
| 42 | 2676.2 | [2662.4, 2693.2] |
| 56 | 2718.3 | [2666.4, 2739.8] |
| 84 | 2562.4 | [2499.5, 2630.8] |

Δ(28 − 42) = **+43.7**, CI [+22.3, +72.6], P(>0) = 1.000 — resolved. But Δ(42 − 56) = −42.1,
CI [−66.4, +11.3], P(>0) = 0.057, and Δ(28 − 56) = +1.6, CI [−24.4, +55.3]: tau = 28 and
tau = 56 are indistinguishable. **The correct description is a plateau across tau = 28–56 with
a resolved local minimum at tau = 42 that is smaller than one grid step** (grid spacing is 146.7
N-units at tau = 28 and 131.8 at tau = 42), resting entirely on interpolation between bracketing
cells. The N* = N0 discontinuity sits at **24.9765 d**, just left of tau = 28, so the plateau's
left edge is the first cell after the feasibility switch. No mechanism is claimed for the
minimum without a finer m-grid at tau = 42.

**R3.6 The baseline regime is feasibility, not a growth-rate mismatch.**
The day-0 action is `on = (N0 ≥ thr)`, so the on/off switch is exactly the N* = N0 crossing.
`brentq` on `nstar_guard(tau) − N0` gives **tau_c = 24.9765 d**, identical to the free-growth
transit N0 → 1.2·N0, as it must be. At tau = 7 the policy **starts off drug** with 29.3% of the
band still available — not because a 7-day window is too coarse to arrest growth.

| tau | N*(tau) | N*/N0 | headroom (1.2N0−N*)/(0.2N0) | day-0 action at m = 1.00 |
|---:|---:|---:|---:|---|
| 7 | 3425.2 | 1.1413 | 29.3% | **OFF** |
| 14 | 3255.3 | 1.0847 | 57.6% | **OFF** |
| 28 | 2933.5 | 0.9775 | 111.2% | ON |
| 42 | 2635.3 | 0.8781 | 160.9% | ON |
| 56 | 2359.7 | 0.7863 | 206.9% | ON |
| 84 | 1872.9 | 0.6241 | 288.0% | ON |

**The crossing is not sampled.** The grid brackets it between tau = 14 (OFF) and tau = 28 (ON),
a 14-day-wide bracket around 24.98 d. Every statement about behaviour *at* the crossing is
interpolation across two cells that differ in day-0 action.

---

## R4. WITHDRAWN-AND-CORRECTED INVENTORY — entry, 2026-08-26

**WITHDRAWN: "the analytic threshold is wrong in both directions; the discrepancy changes sign
with tau; spatial structure makes the threshold anti-conservative at short tau and conservative
at long tau."** Carried in conversation, used as the headline, and **put to four outside
reviews in that form.** All four reviews were conducted on the withdrawn framing.

**TRIGGER.** A defensive check on whether the ODE control arm was *degenerate* — whether
p_prem = 0 at m = 1.00 at every tau meant "protective" or "this model cannot generate the
estimand." The arm is not degenerate: it produces the event at every tau. The check found
something else instead — **the control arm's own boundary is not at m = 1.00**, so the reference
point the whole discrepancy was measured against was mis-located, invalidating the headline's
attribution rather than the arm.

**WHAT SURVIVES.** Both underlying measurements. Every p_prem, every safe threshold, every
bootstrap interval in this filing is unchanged. What is withdrawn is the *attribution* of the
sign change to spatial structure and the single-finding presentation.

**WHAT REPLACES IT.** R0's two findings, with R0.1's tau ≤ 28 floor caveat attached wherever
the split is used.

**CLASS.** Same class as the 2026-08-20 MECH entry and the 2026-08-09 Step-2f withdrawals: a
claim that matched the numbers it was built on and was never tested against a quantity that
could contradict it. Here the untested quantity was the control arm's own boundary, which cost
nothing to compute and was never computed until asked for defensively.

**NOTE ON SCOPE OF PROPAGATION.** The withdrawn framing was searched for repo-wide on
2026-08-26. It appears **only** in this file (R1.2's original last sentence, corrected in place
with the original text quoted) and in the 2026-08-26 recalibration pre-registration, which
already carries the corrected version. It never reached `research_log.md`'s Results Log or
`ode/WRITEUP_ASSETS.md` — the threshold-transfer run has **no Results Log row at all**, which is
its own filing gap and is closed by the row added 2026-08-26. `RESULTS…:70`'s "wrong in both
directions" is a **different and still-correct** statement, about the pre-run degeneracy
prediction, and is left standing.
