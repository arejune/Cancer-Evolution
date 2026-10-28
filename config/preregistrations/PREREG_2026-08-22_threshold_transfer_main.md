# PRE-REGISTRATION — Threshold-transfer main run
**Filed 2026-08-22, before block 1200000 was read. Supersedes nothing; amends the
2026-08-21 threshold-transfer pre-registration.**

**Trigger.** The block-1100000 pilot. This filing exists because the pilot showed the
registered grid could not bracket the safe threshold; it is filed under the bracketing
clause of the parent pre-registration, which requires any extension to carry its own
dated pre-registration before it runs.

---

## A1. Amendment 1 — multiplier grid extension

The registered grid m in {0.90, 1.00, 1.10} did not bracket the safe threshold at any
tau. p_prem was 0.0000 at m = 0.90 for all six tau; it was 0.975 / 0.965 / 0.475 at
m = 1.00 for tau = 7 / 14 / 28, and 0.035 / 0.015 / 0.005 at m = 1.10 for tau = 42 /
56 / 84. The crossing of q = 0.05 therefore lies inside (0.90, 1.00) at short tau and
strictly above 1.10 at long tau. One uniform grid cannot cover both.

**Adopted: option (a), per-tau upper extension within m_max(tau).**

Uniform core m in {0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10}, plus per-tau upper
extension, every entry required to satisfy m < m_max(tau) strictly:

| tau | m grid | n cells |
|----:|---|---:|
| 7  | 0.80 0.85 0.90 0.95 1.00 1.05 | 6 |
| 14 | 0.80 0.85 0.90 0.95 1.00 1.05 1.10 | 7 |
| 28 | 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 | 9 |
| 42 | 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 | 12 |
| 56 | 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 | 15 |
| 84 | 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.20 1.30 1.40 1.50 1.60 1.70 1.80 1.90 | 15 |

64 adaptive-therapy cells + 1 continuous-therapy control = 65 cells, R = 500.

**Why (a) and not (b).** Option (b) — keeping the uniform grid and reporting
"above 1.10, not localised" at tau = 42/56/84 — would leave three of six tau
unresolved, and those three are exactly the tau at which the pilot says the threshold
TRANSFERS. The primary estimand is the tau-dependence of the safe threshold. Under
(b), half the tau axis carries a censored bound rather than a point estimate, and no
statement about how the discrepancy varies with tau could be made from a grid that
localises the crossing on only one half of it. (b) would convert the headline from a
measured tau-trend into an interval-versus-point comparison. That cost is larger than
(a)'s cost, which is 23 additional cells of compute inside a budget with headroom.

**Cost of (a), stated in advance.** The m axis is a multiplier of a tau-dependent
quantity N*(tau), so it is already tau-relative and is NOT a common absolute dose axis;
extending it per tau does not make it less comparable than it already was. What (a)
does forfeit is a balanced factorial: cells per tau range from 6 to 15, so any
tau x m interaction test on the full rectangle is unavailable. The primary estimand
(m_safe(tau) by interpolation) does not require balance. No interaction test on the
rectangle will be reported.

## A2. Amendment 2 — structurally degenerate cells excluded

**Criterion, filed as general:** a cell is structurally degenerate if
m >= m_max(tau) = 1.2 * N0 / N*(tau), because the therapy threshold then sits at or
above the progression limit, the drug can never switch on, and the run measures an
untreated tumour rather than a policy. Such a cell returns p_prem = 1 by construction
and is not a measurement.

| tau | N*(tau) (guard-exact, n_max = L^2) | m_max(tau) | N*/N0 |
|----:|---:|---:|---:|
| 7  | 3425.183 | 1.0514 | 1.1413 |
| 14 | 3255.281 | 1.1063 | 1.0847 |
| 28 | 2933.527 | 1.2276 | 0.9775 |
| 42 | 2635.260 | 1.3665 | 0.8781 |
| 56 | 2359.677 | 1.5261 | 0.7863 |
| 84 | 1872.855 | 1.9228 | 0.6241 |

N0 = 3001, progression limit 1.2*N0 = 3601.2, L = 208, n_max = L^2 = 43264.

**tau = 7, m = 1.10 is excluded** (threshold 3767.7 >= 3601.2). It is the only cell of
the registered grid that the criterion removes. **tau = 7, m = 1.05 is retained and
flagged marginal**: threshold 3596.4, a margin of 4.8 cells. The deterministic ODE arm
progresses prematurely in that cell (t = 26 d), which is the criterion behaving as
intended near its boundary, not a stochastic effect.

## A3. Amendment 3 — horizon cut from 1825 d to 750 d

The pilot's gate (v) required 0.10 <= fraction progressing <= 1.00 in every cell. It
read exactly 1.0000 in all 19 cells, so it could only ever have failed from the low
side and carried no information about whether the horizon was adequate. It is replaced
here by the direct quantity: **the observed maximum time to model-defined progression
in the pilot was 547 d**, against a horizon of 1825 d.

**Known consequence, filed in advance.** The pilot did not sample m = 0.80 or 0.85, and
lower m means more treatment and later progression: pilot mean times to progression
were 419.9 d at (tau = 7, m = 0.90) against 138.3 d at (tau = 7, m = 1.00). Some
low-m cells will therefore be censored at 750 d. **Censoring is reported per cell.**
p_prem is computed as (premature events)/R, which is conservative under censoring — a
censored run contributes 0. Because p_prem is already 0.0000 throughout the low-m
region, censoring cannot move the interpolated safe threshold upward, and the direction
of any residual bias is toward under-stating the failure. The conditional form
(premature events)/(runs that progressed) is reported alongside wherever censoring is
non-zero.

## A4. Amendment 4 — the well-mixed CTMC arm

Identical tau x m grid on `ode/wellmixed_ctmc2.py`'s process at matched N0
(2971 S + 30 Ra + 0 Escape), matched growth law, matched n_max = L^2, matched q, and
the same premature-progression estimand.

**Filed reading, not revised after the arm runs:**
- Most of the failure in ODE -> CTMC => the result is about stochasticity, not space.
  Report it in those words; the spatial claim is withdrawn.
- Most of the failure in CTMC -> ABM => the result is spatial. Report the split as a
  fraction of the total discrepancy at each tau.
- Both material => report both magnitudes and do not privilege the spatial one.

### A4.1 Declared before the arm ran: the drug field is NOT matched by the nominal C_eff

The CTMC applies a constant effective concentration C_eff = kappa/(kappa+kappa_cell) =
0.6667. A diagnostic on block 9,900,000 (NOT a registered experiment block) measured
the ABM's mean occupied-site C_A directly and it does not equal that constant:

- it starts at 0.102 on the first day of treatment and reaches 0.68 only by day 20;
- it plateaus at 0.72-0.75, ABOVE the well-mixed asymptote, because the tumour is small
  relative to the diffusion length and draws drug from the saturated surrounding field;
- it DECAYS with a ~7.4 d time constant after treatment stops, so an early holiday is
  not drug-free (0.206 at 7 days off, 0.082 at 14 days off).

At tau = 7 the field never rises above ~0.46 and never falls below ~0.21. Unmatched,
this dose difference would be attributed to "space" by the filed reading, when it is
drug transport rather than contact limitation.

**Instrument:** a scalar first-order pharmacokinetic state
dC/dt = kappa_pk * u - k_el * C, integrated exactly per day, with the kill rate using
the daily mean concentration. Fitted by weighted least squares to the measured field:
kappa_pk = 0.10004 /day, k_el = 0.13420 /day (plateau 0.7454, time constant 7.45 d).
Residuals: ON branch max 0.0331 / rms 0.0136 over 56 ages; OFF branch max 0.0032 /
rms 0.0013 over 90 ages.

**Two CTMC arms are therefore run: CTMC-nominal (constant C_eff) and CTMC-pk.** The
filed two-bucket reading in A4 is evaluated on **CTMC-nominal**, which is the
matched-parameter CTMC the reading refers to. The finer three-way split
ODE -> CTMC-nom (stochasticity) -> CTMC-pk (drug transport) -> ABM (contact limitation
and lattice discreteness) is reported alongside as a declared secondary decomposition.
It is declared here, before block 1200000 was read.

## A5. Correction carried forward

The parent pre-registration's SYSTEMATIC SHIFT reading stated the discrepancy would
grow with tau. The pilot showed the opposite direction: the failure is worst at short
tau and vanishes by tau = 56. The correction stands as filed and is not re-argued here.

## A6. Power limitation carried forward

The occupancy fix (N0 = 120 -> 3001) costs 38% of the tau leverage: the grid spans
1.829x in tau against 2.928x in the historical configuration. A null on tau-dependence
must be reported as "not resolved over a 1.83x span", never as "the shift does not
depend on tau". Filed before the pilot, restated here.

## A7. Analysis plan, fixed before the read

- **Primary estimand:** m_safe(tau), the largest multiplier with p_prem <= q, obtained
  by linear interpolation in p_prem between the two grid points bracketing q. p_prem is
  forced monotone in m by a running maximum before interpolation; monotonicity is a
  property of the policy (a higher threshold can only delay treatment), and this is
  declared so that bootstrap resamples cannot produce non-monotone interpolants.
- **Relative discrepancy from the analytic prediction:** m_safe(tau) - 1.00.
- **Uncertainty:** 4000-resample bootstrap over the 500 SEEDS, not over runs. Seed is
  the cluster unit; CRN means seed s is the same initial condition in every cell, so a
  resample perturbs all cells of a tau coherently. Resamples in which q is not bracketed
  are reported as a percentage rather than dropped silently.
- **Reported per cell:** p_prem with Wilson 95% intervals, censoring count, mean and
  median time to model-defined progression.
- q = 0.05 is an **arbitrary model-reliability level with no clinical justification**.

## A8. Constraints in force

Block 1200000 is read ONCE. Block 1700000 stays untouched and is reserved for a
confirmation read that is a separate decision. `quad2/simulator2.py`,
`ode/growth_forms.py` and `ode/wellmixed_ctmc2.py` are not edited; n_max is threaded
from experiment code in every arm. The 8-hour wall ceiling and the stop rule stand
unchanged: over budget cuts to the reduced grid first, still over reports that the
design cannot resolve the safe threshold and stops. Neither may be relaxed after the
variance is seen.
