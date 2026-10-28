# PRE-REGISTRATION — State sufficiency, Step 2 (policy test)
**Filed 2026-09-23, before any seed of blocks 2400000 / 2500000 / 2600000 was read.
Status: SIGNED OFF 2026-09-23 — P11 items 1–13 approved with the changes in "Sign-off amendments" at the end, filed before any Step 2 seed was read.**

Parent: `PREREG_2026-09-23_state_sufficiency_logging.md` (Step 1), S12. Step 1 read GAIN on both
questions (research_log §7, 2026-09-23). The question here is whether better prediction buys a better
rule. Prediction and control are different claims; nothing in Step 1 licenses this one.

---

## P1. Seed blocks [VERIFIED clean 2026-09-23]

| block | role | seeds | read |
|---|---|---|---|
| **2600000** + s | GATE + COST PILOT (small, not an experiment block; outcomes feed no estimate) | s = 0..9 | once |
| **2400000** + s | TUNE | s = 0..499 | once per setting, by a pre-declared adaptive search (P5) |
| **2500000** + s | SCORE | s = 0..999 | once |

Repo-wide scan for every integer in [2,000,000, 2,999,999]: only 2000000–2300000 (spent) and the literal
`2,999,999` from the Step 1 filing's scan-range text. **Nothing in Step 2 reads 2100000, 2200000 or
2300000.** The block-1200000 main run is used for DESIGN only (the size-only search bracket, P5), as in Step 1.
CRN: seed s is the same IC in every arm, setting and tau of a block. The cluster unit is the seed index.

**Flag, deviation from the brief:** the brief names two blocks. A third, 10-seed GATE block is proposed so
that the bit-identity gate and cost pilot do not read TUNE or SCORE seeds under a second role.

## P2. Models — frozen, no refitting

The Step 1 primary (5-knot) fits in `ode/results/p_s_models.pkl`, key `F5`, verified at load against the
SHA-256 hashes in `ode/results/p_s_fit.json` (`hash5`). Any mismatch aborts.

| rule | Step 1 model | features | status |
|---|---|---|---|
| **R-CLIN** | `CLIN` | N-spline + history (tsf, never_on, g, log_non, nadir) + image geometry (B, K, Rg, Rmax). No R, no C. | **ACTIONABLE — primary** |
| R-HIST | `A1` | N-spline + history | secondary |
| R-FULL | `S1F` | N-spline + history + C + R + B, B_R, K, Rg, Rmax | attribution only, NOT actionable |
| **SIZE-ONLY** | — | hold iff N < thr (the `p_t_pilot.run` convention `on = N >= thr`) | comparator |

## P3. The threshold covariate — well-defined, and the choice does not matter

Every Step 1 model has one-hot threshold-cell indicators. A rule has no threshold, so the indicators are
**fixed at the cell nearest the main-run N_safe: 3120 / 2950 / 2740 / 2740 at tau = 7 / 14 / 28 / 56.**

**This choice cannot change any rule.** The model is logistic-linear in standardised columns, so fixing
the indicators adds one constant c to the logit at every epoch: logit p̂ = f(N, features) + c. The rule
"hold iff p̂ ≤ p*" is "hold iff f ≤ logit(p*) − c", and p* is tuned by running (P5). Any other fixed cell
gives the identical rule family at a shifted p*. The filed cell only fixes the scale on which p* values are
reported. **The covariate does not make the models ill-defined.**

## P4. Rule definition and the two fixes it needs

At decision epoch e (day e·tau), from the same epoch-state reads as Step 1 (`p_s_log.epoch_reads`) and
history recomputed from the rule's own epoch record:

**R-x holds iff e ≥ 1 and N ≤ N_top(tau) and p̂_x ≤ p\*_x(tau); otherwise it treats.**

Two fixes, both needed, both flagged:
1. **Coverage cap, N_top(tau) = 3300 / 3100 / 2950 / 2980** (the top Step 1 band threshold). The models
   were trained only on hold epochs, and those never have N ≥ N_top. Above that the model is pure
   extrapolation (S11 of the Step 1 filing), so the rule treats.
2. **Epoch 0.** The feature g = log(N/N_prev) is undefined at e = 0, and Step 1 excluded e = 0 from
   training. At e = 0 every rule takes **the tuned SIZE-ONLY decision for that tau** (treat iff
   N0 = 3001 ≥ thr_size\*(tau)). This makes epoch 0 identical across arms, so any arm difference comes
   from state information at e ≥ 1, where the models are defined. SIZE-ONLY is therefore tuned first at
   each tau.

Runtime prediction goes through the unmodified `p_s_analyse.design` + `predict_proba` (measured 0.11 ms per
decision). No RNG is touched, and the simulator is only read.

**Estimand, generalised without change of meaning:** prem = 1 iff model-defined progression (N ≥ 3601.2)
occurs inside an interval that began with a HOLD decision. For SIZE-ONLY this is exactly the Step 1 /
main-run `prem`.

## P5. Tuning on TUNE, by running (R = 500 per setting)

**Criterion (every arm, every tau):** the most permissive setting with **point p_prem ≤ 0.05** on TUNE.
The Wilson 95% upper bound is reported at the selected setting.

**Search, deterministic bisection on CRN seeds (8 evaluations per arm × tau):**
- *SIZE-ONLY:* θ = thr. Bracket [N_safe_main − 200, N_top] = [2915, 3300] / [2746, 3100] / [2520, 2950] /
  [2518, 2980]. That gives 2 bracket evaluations + 6 bisections, a resolution of 6–7 cells.
- *Rules:* θ = log10 p\*. Bracket [−5, 0]. At θ = 0 (p\* = 1) the rule holds whenever N ≤ N_top. That is
  SIZE-ONLY at N_top, which exceeds 0.05 at every tau in Step 1, so the upper bracket should fail by
  construction. 2 bracket evaluations + 6 bisections gives a resolution of 0.078 decades (~20 % in p\*).
- Bisection step: if p_prem(mid) ≤ 0.05 then lo ← mid, else hi ← mid. **Selected = final lo.** Every
  failed evaluation lies above it, so the safe set is a prefix of the evaluated points by construction.
- **Bracket failure, fixed now:** lower end already > 0.05 → extend by −200 cells / −2 decades, at most
  twice. Still failing → the arm is **SAFETY-INFEASIBLE** at that tau, reported as a finding, and not scored
  for GAIN. Upper end ≤ 0.05 → select it, flagged **CAP-LIMITED**.
- **Order:** SIZE-ONLY first at each tau (P4 fix 2), then the three rules.
- **Matching check:** the TUNE p_prem at each selected setting is reported for all four arms. If R-CLIN's
  and SIZE-ONLY's differ by more than 0.010 at a tau, that tau is flagged **UNMATCHED**, and its exposure
  difference is reported with that flag.
- Winner's curse is expected: the selected setting's true p_prem is biased above its TUNE estimate. All
  four arms are tuned by the identical one-parameter search, so the bias is shared, not favoured.
- Selected settings are written with a SHA-256 hash before any SCORE seed is read.

## P6. Scoring on SCORE (R = 1000), all four arms at their tuned settings, tau = 7, 14, 28, 56

Per arm × tau:
- p_prem with Wilson 95% CI.
- **Drug exposure:** drug days D (days with u_A on); run length T = min(t_prog, 750); **on-drug fraction
  f = D / T**; total drug days D.
- **Time to model-defined progression:** RMST₇₅₀ = mean of min(t_prog, 750); censoring count reported.
- Paired differences vs SIZE-ONLY on the same seed (Δf, ΔD, ΔRMST, Δp_prem), with a 4,000-resample
  seed-cluster bootstrap and percentile 95% CIs. The relative form is Z = mean Δ / mean(SIZE-ONLY).

**PRIMARY: R-CLIN vs SIZE-ONLY on on-drug fraction f.** *Choice flagged:* f, not total drug days, is
primary. A rule that delays progression accumulates more drug days simply by running longer, so D confounds
exposure with survival. D is reported as secondary. R-HIST and R-FULL vs SIZE-ONLY are secondary.

**Pooled:** the equal-weight mean over the 4 tau of per-tau paired Δ, bootstrapped over seed index (seeds
are common across tau).

## P7. Readings — filed before compute

**Safety, per arm × tau (evaluated first):**
- **SAFETY FAILURE:** the SCORE p_prem Wilson **lower** bound is > 0.05.
- **Deviation from the brief, reason:** the brief's rule (upper Wilson > 0.05 AND point > 0.05) fires
  about half the time on pure sampling noise for any setting tuned to 0.05, SIZE-ONLY included. A setting
  whose true p_prem is 0.050 has a SCORE point above 0.05 with probability ≈ 0.5, and then its upper bound
  is always > 0.05. The brief's criterion is still reported for every arm as **NOMINAL EXCEEDANCE**, with
  SIZE-ONLY as the reference rate for what tuning noise alone produces.
- **SAFETY SHIFT:** paired Δp_prem (arm − SIZE-ONLY) has CI lower bound > 0. Running the rule raised risk
  relative to the comparator. This is the distribution-shift finding.

**Control, R-CLIN vs SIZE-ONLY, per tau.** A GAIN requires no SAFETY FAILURE for R-CLIN **and** risk
non-inferiority: the paired Δp_prem CI upper bound ≤ **+0.02** (margin arbitrary, stated now; *deviation*,
added so that "at matched p_prem" is tested rather than assumed).
- **GAIN (exposure):** Δf CI < 0, and ΔRMST CI not entirely < 0.
- **GAIN (time):** ΔRMST CI > 0, and Δf CI not entirely > 0.
- **TRADE-OFF:** one improves and the other worsens, both with CIs excluding 0. Reported as such, not as GAIN.
- **NO GAIN:** safety conditions hold and neither Δf nor ΔRMST improves with CI excluding 0. *Reportable:
  prediction improved but control did not; better prediction does not automatically buy a better rule.*
  If one worsens with CI excluding 0 and nothing improves, append "(worse)".
- **NOT ASSESSABLE:** R-CLIN is SAFETY FAILURE, SAFETY-INFEASIBLE, or fails the non-inferiority margin at
  that tau. Its exposure difference is reported but carries no reading.

**Pooled headline:** stands only if ≥ 2 of 4 tau reach the same reading with the same sign of Δf (the
Step 1 Amendment D rule). Form: *"At matched ≤ 5 % premature progression, R-CLIN changes on-drug fraction
by Z % [CI] vs the best size-only threshold on held-out seeds."* UNMATCHED and CAP-LIMITED flags travel
with it. No pattern across tau is claimed.

## P8. Diagnostics — descriptive, pre-declared, no readings
- **Self-calibration:** O/E of each rule's own p̂ on its SCORE hold epochs. Does the model's risk
  estimate hold under the policy it drives?
- **N-range shift:** the fraction of each rule's SCORE hold epochs with N outside the fitted spline's
  boundary knots (stored in `p_s_models.pkl`, so no 2100000 data is read).
- The count of cap-forced treatments (N > N_top) per arm.
- *Not included, flagged:* a full per-feature distribution-shift check would need feature ranges from the
  stored TRAIN shards (block 2100000, features only, no outcomes). It is left out to honour "nothing in
  Step 2 reads 2100000"; say so if you want it.

## P9. Gate (block 2600000), VOID on any failure of 1–4
1. SIZE-ONLY arm of the new runner vs `p_t_pilot.run` at thr = 3120 / 2950 / 2740 / 2740 × 10 seeds:
   burden, t_prog, prem, N_end, N_max exactly equal, 40/40.
2. For each rule arm (at p\* = 10⁻², 40 runs each): final RNG state of the instrumented runner equals that of a
   non-logging replica of the same decisions. Non-vacuity check as in Step 1.
3. **Feature identity:** history recomputed from the logged (N, on) sequence by `p_s_analyse.history`
   equals the runtime history exactly, and runtime p̂ equals `p_s_analyse.predict` on the logged features
   exactly.
4. Label identity: Σ(labels over hold epochs) = prem for every run; no label on a TREAT epoch.
5. Cost pilot: SIZE-ONLY and each rule arm timed under current load. TUNE and SCORE re-projected (P10).

## P10. Cost [UNCERTAIN] and launch rules

Basis: contended 650 µs/simulated-day × 1.11 (Step 1 gate overhead) = **722 µs/day**, plus 0.11 ms per rule
decision (≤ 12 ms/run at tau = 7). Run length is unknown for the rules. **Central:** the main-run cell
means near N_safe (410 / 370 / 330 / 305 d, mean 354 d). **Worst:** every run reaches the 750 d horizon.

| stage | runs | central | worst | quiet-machine (Step 1 gate rate, 55 µs/day) |
|---|---:|---:|---:|---:|
| GATE + pilot | ~400 | ~2 min | ~4 min | < 1 min |
| **TUNE** (16 arm×tau × 8 settings × 500; up to 10 settings with extensions) | 64,000 (≤ 80,000) | **~4.6 h** | **~12 h** | ~21 min |
| **SCORE** (4 arms × 4 tau × 1000) | 16,000 | **~68 min** | **~144 min** | ~5 min |

**Launch rules, fixed now:**
- **TUNE:** after gate PASS, the pilot re-projects TUNE at central run length. TUNE launches only if that is
  ≤ **6 h** (proposed; your number). Otherwise stop and report. No unilateral cut to R, settings or arms.
- **SCORE:** after TUNE, the realised TUNE run lengths re-project SCORE. The **150-min rule** applies:
  ≤ 150 → launch without a stop; > 150 → stop and report.
- Execution as in Step 1: shards of 25 seeds, a resuming supervisor, foreground invocations with an
  internal wall budget, and orphan check before each launch. Realised wall is reported against every
  projection.

**Options if the TUNE projection is too heavy, for you to choose, not taken:** (a) drop R-FULL from TUNE and
SCORE, since it is attribution only (−25 %); (b) 6 settings instead of 8 (−25 %, halves the resolution);
(c) R_tune = 300 (−40 %, Wilson half-width at 0.05 grows from ±0.019 to ±0.025).

## P11. Deviations from the brief, collected
1. Third block (GATE 2600000, 10 seeds) (P1).
2. Coverage cap at N_top, and epoch-0 decision tied to the tuned SIZE-ONLY (P4).
3. Threshold covariate fixed at the N_safe cell, with the invariance argument (P3). No fix is needed.
4. SIZE-ONLY holds iff N < thr (code convention), not ≤ (P2).
5. SAFETY FAILURE uses the Wilson lower bound; the brief's rule is kept as NOMINAL EXCEEDANCE (P7).
6. GAIN adds a risk non-inferiority margin (+0.02) and the TRADE-OFF / NOT ASSESSABLE states (P7).
7. Primary exposure metric = on-drug fraction f; total drug days secondary (P6).
8. R_score = 1000 (unstated in the brief); R_tune = 500.
9. Distribution-shift diagnostic limited to N-range + self-calibration, so that no 2100000 data is read (P8).
10. TUNE launch ceiling 6 h (P10).

## P12. Constraints
Closed files stay closed (`quad2/simulator2.py`, `ode/growth_forms.py`, `ode/wellmixed_ctmc2.py`,
`ode/p_t_pilot.py`). Step 1 artifacts are read-only; the models are used as frozen by hash. No paper prose.
Nothing runs until sign-off.

---

## AMENDMENTS — filed 2026-09-23 (second brief), before any seed of 2400000 / 2500000 / 2600000 was read

### Amendment 1 — R-VEL, the minimal actionable rule
- **Model `VEL`:** N-spline (5 knots, TRAIN quantiles) + threshold-cell indicators + dN = N_t − N_{t−1}
  (epoch-observed). Fitted on TRAIN (2100000, already read in Step 1) by the identical Step 1 procedure
  (`p_s_vel.py fit`). **Frozen by SHA-256 in `ode/results/p_s_vel_fit.json` (2026-09-23), before any
  Step 2 block was read.** All four tau converged.
- **Rule R-VEL:** hold iff e ≥ 1, N ≤ N_top(tau), and p̂_VEL ≤ p\*. The same P3 covariate fix, the P4 cap
  and epoch-0 fixes, the P5 tuning and the P6 scoring apply as for the other rules. Reported alongside R-HIST.
- **Added reading (secondary):** recovery ρ = Δf(R-VEL) / Δf(R-HIST), where Δf is the paired mean
  on-drug-fraction difference vs SIZE-ONLY on SCORE, per tau and pooled, with a seed-cluster bootstrap CI.
  - **ρ ≥ 0.80 (point estimate):** the actionable content is size plus growth rate.
  - **Otherwise:** the richer history features matter.
  - **Assessable only** where R-HIST reads GAIN (exposure) and R-VEL passes the P7 safety conditions.
    Otherwise it is reported as NOT ASSESSABLE. The ρ CI is reported descriptively.
- **Prediction-level recovery, on record before the policy run** (frozen fit on the already-read TEST
  data, no simulation; `ode/results/p_s_vel_prediction.json`): R-VEL vs the Q-ACT baseline A0, pooled
  Δ log-loss +17.1% [16.1, 18.0]. R-VEL recovers **35.1% [33.2, 37.0]** of R-HIST's log-loss gain
  (τ = 7: 41.7%, τ = 14: 52.8%, τ = 28: 8.4%, τ = 56: −19.1% [−108.8, 9.5]). Within-bin ΔAUC vs A0:
  +0.227 / +0.280 / +0.093 / −0.046 [−0.180, +0.105].

### Amendment 2 — primary comparison unchanged
**R-CLIN vs SIZE-ONLY on on-drug fraction remains the pre-registered PRIMARY.** R-VEL is secondary and
cannot be promoted to primary after SCORE is seen. *Reconciliation, flagged:* the second brief's
Amendment 1 text says that on ρ ≥ 0.80 "the headline rule is R-VEL"; Amendment 2 forbids promotion.
Filed as: ρ ≥ 0.80 licenses the interpretive statement "the actionable content is size plus growth rate"
as a SECONDARY finding. The pre-registered headline stays the R-CLIN vs SIZE-ONLY comparison.

### Cost update
Five arms instead of four. TUNE becomes 20 arm × tau searches: 80,000 runs central (≤ 100,000 with
bracket extensions), **~5.7 h contended central / ~15 h worst / ~26 min at the Step 1 quiet rate**.
SCORE becomes 5 × 4 × 1000 = 20,000 runs: **~85 min central / ~181 min worst** contended. At worst
case SCORE exceeds 150 min and would stop under the rule. [UNCERTAIN]

---

## SIGN-OFF AMENDMENTS — filed 2026-09-23, before any seed of 2400000 / 2500000 / 2600000 was read

**Approved:** all 13 deviations, with the following changes. These supersede the corresponding text above.

**#7 — risk non-inferiority, restated.** The GAIN readings (exposure or time) require the **upper 95% bound
of the paired difference Δp_prem = p_prem(R-CLIN) − p_prem(SIZE-ONLY)** on SCORE, from the 4,000-resample
seed-cluster bootstrap, to be **≤ +0.02**. The point estimate is always reported. **The 0.02 margin is
arbitrary:** it is a 40% relative increase at the 0.05 target. **If the non-inferiority bound fails, the tau
reads NOT ASSESSABLE, never TRADE-OFF.** Two variants are labelled:
- "(bound fails on width)" when the point Δp_prem is ≤ 0.02;
- "(risk not matched)" when the point Δp_prem is > 0.02.

The same test applies to the secondary rules against SIZE-ONLY.

**#12 — R-VEL is secondary only, as filed** (Amendment 2 reconciliation stands).

**#13 — SCORE launch ceiling raised to 240 min.** All five rules are kept.

**Pre-launch diagnostic (descriptive only, zero simulation).** Two models are fitted on TRAIN by the Step 1
procedure and frozen by hash before any Step 2 block is read:
- **A0R** = A0 + R (N + threshold + resistant count);
- **A1R** = A1 + R (N + threshold + history + resistant count).

Both are scored on the already-read TEST data. Reported: the history gain with R in the baseline (A1R vs
A0R) against the gain without R (A1 vs A0), and within-bin ΔAUC for both. If the history gain largely
disappears once R is conditioned on, observable history is acting as a proxy for the hidden resistant
fraction, and that is said plainly. No reading. Nothing in Step 2 depends on it.
