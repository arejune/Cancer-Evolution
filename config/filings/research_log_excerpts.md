# Excerpts from the research log

The number ledger (`src/threshold_transfer/report/ledger.py`) checks some values against text filed in the project's research log,
which is part of the research repository and not of this release. This file holds those passages verbatim:
each excerpt is copied character for character from `research_log.md` at research-repository commit `ed27e50`
(file MD5 `2342cc95fe5d77cb03378c5a9566226e`). An excerpt may start or end in the middle of a sentence. The ledger ids that
use each excerpt are listed above it.

## recal.nsafe.7

```
**REGISTERED POINT PREDICTION HOLDS:** N*_safe re-measured on an independent block agrees with the main run **to within 0.27 of one grid step at every tau** (3156.5/2994.1/2723.5/2645.8/2706.9/2581.5 vs 3114.6/2946.1/2719.9/2676.2/2718.3/2562.4)
```

## recal.dmed.7, recal.dmed.14, recal.dmed.42, recal.dq05.42, recal.dmed.56, recal.dq05.84

```
Discrepancy vs N_q05 = **−11.18 / −14.44 / −18.30 / −16.44 / −9.77 / −3.60%** against the median-based −11.87/−14.90/−19.10/−17.60/−11.09/−5.94
```

## recal.nq05.7, recal.nq05.14, recal.nq05.28, recal.nq05.42, recal.nq05.56, recal.nq05.84, share.84

```
Inverting the **5th percentile** instead moves N*_recal by just **−0.53% to −2.43%** (3553.7 / 3499.4 / 3333.3 / 3166.5 / 3000.0 / 2677.8)
```

## recal.dq05.7, recal.dq05.14, recal.dq05.28, recal.dq05.56

```
q05 discrepancy −11.18 / −14.44 / −18.30 / −9.77 → −9.97 / −14.15 / −17.23 / −9.76 %.** Every qualitative claim survives (anti-conservative at all four tau; shares move ≤ 1.3 points)
```

## recal.nsafe.14

```
ISTERED POINT PREDICTION HOLDS:** N*_safe re-measured on an independent block agrees with the main run **to within 0.27 of one grid step at every tau** (3156.5/2994.1/2723.5/2645.8/2706.9/2581.5 vs 3114.6/2946.1/2719.9/2676.2/2718.3/2562.4)
```

## recal.nsafe.28

```
POINT PREDICTION HOLDS:** N*_safe re-measured on an independent block agrees with the main run **to within 0.27 of one grid step at every tau** (3156.5/2994.1/2723.5/2645.8/2706.9/2581.5 vs 3114.6/2946.1/2719.9/2676.2/2718.3/2562.4)
```

## recal.dmed.28, recal.dmed.84

```
**READINGS.** **(1) discrepancy largely vanishes: FAILS** — post-recalibration discrepancy is **−5.94% to −19.10%, every tau above the 5% bar**; the spatial-state claim is NOT withdrawn
```

## recal.nsafe.42

```
PREDICTION HOLDS:** N*_safe re-measured on an independent block agrees with the main run **to within 0.27 of one grid step at every tau** (3156.5/2994.1/2723.5/2645.8/2706.9/2581.5 vs 3114.6/2946.1/2719.9/2676.2/2718.3/2562.4)
```

## recal.nsafe.56

```
ION HOLDS:** N*_safe re-measured on an independent block agrees with the main run **to within 0.27 of one grid step at every tau** (3156.5/2994.1/2723.5/2645.8/2706.9/2581.5 vs 3114.6/2946.1/2719.9/2676.2/2718.3/2562.4)
```

## recal.nsafe.84

```
DS:** N*_safe re-measured on an independent block agrees with the main run **to within 0.27 of one grid step at every tau** (3156.5/2994.1/2723.5/2645.8/2706.9/2581.5 vs 3114.6/2946.1/2719.9/2676.2/2718.3/2562.4)
```

## ctmc.7, ctmc.14, ctmc.28, ctmc.42, ctmc.56, ctmc.84

```
**Refined CTMC boundary 3513.2 / 3453.2 / 3321.6 / 3233.9 / 3078.7 / 2998.6** (shifts +101.5 / +98.7 / +49.4 / +13.9 / +14.4 / +110.1)
```

## share.7, share.14, share.28, share.56

```
Spatial share 89.4 / 94.2 / 85.0 / 45.9 → 88.1 / 94.1 / 84.2 / 45.9 %
```

## share.42

```
**REFINED SPATIAL SHARE 89.4 / 94.2 / 85.0 / 78.4 / 45.9 / 99.4%.** **THE SPLITS MOVED MATERIALLY**, upward at five of six and flat at tau=56; stated plainly rather than as vindication — the coarse numbers wer
```

## band.7

```
Variance-share bands: [10.9,11.4] / [6.1,6.9] / [12.9,15.8] / [20.0,22.4] / **[47.1,58.3]** / [2.0,2.0] — **only tau=56 straddles 50%.** **CHECK B — tau=56 IS NOT FLOOR-CONTAMINATED; the flag d
```

## band.14

```
Variance-share bands: [10.9,11.4] / [6.1,6.9] / [12.9,15.8] / [20.0,22.4] / **[47.1,58.3]** / [2.0,2.0] — **only tau=56 straddles 50%.** **CHECK B — tau=56 IS NOT FLOOR-CONTAMINATED; the flag did not miss
```

## band.28

```
Variance-share bands: [10.9,11.4] / [6.1,6.9] / [12.9,15.8] / [20.0,22.4] / **[47.1,58.3]** / [2.0,2.0] — **only tau=56 straddles 50%.** **CHECK B — tau=56 IS NOT FLOOR-CONTAMINATED; the flag did not miss it.** ODE fail
```

## band.56

```
Variance-share bands: [10.9,11.4] / [6.1,6.9] / [12.9,15.8] / [20.0,22.4] / **[47.1,58.3]** / [2.0,2.0] — **only tau=56 straddles 50%.** **CHECK B — tau=56 IS NOT FLOOR-CONTAMINATED; the flag did not miss it.** ODE failure mode at the boundary cell,
```

## recal.cv

```
RESIDUAL LARGELY SURVIVES.** Transit spread is small: median **CV 4.49%**, p5 sits only **4.7 d** below the median, 5–95 band 5–13 d
```

## cred.res

```
fied clean; step 0.0025 = 6.9–9.0 cells; 126 cells, 994.7 s).** One coarse 0.05 m-step was **137–179 cells against ODE→CTMC terms of −143.8 / −175.9 / −112.5 — resolution 91–125% of the term.** Scoped to tau=7/42/84, then widened to ALL SIX on inspection: the q=0.05 crossing was bracketed by a single 0.05 step everywhere (tau=14 0.000→0.718, tau=28
```

## cred.basis

```
New basis — **179.9 ms/run** for a 750-d ABM AT cell at 7 workers session-concurrent, **650 µs/simulated-day** (vs 422 µs idle, contention factor **1.54x** at identical workload), **43.9 ms/run** for a drug-free
```
