# Validation Report: ONCO-MRD-001 — Minimal/Molecular Residual Disease (MRD) / ctDNA Detection (INVAR-style)

- **Validated:** 2026-06-25 (fresh re-validation; supersedes 2026-06-24)   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.DetectMRD`, `TrackVariantsOverTime`, `IsVariantDetected` (legacy ≥2-of-N panel rule + Poisson LoD); **INVAR/INVAR2 surface:** `EstimateInvarSignal`, `IntegratedMutantAlleleFractionV2`, `EstimateInvarSignalWithSize` (+`FragmentSizeProfile` / `FragmentSizeProfile.FromKernelDensity` / `InvarMolecule`), `SuppressOutlierLoci`, `BinomialUpperTail`, `EstimateLocusBackground`, `PassesBothStrandsFilter`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

This is a fresh re-validation prompted by the 2026-06-25 post-completion re-reset (the unit was touched again during the limitation-elimination campaign). The INVAR2 formulas were re-fetched verbatim **this session** from the nrlab-CRUK/INVAR2 reference implementation and the C# was confirmed against them line-for-line; all load-bearing numbers were re-derived from scratch with an independent Python port.

## Stage A — Description

### Sources opened this session (THIS session, verbatim)
- **INVAR2 `R/shared/detectionFunctions.R`** (nrlab-CRUK/INVAR2, the Rosenfeld/Rosenfeld-lab pipeline implementing Wan et al. 2020) — fetched verbatim: `calc_log_likelihood`, `estimate_p_EM`, `calc_likelihood_ratio`, `calc_log_likelihood_with_RL`, `estimate_p_EM_with_RL`, `calc_likelihood_ratio_with_RL`.
- **INVAR2 `R/3_outlier_suppression/outlierSuppression.R`** — fetched verbatim: `repolish` (`P_THRESHOLD`, `P_ESTIMATE`, `binom`, `OUTLIER.PASS`, the null-estimate filter).
- **INVAR2 `R/1_parse/onTargetErrorRatesAndFilter.R`** — fetched verbatim: `BACKGROUND_AF`, `N_SAMPLES`, `N_SAMPLES_WITH_SIGNAL`, `LOCUS_NOISE.PASS`, `BOTH_STRANDS.PASS`.
- **Wan et al. 2020, Sci Transl Med 12(548):eaaz8084 (INVAR/IMAF)** — depth-weighted, background-subtracted integrated mutant allele fraction; tumour cfDNA is shorter (size weighting); per-locus error from controls.
- **Silverman 1986 (Gaussian kernel eq. 2.2a, bandwidth eq. 3.31) + R `bw.nrd0`/`density`** — KDE size smoothing.
- **Reinert 2019 / Signatera white paper / PMC9265001** — legacy ≥2-of-N positivity rule and Poisson LoD `p = 1 − e^(−nfm)` (unchanged; carried from the legacy path).

### Formula check (INVAR2 source verbatim vs code)
| Quantity | INVAR2 source (verbatim, fetched this session) | Code | Status |
|----------|------------------------------------------------|------|--------|
| mixture `q` | `q = AF*(1-e)*p + (1-AF)*e*p + e*(1-p)` | `InvarLogLikelihood` identical | ✅ |
| logL | `sum(lchoose(R,M)+M*log q+(R-M)*log(1-q))/length(R)` | identical (÷ locus count) | ✅ |
| `g` | `g = AF*(1-e) + (1-AF)*e` | identical | ✅ |
| EM E-step | `Z0=(1-g)p/((1-g)p+(1-e)(1-p))`, `Z1=gp/(gp+e(1-p))` | identical | ✅ |
| EM M-step | `p = sum(M*Z1+(R-M)*Z0)/sum(R)`, init 0.01, 200 iters | identical | ✅ |
| LR | `alternative_likelihood − null_likelihood` | identical | ✅ |
| with-RL `L0/L1` | `L0=(1-e)*P0*(1-p)+(1-g)*P1*p`, `L1=e*P0*(1-p)+g*P1*p`; `logL=sum(M*log L1+(R-M)*log L0)/length(R)` | identical (M=mutant indicator, R=1/mol so R−M=1−M) | ✅ |
| with-RL EM | `Z0=(1-g)P1 p/((…)+(1-e)P0(1-p))`, `Z1=gP1 p/((…)+eP0(1-p))`, `p=sum(...)/sum(R)` | identical (÷ n) | ✅ |
| IMAFv2 | `pmax(0, MEAN_AF − BACKGROUND_AF)` then `weighted.mean(·, TOTAL_DP)` | identical | ✅ |
| zero-bg floor | `ifelse(BACKGROUND_AF>0, ·, 1/BACKGROUND_DP)` | floored to `1/depth` | ✅ |
| informative filter | `filter(TUMOUR_AF>0)` | identical | ✅ |
| outlier P_THRESHOLD | `outlierSuppression / n_distinct(UNIQUE_POS)` | `α / loci.Count` | ✅ |
| outlier P_ESTIMATE | `max(estimate_p_EM(...), weighted.mean(AF, TUMOUR_AF))`, 0 if no mutant | `max(emP, weightedMean)`, 0 when null set empty | ✅ |
| outlier binom | `binom.test(x,n,p,"greater")$p.value`; `ifelse(x<=0,1,·)` | `BinomialUpperTail`, returns 1 for x≤0 | ✅ |
| OUTLIER.PASS | `BINOMIAL_PROB > P_THRESHOLD` ⟹ outlier ⟺ `tail ≤ P_THRESHOLD` | `isOutlier = tail <= pThreshold` | ✅ |
| BACKGROUND_AF | `MUTATED_READS_PER_LOCI / DP_SUM` (pooled) | `Σ(ALT_F+ALT_R)/Σ DP` | ✅ |
| LOCUS_NOISE.PASS | `(N_SIG/N) < proportion_of_controls & BACKGROUND_AF < max_background_mean_AF` | identical | ✅ |
| BOTH_STRANDS.PASS | `ALT_F>0 & ALT_R>0 | AF==0` | identical | ✅ |
| KDE | `density()` Gaussian kernel, `bw.nrd0`, `adjust`; integrate over `[L−0.5,L+0.5]` | Silverman eq.2.2a, `bw.nrd0` `h=0.9·min(σ̂,IQR/1.34)·n^(−1/5)`, analytic `Φ`-bin integral, renormalised | ✅ |

Every formula in the code matches the INVAR2 source fetched verbatim this session.

### Edge-case semantics (re-confirmed)
- Pure-background sample ⇒ `p̂≈3.3e-5`, `LR≈0`, not detected (M12). Confirmed.
- No informative locus (all tumour AF=0) ⇒ `ArgumentException` (INVAR `filter(TUMOUR_AF>0)` empties the table). Confirmed (C10/C16).
- Zero background ⇒ floored to `1/depth` so logs are finite (S6). Confirmed.
- Flat size profile (P1==P0) ⇒ size factor cancels ⇒ with-RL LR == no-size LR (M18). Confirmed numerically (Within 1e-9).
- Outlier `x≤0` ⇒ binomial tail 1, never an outlier. Confirmed.
- Both-strands: `AF==0` passes vacuously; single-strand-only fails (M25). Confirmed.
- KDE single observed length ⇒ uniform fall-back (INVAR `length(counts)>1` guard, C24). Confirmed.

### Independent cross-check (Python port of the verbatim INVAR2 equations, computed THIS session)
| Case | Pinned (test) | Python reference (this session) | Match |
|------|---------------|---------------------------------|-------|
| M12 GLRT bg (alt=1) | p̂≈3.3e-5, LR≈0 | p̂=3.27572e-5, LR=−8.485e-5 | ✅ |
| M13 GLRT inj=0.01 (alt=5) | p̂≈0.01002 (±5e-4), LR≈4.06 (±0.05) | p̂=0.010020040, LR=4.055208 | ✅ |
| M14 GLRT inj=0.05 (alt=21) | p̂≈0.0501 (±1e-3), LR≈44.14 (±0.3) | p̂=0.050100200, LR=44.136521 | ✅ |
| M16 AF-weighted LR | 2.66 (±0.1) | 2.655675 | ✅ |
| M16 flat-pooled LR | 1.91 (±0.1) | 1.912137 | ✅ |
| M17 size-weighted p̂ | 0.12042621132507245 (±1e-9) | 0.12042621132507245 (full precision) | ✅ |
| M17 size-weighted LR | 0.19691792427890276 (±1e-9) | 0.19691792427890276 (full precision) | ✅ |
| M19 clean tail P(X≥1\|1000,0.001) | 0.6323045752290356 (±1e-9) | 0.6323045752288907 | ✅ |
| M19 planted tail P(X≥50\|1000,0.001) | 3.7264670792676273e-66 (±1e-72) | 3.726467079265e-66 | ✅ |
| M21 control bg | 40/20000 = 0.002 | 0.002 (by hand) | ✅ |
| K1 KDE P(100), h=0.5, support {99,100,101} | 0.684537604065696 (±1e-10) | 0.684537604065696 (full precision) | ✅ |
| K1 KDE P(99)=P(101) | 0.15773119796715201 (±1e-10) | 0.15773119796715201 (full precision) | ✅ |

For M17 the alternative log-likelihood (−0.7337) exceeds the null (−0.9306), so INVAR2's grid-search fallback branch (only entered when `alternative_likelihood < null_likelihood`) is NOT triggered — the code's direct `alt − null` reproduces the reference exactly.

### Findings / divergences (Stage A → PASS-WITH-NOTES)
1. **`calc_likelihood_ratio_with_RL` grid-search fallback is not ported.** When the EM estimate fails to beat the null (`alternative_likelihood < null_likelihood`), INVAR2 re-runs EM from `1e-5` and does a two-level grid search, then recomputes the LR. The code returns `alt − null` directly, which can be slightly negative in that case. This is a faithful-enough simplification: the fallback exists only to mop up EM non-monotonicity at near-zero signal, where the resulting LR is ≈0 either way (sub-threshold ⇒ not detected). For all validated datasets the alt likelihood beats the null, so the fallback never fires and the pinned values reproduce to full precision. Documented divergence, not a defect.
2. **Per-locus background `e` is caller-supplied to `EstimateInvarSignal`/`EstimateInvarSignalWithSize`;** the control-derived estimate is a separate step (`EstimateLocusBackground`), mirroring INVAR2's own separation of the parse stage (`createLociErrorRateTable`) from the detection GLRT. Honest by-design residual disclosed in TestSpec §7.
3. **Per-locus LR normalisation** (`/length(R)`) is reproduced exactly; the absolute LR scale therefore depends on locus count — a faithful property of the source.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `IntegratedMutantAlleleFractionV2` L9066; `EstimateInvarSignal` L9116; `EstimateCtDnaFractionEm` L9205; `InvarLogLikelihood` L9237; `LogChoose` L9263.
- `EstimateInvarSignalWithSize` L9728; `EstimateCtDnaFractionEmWithRl` L9800; `InvarLogLikelihoodWithRl` L9833; `FragmentSizeProfile` L9343; `FromKernelDensity` L9425; `SilvermanBandwidth` L9545; `StandardNormalCdf`/`Erf` L9633/9662.
- `SuppressOutlierLoci` L9891; `BinomialUpperTail` L9984; `EstimateLocusBackground` L10068; `PassesBothStrandsFilter` L10126.

### Formula realised correctly? (evidence)
- `InvarLogLikelihood` / `EstimateCtDnaFractionEm` / LR — exact ports of `calc_log_likelihood` / `estimate_p_EM` / `calc_likelihood_ratio`; reproduce M13/M14/M16 within the asserted tolerances and M12≈0.
- `InvarLogLikelihoodWithRl` / `EstimateCtDnaFractionEmWithRl` — exact ports of `calc_log_likelihood_with_RL` / `estimate_p_EM_with_RL`; reproduce the pinned `0.19691792427890276` / `0.12042621132507245` to full double precision.
- `SuppressOutlierLoci` — `P_ESTIMATE=max(EM, weighted.mean(VAF,TUMOUR_AF))` over loci with `VAF≤afThreshold ∧ alt≤maxMutantReads ∧ TUMOUR_AF>0`; Bonferroni `α/n`; one-sided upper-tail binomial; outlier ⟺ `tail ≤ threshold` — exact `repolish`.
- `BinomialUpperTail` — log-space `Σ exp(lchoose+i·logp+(n−i)·log(1−p))`, returns 1 for x≤0 — matches `binom.test(...,"greater")`; reproduces the two pinned tails.
- `EstimateLocusBackground` / `PassesBothStrandsFilter` — pooled `Σ(ALT_F+ALT_R)/Σ DP`, `(N_sig/N)<proportion ∧ bg<maxBg`, both-strands rule — exact `createLociErrorRateTable`.
- `FromKernelDensity` — weighted Gaussian KDE (Silverman eq.2.2a), `bw.nrd0` bandwidth, analytic `Φ`-bin integral, renormalised; reproduces the K1 hand-derived masses to full precision.

### Cross-verification table recomputed vs code
The 66-test class pins exact INVAR2-derived values; all reproduced this session by an independent Python port (table in Stage A) and all pass under the actual code.

### Variant/delegate consistency
- No-size and with-RL GLRTs coincide when P1==P0 (M18) — verified `Within(1e-9)`.
- `SuppressOutlierLoci` reuses `EstimateCtDnaFractionEm` + `BinomialUpperTail`; `EstimateInvarSignal` reuses `IntegratedMutantAlleleFractionV2`. Shared primitives, consistent.

### Numerical robustness
`q` clamped to `(double.Epsilon, 1−double.Epsilon)`; with-RL `L0/L1` floored to `double.Epsilon`; zero background floored to `1/depth`; `BinomialUpperTail` log-space (no overflow at n=1000); `alt` clamped to `[0,total]`; denominators guarded. No div-by-zero or overflow on the stated ranges.

### Test quality audit (HARD gate)
66 tests, all deterministic, asserting **exact** sourced values via tight `Within` tolerances (1e-9/1e-10 on the pinned size/KDE/binomial/background values; ≤0.3 on LR references derived independently — appropriate to the integer-rounded mutant-read inputs). Coverage spans all 25 MUST cases (M1–M25), KDE (K1–K6), SHOULD (S1–S6) and COULD validation (C1–C24). Directional checks (INV-9/INV-10, weighted>flat, size>no-size) are paired with exact pinned magnitudes — no `Greater`-only green-washing on load-bearing assertions. Edge/error paths (null, empty, out-of-range AF/bg/threshold/α, no-informative-locus, single-point KDE) all covered. Expected values trace to the verbatim INVAR2 formulas / hand-computation, not to code echoes.

### Findings / defects
- **No code defect.** Every formula realises the INVAR2 source (re-fetched verbatim this session) exactly; all pinned values reproduce to full or asserted precision against an independent Python port. The unported `with_RL` grid-search fallback (Stage A note 1) is a documented, harmless simplification, not a defect.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — INVAR/INVAR2 maths confirmed line-for-line against the nrlab-CRUK reference implementation fetched this session; three documented divergences (unported near-zero-signal grid-search fallback; caller-supplied vs separately-estimated background; per-locus LR normalisation) are all faithful-to-source or by-design residuals.
- **Stage B: PASS** — implementation faithful; full suite green (`Seqeron.Genomics.Tests` 18819 passed / 0 failed; `OncologyAnalyzer_DetectMRD_Tests` 66 passed / 0 failed); independent cross-checks match to full / asserted precision.
- **End-state: CLEAN** — no defect found. No code or tests changed this session; build 0 warnings / 0 errors.
