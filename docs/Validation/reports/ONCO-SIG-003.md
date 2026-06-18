# Validation Report: ONCO-SIG-003 — Signature Exposure Estimation: Bootstrap Confidence Intervals

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates, confidence, seed)`; internals `Percentile` (type-7 quantile), `MultinomialResample`, `SampleBinomial`, `Mean`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN

> Note: the orchestrator's session brief guessed this unit was "signature-aetiology mapping
> (SBS1/SBS2/SBS4…)". That is **ONCO-SIG-004** (Mutational Process Classification). ONCO-SIG-003 is
> the bootstrap confidence-interval algorithm on NNLS signature exposures, per its Registry row,
> TestSpec, and Evidence. This report validates the actual unit.

## Stage A — Description

### Sources opened & what they confirm
- **Senkin S. (2021), MSA, BMC Bioinformatics 22:540** (WebFetched PMC8567580) — confirmed verbatim:
  "For each bootstrap sample, NNLS attribution is applied to derive the vector of signature
  activities"; "95% confidence intervals … taking [2.5%, 97.5%] percentiles of the resulting
  bootstrap activities". The MSA resampling is described as accumulating mutations following
  **Poisson** distributions per mutation class (variable total burden).
- **sigminer `sig_fit_bootstrap`** (per Evidence) — resamples via
  `sample(K, total, replace=TRUE, prob=catalog/sum(catalog))`, i.e. a **fixed-N multinomial** draw
  with N = Σ catalog, pₖ = catalogₖ/N, then NNLS-refits each resample. Cites Huang 2018 as the method.
- **Huang, Wojtowicz & Przytycka (2018), Bioinformatics 34(2):330–337** — bootstrap of the
  mutation-count vector to assess decomposition confidence (`bootstrapSigExposures`, fixed-N multinomial).
- **Efron (1979) percentile method** — 95% CI = [2.5%, 97.5%] percentiles of the bootstrap distribution;
  generally [½(1−c), 1−½(1−c)].
- **Hyndman & Fan (1996) type-7 quantile** (WebFetched Wikipedia "Quantile") — R/NumPy default linear
  interpolation; 0-based rank h = p·(n−1), Q(p) = x₍⌊h⌋₎ + (h−⌊h⌋)(x₍⌈h⌉₎ − x₍⌊h⌋₎).

### Formula check
- Resample: Multinomial(N, p), N = Σ catalog, pₖ = catalogₖ/N. Matches sigminer/Huang exactly.
- Percentile interval: Lower = Q((1−c)/2), Upper = Q(1−(1−c)/2) → 2.5/97.5 for c=0.95. Matches Efron/Senkin.
- Quantile estimator: 0-based h = p·(n−1) linear interpolation = type-7. Matches Hyndman & Fan.

### Edge-case semantics check
- **N = 0** → every resample all-zero → all replicate exposures 0 → interval [0,0]. Sourced corner case. ✅
- **Single non-zero channel** → multinomial collapses deterministically → every replicate = point estimate. ✅
- **R = 1** → percentile of a one-element sample is that element → lower = upper = mean. ✅

### Independent cross-check (numbers)
Type-7 worked values, **independently recomputed with `numpy.quantile(method='linear')`** this session:

| Sample | p | type-7 Q(p) |
|--------|------|------|
| [0,1,2,3,4] | 0.025 | 0.1 |
| [0,1,2,3,4] | 0.5 | 2.0 |
| [0,1,2,3,4] | 0.975 | 3.9 |
| [2,4,6,8] | 0.025 | 2.15 |
| [2,4,6,8] | 0.5 | 5.0 |
| [2,4,6,8] | 0.975 | **7.85** |

These match the Evidence artifact's hand-derived table exactly. (A WebFetch page-summarizer reported the
last value as 9.85 — an arithmetic slip using x₍4₎=8 as the interpolation base instead of x₍3₎=6; the
correct type-7 value is 7.85, confirmed by NumPy.)

### Findings / divergences
- **NOTE (PASS-WITH-NOTES):** Senkin's MSA uses an unconditioned **Poisson** resample (variable N), while
  this unit (and its spec) use sigminer's **fixed-N multinomial**. Both are published, peer-reviewed
  approaches; fixed-N multinomial is exactly sigminer `sig_fit_bootstrap` and SignatureEstimation
  `bootstrapSigExposures` (Huang 2018). The divergence is documented in the Registry scope note and is a
  legitimate, source-grounded modelling choice — not a defect.

## Stage B — Implementation

### Code path reviewed
- `BootstrapExposures` — `OncologyAnalyzer.cs:2914–2997`
- `MultinomialResample` (sequential conditional-binomial) — `:3007–3043`
- `SampleBinomial` (sum of Bernoulli) — `:3050–3072`
- `Percentile` (type-7) — `:3094–3116`; `Mean` — `:3075–3084`

### Formula realised correctly?
- Multinomial via conditional binomials: channel k draws Binomial(remaining, pₖ/Σ_{i≥k}pᵢ). Standard,
  exact construction; total conserved; N=0 → all zeros. ✅
- Lower/Upper probabilities `(1−c)/2` and `1−(1−c)/2`. ✅
- `Percentile`: 0-based h = p·(n−1), linear interpolation; n=1 returns the element. This is type-7. ✅
- PointEstimate = `FitSignatures(observed).Exposures` (un-resampled). Matches Senkin/Huang. ✅

### Cross-verification table recomputed vs code
The private `Percentile` was driven directly (reflection) on the two NumPy-confirmed datasets; every
bound matched type-7 to 1e-12 (new test `Percentile_Type7_NonConstantSamples_MatchesHandDerivedValues`).
Degenerate API cases (single-channel collapse, N=0, R=1, two-channel split) reproduce the sourced
point/mean/lower/upper exactly.

### Variant/delegate consistency
Single public method; no `*Fast`/delegate variants. Defaults: replicates 1000, confidence 0.95, seed 42
— all documented and source-attributed.

### Test quality audit (HARD gate)
- **Sourced expectations:** exact `Is.EqualTo(...).Within(1e-10/1e-12)` for all deterministic cases.
- **Defect found & fixed:** the type-7 interpolation — the *only* place exact bound values are set — was
  exercised **only on constant/degenerate distributions** (M1/M2/M3/M8/M9/S2). On a constant sample every
  quantile rule returns the constant, so M8 ("Type-7 median is exact") would still pass against a wrong
  estimator (nearest-rank, 1-based offset, etc.). That is a code-echo defect per the gate. **Fixed** by
  adding `Percentile_Type7_NonConstantSamples_MatchesHandDerivedValues`, asserting exact type-7 values on
  [0,1,2,3,4] and [2,4,6,8] at p∈{0,0.025,0.5,0.975,1}, traced to Hyndman & Fan 1996 and cross-checked
  against `numpy.quantile(method='linear')` — it fails against any non-type-7 rule. Added an honest NOTE
  on M8 clarifying it cannot distinguish the interpolation rule.
- **No green-washing:** no assertion weakened, no tolerance widened, no test skipped, no expected value
  bent to actual output. No code change was needed (the implementation already computes type-7 correctly).
- **Coverage:** all 11 MUST/SHOULD cases (M1–M9, S1–S2), failure modes (null catalog/signatures, empty
  signatures, length mismatch, negative count, replicates<1, confidence∉(0,1)), invariants INV-1..5, and
  now the type-7 interpolation on non-constant samples.
- **Honest green:** full **unfiltered** suite = 6641 passed, 0 failed, 1 skipped (pre-existing
  `MFE_Benchmark_AllScenarios`); `dotnet build` 0 errors; changed test file warning-free.

### Findings / defects
- One test-quality defect (type-7 only tested on constant distributions) — completely fixed this session
  (logged FINDINGS_REGISTER A46). No implementation or description defect.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — formula, percentile method, and edge cases are correct and externally
  sourced; the multinomial-vs-Poisson resample choice is a documented, reference-implementation-backed
  divergence from Senkin's Poisson variant.
- **Stage B: PASS-WITH-NOTES** — code faithfully realises the validated formula (type-7 confirmed against
  NumPy); the one test-quality gap (interpolation tested only on constant samples) was fixed by a new
  source-locked non-constant type-7 test.
- **End-state: ✅ CLEAN** — algorithm fully functional; defect completely fixed; full suite green.
- Follow-up logged: FINDINGS_REGISTER A46.
