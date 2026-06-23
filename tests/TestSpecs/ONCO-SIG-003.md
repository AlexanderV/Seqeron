# Test Specification: ONCO-SIG-003

**Test Unit ID:** ONCO-SIG-003
**Area:** Oncology
**Algorithm:** Signature Exposure Estimation — Bootstrap Confidence Intervals
**Status:** ☐ Pending re-validation (Poisson variant added 2026-06-23; verification reset across checklists)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Senkin S. (2021). MSA. *BMC Bioinformatics* 22:540 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/ | 2026-06-14 |
| 2 | Huang X., Wojtowicz D., Przytycka T.M. (2018). *Bioinformatics* 34(2):330–337 | 1 | https://academic.oup.com/bioinformatics/article/34/2/330/4209996 | 2026-06-14 |
| 3 | sigminer `sig_fit_bootstrap` (Wang S. et al.) | 3 | https://raw.githubusercontent.com/ShixiangWang/sigminer/master/R/sig_fit_bootstrap.R | 2026-06-14 |
| 4 | Efron B. (1979). *Annals of Statistics* 7(1):1–26 (percentile method) | 1 | https://doi.org/10.1214/aos/1176344552 | 2026-06-14 |
| 5 | Hyndman R.J., Fan Y. (1996). *The American Statistician* 50(4):361–365 (type-7) | 1 | https://doi.org/10.1080/00031305.1996.10473566 | 2026-06-14 |
| 6 | Senkin S. (2021), MSA — Poisson resampling variant (re-retrieved verbatim quotes) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/ | 2026-06-23 |
| 7 | Knuth D.E. (1997). TAOCP Vol. 2 §3.4.1 (Poisson deviate generation) | 1 | https://en.wikipedia.org/wiki/Poisson_distribution#Random_variate_generation | 2026-06-23 |

### 1.2 Key Evidence Points

1. Parametric bootstrap resamples the catalog as a multinomial draw of N = Σ catalog mutations with pₖ = catalogₖ / N — Senkin 2021; sigminer `sample(K, N, replace=TRUE, prob=catalog/sum(catalog))`.
2. "For each bootstrap sample, NNLS attribution is applied to derive the vector of signature activities." — Senkin 2021.
3. "95% confidence intervals … taking [2.5%, 97.5%] percentiles of the resulting bootstrap activities." — Senkin 2021; Efron percentile method.
4. 1000 replicates is the standard count; ≥ 100 recommended — Senkin 2021 (Fig. 2); sigminer doc.
5. The point estimate is the NNLS exposure of the observed (un-resampled) catalog — Senkin 2021; Huang 2018.
6. Type-7 quantile Q(p) = x₍⌊h⌋₎ + (h−⌊h⌋)(x₍⌊h⌋₊₁₎ − x₍⌊h⌋₎), h = p·(n−1) — Hyndman & Fan 1996.
7. **Poisson variant (Senkin 2021):** "mutations are accumulated following Poisson distributions for each mutation class"; the modified scheme draws "counts from independent binomial distributions, so that the total mutational burden is no longer fixed" and "for any given mutation category … the distribution of bootstrapped mutation counts follows a Poisson distribution." Construction: each channel k resampled as Poisson(catalogₖ), refit by NNLS, same percentile CI. The default scheme remains the fixed-N multinomial (sigminer); the Poisson scheme is opt-in via the `resampling` parameter.

### 1.3 Documented Corner Cases

- **N = 0 (zero-mutation catalog):** every resample all-zero ⇒ all replicate exposures 0 ⇒ interval [0, 0] (Senkin 2021 resampling definition).
- **Single non-zero channel:** multinomial draw is deterministic ⇒ resample equals observed ⇒ each replicate exposure equals the point estimate.
- **Single replicate (R = 1):** percentile of a one-element sample is that element ⇒ lower = upper = mean.
- **Poisson variant — zero-count channel:** Poisson(0) = 0 ⇒ that channel/signature is 0 in every replicate.
- **Poisson variant — single non-zero channel:** Poisson(λ) fluctuates (variance = mean) ⇒ interval has positive width (NOT the deterministic collapse of the multinomial scheme).

### 1.4 Known Failure Modes / Pitfalls

1. Classic (non-parametric) bootstrap is inappropriate here; the parametric multinomial resample of the count catalog is required — Senkin 2021.
2. Proportions cannot be resampled — the integer total N is needed for the multinomial sample size (sigminer uses `sum(catalog)`).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `BootstrapExposures(catalog, signatures, replicates, confidence, seed, resampling)` | OncologyAnalyzer | Canonical | Multinomial (default) or Poisson bootstrap + NNLS refit + percentile CI per signature. |
| `BootstrapResampling` (enum) | OncologyAnalyzer | Canonical | Selects `Multinomial` (default, fixed N) or `Poisson` (Senkin 2021, per-channel Poisson). |
| `Percentile(values, probability)` (internal) | OncologyAnalyzer | Internal | Type-7 quantile; tested directly via reflection and indirectly via deterministic CI cases. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | All bootstrap exposures and interval bounds ≥ 0 | Yes | NNLS x ≥ 0 (Lawson & Hanson 1974); resampled counts ≥ 0 |
| INV-2 | Lower ≤ Upper for every signature | Yes | percentile ordering, ½(1−c) < 1−½(1−c) (Efron 1979) |
| INV-3 | Lower ≤ Mean ≤ Upper for c such that bounds bracket central mass; Mean ≥ 0 | Yes | mean lies within [min,max] ⊇ [2.5%,97.5%]-adjacent; for the deterministic cases equals the bounds |
| INV-4 | Determinism: same (catalog, signatures, replicates, confidence, seed, scheme) ⇒ identical result | Yes | fixed RNG seed |
| INV-5 | One interval per signature, in signature order; PointEstimate = `FitSignatures(observed).Exposures[j]` (scheme-independent) | Yes | contract (Senkin 2021 point estimate) |
| INV-6 | Poisson scheme: observed count 0 ⇒ that channel resamples to 0 every replicate | Yes | Poisson(0) = 0 (Senkin 2021) |
| INV-7 | Default (no `resampling` arg) is byte-for-byte equal to explicit `Multinomial` | Yes | backward compatibility |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Single-channel collapse | catalog [10], signatures [[1.0]], 100 reps, seed 42 | point=mean=lower=upper=10 | Senkin 2021 (multinomial), NNLS fit |
| M2 | Single-channel collapse, 95% bounds equal | same as M1, confidence 0.95 | lower = upper = 10 exactly | deterministic multinomial collapse |
| M3 | Zero-mutation catalog | catalog [0,0,0], signatures [[1,0,0],[0,1,0]] | every signature: point=mean=lower=upper=0 | N=0 corner case |
| M4 | Point estimate equals observed NNLS fit | non-degenerate catalog/signatures | result[j].PointEstimate == FitSignatures(observed).Exposures[j] | Senkin 2021 point estimate; INV-5 |
| M5 | Interval ordering & non-negativity | two-signature non-degenerate catalog, seed 42 | lower ≤ upper, all bounds ≥ 0, lower ≤ mean ≤ upper | INV-1, INV-2, INV-3 |
| M6 | Determinism | two identical calls (seed 42) | element-wise identical intervals | INV-4 |
| M7 | One interval per signature in order | k signatures | result.Count == k | INV-5 |
| M8 | Percentile type-7 [0,1,2,3,4] | via deterministic constant per-rep distribution; verify median bound | exact type-7 value 2.0 at p=0.5 | Hyndman & Fan 1996 |
| M9 | Single replicate R=1 | catalog [10], [[1.0]], replicates 1 | lower=upper=mean=10 | R=1 corner case |
| P1 | Poisson single-channel matches independent Poisson draws | catalog [12], [[1.0]], 300 reps, seed 42, Poisson | mean/bounds == mean/[2.5%,97.5%] type-7 of independent Poisson(12) draws (Knuth reference, same seeded order) | Senkin 2021 Poisson construction; INV (variance=mean) |
| P2 | Poisson zero-count channel stays 0 | catalog [0,20] over [[1,0],[0,1]], Poisson | sig0: point=mean=lower=upper=0; sig1 point=20 | Poisson(0)=0; INV-6 |
| P3 | Poisson single-channel non-degenerate | catalog [25], [[1.0]], 400 reps, Poisson | Upper − Lower > 0 | variance=mean (NOT multinomial collapse) |
| P4 | Poisson determinism | two identical Poisson calls, seed 42 | element-wise identical intervals | INV-4 |
| P5 | Default == Multinomial, ≠ Poisson | same inputs, default vs explicit Multinomial vs Poisson | default == Multinomial byte-for-byte; Poisson differs | INV-7; backward compatibility |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Confidence affects width | wider c ⇒ wider-or-equal interval on a non-degenerate catalog | upper(0.99) ≥ upper(0.5); lower(0.99) ≤ lower(0.5) | percentile monotonicity |
| S2 | Two-channel deterministic split | catalog [0, 7] over signatures [[1,0],[0,1]] | sig0 interval [0,0,0,0]; sig1 [7,7,7,7] | multinomial collapse on the single non-zero channel |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Larger replicate count stability | 1000 reps default still brackets point estimate | lower ≤ point ≤ upper for the dominant signature | informal stability |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior tests exist for `BootstrapExposures` (new method; grep of `src/`, `tests/` confirmed no Oncology bootstrap). Sibling fixture `OncologyAnalyzer_FitSignatures_Tests.cs` (ONCO-SIG-002) establishes conventions.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new |
| M2 | ❌ Missing | new |
| M3 | ❌ Missing | new |
| M4 | ❌ Missing | new |
| M5 | ❌ Missing | new |
| M6 | ❌ Missing | new |
| M7 | ❌ Missing | new |
| M8 | ❌ Missing | new |
| M9 | ❌ Missing | new |
| S1 | ❌ Missing | new |
| S2 | ❌ Missing | new |
| P1 | ❌ Missing | Poisson variant (2026-06-23) |
| P2 | ❌ Missing | Poisson variant (2026-06-23) |
| P3 | ❌ Missing | Poisson variant (2026-06-23) |
| P4 | ❌ Missing | Poisson variant (2026-06-23) |
| P5 | ❌ Missing | Poisson variant (2026-06-23) |
| Failure modes (null/ragged/mismatch/negative/replicates/confidence/undefined-scheme) | ❌ Missing | undefined-scheme added 2026-06-23 |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_BootstrapExposures_Tests.cs` — all ONCO-SIG-003 tests.
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_BootstrapExposures_Tests.cs | Canonical ONCO-SIG-003 fixture | 25 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | ❌ Missing | implemented (deterministic constant distribution) | ✅ Done |
| 9 | M9 | ❌ Missing | implemented | ✅ Done |
| 10 | S1 | ❌ Missing | implemented | ✅ Done |
| 11 | S2 | ❌ Missing | implemented | ✅ Done |
| 12 | Failure modes | ❌ Missing | implemented (7 tests) | ✅ Done |
| 13 | P1 | ❌ Missing | implemented (Poisson vs independent Knuth Poisson reference) | ✅ Done |
| 14 | P2 | ❌ Missing | implemented | ✅ Done |
| 15 | P3 | ❌ Missing | implemented | ✅ Done |
| 16 | P4 | ❌ Missing | implemented | ✅ Done |
| 17 | P5 | ❌ Missing | implemented | ✅ Done |
| 18 | Undefined-scheme failure mode | ❌ Missing | implemented | ✅ Done |

**Total items:** 18
**✅ Done:** 18 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | BootstrapExposures_SingleChannelCatalog_CollapsesToPointEstimate |
| M2 | ✅ Covered | BootstrapExposures_SingleChannel95Percent_BoundsEqualPointEstimate |
| M3 | ✅ Covered | BootstrapExposures_ZeroMutationCatalog_AllIntervalsZero |
| M4 | ✅ Covered | BootstrapExposures_PointEstimate_EqualsObservedNnlsFit |
| M5 | ✅ Covered | BootstrapExposures_NonDegenerate_IntervalsOrderedAndNonNegative |
| M6 | ✅ Covered | BootstrapExposures_SameSeed_IsDeterministic |
| M7 | ✅ Covered | BootstrapExposures_ReturnsOneIntervalPerSignatureInOrder |
| M8 | ✅ Covered | BootstrapExposures_Type7Median_OnConstantSplit_IsExact |
| M9 | ✅ Covered | BootstrapExposures_SingleReplicate_BoundsEqualMean |
| S1 | ✅ Covered | BootstrapExposures_WiderConfidence_GivesWiderInterval |
| S2 | ✅ Covered | BootstrapExposures_TwoChannelDeterministicSplit_ExactExposures |
| P1 | ✅ Covered | BootstrapExposures_Poisson_SingleChannel_MatchesIndependentPoissonDraws |
| P2 | ✅ Covered | BootstrapExposures_Poisson_ZeroCountChannel_StaysZero |
| P3 | ✅ Covered | BootstrapExposures_Poisson_SingleChannel_HasNonZeroSpread |
| P4 | ✅ Covered | BootstrapExposures_Poisson_SameSeed_IsDeterministic |
| P5 | ✅ Covered | BootstrapExposures_DefaultEqualsMultinomial_AndDiffersFromPoisson |
| Failure modes | ✅ Covered | 8 validation tests (null/empty/ragged/mismatch/negative/replicates/confidence/undefined-scheme) |

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Type-7 quantile interpolation (R/NumPy default; Hyndman & Fan 1996) realizes the percentile method | M8, M5, S1, P1 bound values |
| 2 | Fixed default RNG seed 42 for reproducibility (sources do not prescribe a seed) | M5, M6, S1, P1, P4 (randomized cases); deterministic cases M1–M3, M9, S2, P2 are seed-independent |
| 3 | Poisson deviates generated by Knuth multiplication-of-uniforms (TAOCP §3.4.1); the source fixes the Poisson distribution but not the generation algorithm. Non-correctness-affecting for the distribution (any exact Poisson sampler yields the same distribution); the test reproduces this exact sampler only to derive deterministic seeded expected values | P1 |

---

## 7. Open Questions / Decisions

1. **Decision:** the registry titled ONCO-SIG-003 "Signature Exposure Estimation" with `EstimateExposures`, but point exposure estimation (NNLS) was already delivered by ONCO-SIG-002 (`FitSignatures`). The genuine next mutational-signature piece — and the registry's second listed method `BootstrapConfidenceIntervals(spectrum)` — is **bootstrap confidence intervals on exposures**, implemented here as `BootstrapExposures`. External evidence (Senkin 2021; Huang 2018) supersedes the original by-area title; recorded in the registry scope note.
2. **Decision (2026-06-23, limitation fix):** LIMITATIONS.md §2 recorded that only the sigminer fixed-N multinomial resample was implemented, not Senkin's Poisson variant. The Poisson scheme is now added behind a `resampling` enum defaulting to `Multinomial`, so existing behaviour is byte-for-byte unchanged (INV-7) while the Senkin 2021 Poisson construction (per-channel Poisson(observedₖ), non-fixed N) is available opt-in. The limitation row is removed. Because the algorithm changed, prior cross-checklist verification for ONCO-SIG-003 is stale and was reset (☑→☐) pending re-validation; this TestSpec's own status is likewise ☐ pending re-validation.
