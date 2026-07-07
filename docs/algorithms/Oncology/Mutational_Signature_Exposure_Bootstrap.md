# Mutational Signature Exposure Bootstrap Confidence Intervals

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-SIG-003 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-23 |

## 1. Overview

This algorithm quantifies the uncertainty of mutational-signature exposures (activities) by the **parametric bootstrap**. Given an observed integer mutational catalog and a set of caller-supplied reference signatures, it repeatedly resamples the catalog, refits each resample to the signatures by non-negative least squares (NNLS), and reports a two-sided percentile confidence interval for each signature's exposure [1][3]. Two resampling schemes are supported, selected by the `resampling` parameter: the **fixed-N multinomial** scheme (default — a draw of N = Σ catalog mutations from Multinomial(N, p), the sigminer scheme [3]) and the **Poisson** scheme (each channel count drawn independently as Poisson(observedₖ), the Senkin 2021 MSA Poisson variant where the total burden is not fixed [1]). It is a probabilistic (Monte-Carlo) method: outputs depend on a random resampling that is made reproducible with a fixed seed. It complements the point exposure estimate of ONCO-SIG-002 (`FitSignatures`) by answering "how stable is each exposure?" — the standard confidence measure for signature attribution [1][2].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A single-base-substitution (SBS) mutational catalog d is a vector of per-channel mutation counts (e.g. 96 trinucleotide channels). Signature *refitting* decomposes d into non-negative contributions x (exposures) of known reference signatures S: minₓ‖S·x − d‖², x ≥ 0 (ONCO-SIG-002). A point estimate alone does not convey how robust each exposure is to sampling noise in a finite mutation count; bootstrap resampling provides that confidence measure [1][2].

### 2.2 Core Model

Let d ∈ ℤ≥0ⁿ be the observed catalog with total N = Σₖ dₖ and per-channel probabilities pₖ = dₖ / N. The parametric bootstrap models the mutations as accumulating by Poisson/multinomial processes per channel [1]. For each replicate b = 1…R, draw a resampled catalog d⁽ᵇ⁾ by the selected scheme, then refit x⁽ᵇ⁾ = argminₓ‖S·x − d⁽ᵇ⁾‖², x ≥ 0 (NNLS) [1]:

- **Multinomial** (default): d⁽ᵇ⁾ ~ Multinomial(N, p) — distribute N mutations over the channels with probabilities pₖ [1][3]. (sigminer realises this as `sample(K, N, replace=TRUE, prob=d/Σd)` then tabulating [3].) The total burden is held fixed at N.
- **Poisson** (Senkin 2021 variant [1]): each channel is drawn independently d⁽ᵇ⁾ₖ ~ Poisson(dₖ), so "for any given mutation category … the distribution of bootstrapped mutation counts follows a Poisson distribution" and "the total mutational burden is no longer fixed" [1]. By the Poisson property the variance of each channel equals its mean (the observed count). Senkin notes the two schemes are linked because "the conditional distribution of a vector of independent Poisson variables is equivalent to multinomial distribution" [1].

For signature j, the bootstrap distribution {x⁽ᵇ⁾ⱼ}ᵇ gives the two-sided confidence interval at level c by the **percentile method** [4]:

```
lower_j = Q_j( (1−c)/2 )
upper_j = Q_j( 1 − (1−c)/2 )
```

where Q_j is the empirical quantile of the replicate exposures. For c = 0.95 these are the 2.5th and 97.5th percentiles [1][4]. The point estimate is xⱼ from the NNLS fit of the un-resampled observed catalog d [1][2].

The empirical quantile Q(p) over a sorted sample x₍₀₎ ≤ … ≤ x₍ₙ₋₁₎ uses the type-7 convention with 0-based rank h = p·(n−1) [5]:

```
Q(p) = x₍⌊h⌋₎ + (h − ⌊h⌋)·(x₍⌊h⌋₊₁₎ − x₍⌊h⌋₎)
```

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Mutations accumulate as independent Poisson/multinomial processes per channel [1] | Resampling variance misestimates true uncertainty; intervals miscalibrated |
| ASM-02 | The supplied reference signatures span the true generating process | NNLS misattributes; bootstrap intervals reflect model misspecification, not just noise |
| ASM-03 | The catalog is given as integer counts (not proportions) | The multinomial sample size N is undefined; resampling cannot proceed |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | All exposures and bounds ≥ 0 | NNLS enforces x ≥ 0 [6]; multinomial/Poisson counts ≥ 0 |
| INV-02 | lower_j ≤ upper_j | (1−c)/2 < 1−(1−c)/2 and Q is monotone in p [4] |
| INV-03 | Determinism for fixed (d, S, R, c, seed, scheme) | RNG seeded from a fixed value |
| INV-04 | One interval per signature, in signature order; PointEstimate_j = NNLS exposure of d (scheme-independent) | by construction [1] |
| INV-05 | N = 0 ⇒ every interval is [0,0] with point 0 (both schemes) | empty multinomial / all Poisson means 0; NNLS of zero vector is zero |
| INV-06 | Poisson scheme: a channel with observed count 0 resamples to 0 every replicate | Poisson(0) = 0 [1] |

### 2.5 Comparison with Related Methods

| Aspect | Multinomial bootstrap (this) | Bayesian credible intervals (sigfit/signeR) |
|--------|------------------------------|---------------------------------------------|
| Uncertainty source | Resampling the observed count vector | Posterior over exposures given a generative model |
| Refit engine | NNLS per replicate | MCMC / variational inference |
| Output | Percentile CI per signature | Posterior credible interval per signature |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| catalog | `IReadOnlyList<int>` | required | Observed per-channel mutation counts | non-null, length = channel count, each ≥ 0 |
| signatures | `IReadOnlyList<IReadOnlyList<double>>` | required | Reference signatures (one channel-vector each) | non-null, ≥ 1, equal-length, non-empty |
| replicates | `int` | 1000 | Number of bootstrap resamples | ≥ 1 |
| confidence | `double` | 0.95 | Two-sided confidence level | in (0, 1) |
| seed | `int` | 42 | RNG seed for resampling | any int |
| resampling | `BootstrapResampling` | `Multinomial` | Resampling scheme: fixed-N `Multinomial` or per-channel `Poisson` (Senkin 2021) | a defined enum value |

### 3.2 Output / Return Value

`IReadOnlyList<ExposureConfidenceInterval>`, one per signature in signature order:

| Field | Type | Description |
|-------|------|-------------|
| PointEstimate | `double` | NNLS exposure of the observed (un-resampled) catalog |
| Mean | `double` | Mean of the signature's exposure over replicates |
| Lower | `double` | (1−c)/2 percentile of the replicate exposures |
| Upper | `double` | 1−(1−c)/2 percentile of the replicate exposures |
| Confidence | `double` | The confidence level c used |

### 3.3 Preconditions and Validation

Counts are 0-based channel-indexed integers. Validation: null `catalog`/`signatures` (or a null signature vector) → `ArgumentNullException`; empty/ragged signatures, catalog length ≠ channel count, or a negative count → `ArgumentException`; `replicates < 1`, `confidence ∉ (0,1)`, or an undefined `resampling` enum value → `ArgumentOutOfRangeException`. Channel count is taken from the first signature; all signatures must match it.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs (including the `resampling` enum value); compute the observed count vector and N = Σ catalog.
2. Compute the point estimate by NNLS on the observed catalog (`FitSignatures`) — scheme-independent.
3. For each of R replicates: draw d⁽ᵇ⁾ by the selected scheme (Multinomial(N, p) or per-channel Poisson(dₖ)), refit by NNLS, record per-signature exposures.
4. For each signature: compute the mean and the [(1−c)/2, 1−(1−c)/2] type-7 percentiles of its replicate exposures.
5. Return one interval per signature.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Replicate count:** default 1000 (Senkin 2021, Fig. 2); ≥ 100 recommended (sigminer) [1][3].
- **Percentile bounds:** 2.5%/97.5% for the default 95% level [1][4].
- **Quantile estimator:** type-7 linear interpolation, h = p·(n−1) [5].
- **Multinomial draw:** sequential conditional-binomial construction — channel k receives Binomial(remaining, pₖ / Σ_{i≥k} pᵢ); Binomial sampled as a sum of Bernoulli trials.
- **Poisson draw:** each channel k drawn independently as Poisson(dₖ) [1], realised by Knuth's multiplication-of-uniforms generator [7]; λ = 0 returns 0.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| BootstrapExposures | O(R·(N + NNLS(n,k))) | O(R·k + n) | R replicates; each resamples (Multinomial: ~N Bernoulli trials; Poisson: ~Σλ = N uniform draws total) and runs one NNLS fit (n channels, k signatures) |
| Percentile (per signature) | O(R log R) | O(R) | one sort of the replicate distribution |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates, confidence, seed, resampling)`: the public entry point; returns one `ExposureConfidenceInterval` per signature. The `resampling` parameter defaults to `Multinomial`, preserving the historical behaviour byte-for-byte.
- `OncologyAnalyzer.BootstrapResampling`: enum selecting the `Multinomial` (default) or `Poisson` scheme.
- `OncologyAnalyzer.ExposureConfidenceInterval`: the per-signature result record.
- Internal: `MultinomialResample`, `SampleBinomial`, `PoissonResample`, `SamplePoisson`, `Percentile`, `Mean` (private helpers); reuses `FitSignatures` / NNLS from ONCO-SIG-002.

### 5.2 Current Behavior

Reference signatures are **not** hardcoded — they are caller-supplied (Framework status), consistent with ONCO-SIG-001/002. The point estimate is exactly `FitSignatures(observed).Exposures` and is independent of the resampling scheme. The multinomial resample is implemented by the conditional-binomial construction rather than R's `sample`+`table`; both produce a Multinomial(N, p) draw, so the distribution is identical (only the RNG stream differs). The Poisson resample draws each channel independently as Poisson(dₖ) via Knuth's generator [7]; the default scheme remains multinomial so existing callers are unaffected. This is **not** a search/matching unit, so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Multinomial resampling of the count catalog with N = Σ catalog, pₖ = dₖ/N [1][3].
- Poisson resampling: each channel drawn independently as Poisson(dₖ); total burden not fixed — the Senkin 2021 MSA Poisson variant [1].
- NNLS refit per replicate [1].
- Percentile confidence interval [(1−c)/2, 1−(1−c)/2]; 2.5%/97.5% at 95% [1][4].
- Point estimate = NNLS fit of the observed catalog [1][2].
- Default 1000 replicates [1].

**Intentionally simplified:**

- Binomial and Poisson draws use exact definitional samplers (Bernoulli-sum for Binomial; Knuth multiplication-of-uniforms for Poisson [7]) rather than O(1) inversion/BTPE/PTRS samplers; **consequence:** per-replicate resampling cost grows with N / Σλ — fine for catalog magnitudes, slower for very large N.

**Not implemented:**

- The MSA Gaussian noise model (default in MSA, σ = 10% of mutational burden) [1]; **users should rely on:** the two count-based schemes here (multinomial, Poisson), which are the discrete-count resamplers; the Gaussian model is a continuous-noise alternative outside the count-resampling scope of this unit.
- Bayesian credible intervals and p-value tests of signature presence; **users should rely on:** external tools (sigfit, signeR, sigminer `report_bootstrap_p_value`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Type-7 quantile convention | Assumption | Sets exact bound values | accepted | R/NumPy default; ASM not in §2.3 (it is a quantile-estimator choice), see Evidence Assumption 1 |
| 2 | Fixed default seed 42 | Assumption | Reproducibility | accepted | sources prescribe no seed; deterministic cases are seed-independent |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| N = 0 (all-zero catalog) | All intervals [0,0], point/mean 0 | empty multinomial / all Poisson means 0 ⇒ zero resamples ⇒ NNLS(0)=0 (INV-05) |
| Single non-zero channel (Multinomial) | Every replicate exposure = point estimate; lower=upper=mean=point | multinomial over one outcome is deterministic |
| Single non-zero channel (Poisson) | Interval has positive width centred near the point estimate | Poisson(λ) fluctuates (variance = mean) |
| Zero-count channel (Poisson) | That channel/signature is 0 in every replicate | Poisson(0) = 0 (INV-06) |
| replicates = 1 | lower = upper = mean = the single replicate exposure | percentile of a 1-element sample |
| null / ragged / mismatched / negative input / undefined scheme | Argument*Exception per §3.3 | input contract |

### 6.2 Limitations

Confidence intervals are conditional on the supplied signature set (ASM-02): an incomplete or wrong reference panel yields tight but biased intervals. The parametric bootstrap assumes per-channel Poisson/multinomial independence (ASM-01); over-dispersed or correlated processes are not modelled. The method requires integer counts, not proportions (ASM-03). Very large N makes the exact Bernoulli-sum sampler slow.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
int[] catalog = { 30, 10, 0, 5 };
var signatures = new IReadOnlyList<double>[]
{
    new double[] { 1, 0, 0, 0 },   // signature A
    new double[] { 0, 1, 0, 1 },   // signature B
};

// Default fixed-N multinomial scheme:
var intervals = OncologyAnalyzer.BootstrapExposures(
    catalog, signatures, replicates: 1000, confidence: 0.95, seed: 42);

// Senkin 2021 Poisson variant (per-channel Poisson(observed); total burden not fixed):
var poissonIntervals = OncologyAnalyzer.BootstrapExposures(
    catalog, signatures, replicates: 1000, confidence: 0.95, seed: 42,
    resampling: OncologyAnalyzer.BootstrapResampling.Poisson);

// intervals[0].PointEstimate == FitSignatures(catalog).Exposures[0]
// intervals[0].Lower <= intervals[0].Upper, both >= 0
```

**Numerical / biological walk-through:**

For catalog `[10]` and signature `[[1.0]]`: N = 10, p = [1]. Every multinomial draw assigns all 10 mutations to the single channel ⇒ d⁽ᵇ⁾ = [10] for all b. NNLS of [10] onto [[1.0]] gives exposure 10. The bootstrap distribution is the constant {10}, so lower = upper = mean = point = 10 at any confidence level.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_BootstrapExposures_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_BootstrapExposures_Tests.cs) — covers INV-01…INV-05
- Evidence: [ONCO-SIG-003-Evidence.md](../../../docs/Evidence/ONCO-SIG-003-Evidence.md)
- Related algorithms: [Mutational_Signature_Fitting](./Mutational_Signature_Fitting.md) (ONCO-SIG-002); [SBS96_Trinucleotide_Context_Catalog](./SBS96_Trinucleotide_Context_Catalog.md) (ONCO-SIG-001)

## 8. References

1. Senkin S. 2021. MSA: reproducible mutational signature attribution with confidence based on simulations. *BMC Bioinformatics* 22:540. https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/
2. Huang X., Wojtowicz D., Przytycka T.M. 2018. Detecting presence of mutational signatures in cancer with confidence. *Bioinformatics* 34(2):330–337. https://academic.oup.com/bioinformatics/article/34/2/330/4209996
3. Wang S. et al. sigminer `sig_fit_bootstrap`. https://shixiangwang.github.io/sigminer-doc/sigfit.html ; source https://raw.githubusercontent.com/ShixiangWang/sigminer/master/R/sig_fit_bootstrap.R
4. Efron B. 1979. Bootstrap Methods: Another Look at the Jackknife. *Annals of Statistics* 7(1):1–26. https://doi.org/10.1214/aos/1176344552
5. Hyndman R.J., Fan Y. 1996. Sample Quantiles in Statistical Packages. *The American Statistician* 50(4):361–365. https://doi.org/10.1080/00031305.1996.10473566
6. Lawson C.L., Hanson R.J. 1974. Solving Least Squares Problems, Ch. 23 (NNLS). https://en.wikipedia.org/wiki/Non-negative_least_squares
7. Knuth D.E. 1997. *The Art of Computer Programming*, Vol. 2: Seminumerical Algorithms, §3.4.1 (Poisson deviate generation). https://en.wikipedia.org/wiki/Poisson_distribution#Random_variate_generation
