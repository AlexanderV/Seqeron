# Evidence Artifact: ONCO-SIG-003

**Test Unit ID:** ONCO-SIG-003
**Algorithm:** Signature Exposure Estimation — Bootstrap Confidence Intervals
**Date Collected:** 2026-06-14

---

## Online Sources

### MSA: reproducible mutational signature attribution with confidence based on simulations (Senkin 2021)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed: *BMC Bioinformatics* 22:540)

**Key Extracted Points:**

1. **Resampling distribution:** The method uses a parametric bootstrap in which "mutations are accumulated following Poisson distributions for each mutation class, such as a specific trinucleotide context." The catalog is resampled with total mutational burden M and per-channel probabilities pᵢ = catalogᵢ / M (normalized mutation counts). (Retrieved by fetching the PMC article and asking for the bootstrap methodology; the page returned the quoted sentence and the M / pᵢ parameterization.)
2. **Refit per replicate:** "For each bootstrap sample, NNLS attribution is applied to derive the vector of signature activities." — the same NNLS fit is re-run on every resampled catalog.
3. **Confidence interval:** "95% confidence intervals are then derived for each signature attribution by taking [2.5%, 97.5%] percentiles of the resulting bootstrap activities." — the percentile method.
4. **Replicate count:** "Simulated sample with 1000 bootstrap variations" (Figure 2 caption) — 1000 replicates is the standard count used.

### sigminer `sig_fit_bootstrap` (reference implementation; Wang S. et al.)

**URL (doc):** https://shixiangwang.github.io/sigminer-doc/sigfit.html
**URL (source):** https://raw.githubusercontent.com/ShixiangWang/sigminer/master/R/sig_fit_bootstrap.R
**Accessed:** 2026-06-14
**Authority rank:** 3 (established bioinformatics package; cites Huang et al. 2018 as the primary method)

**Key Extracted Points:**

1. **Multinomial resample (source code):** the catalog is resampled by
   `sampled <- sample(seq(K), total_count, replace = TRUE, prob = catalog / sum(catalog))` then tabulated into a resampled catalog vector. Drawing `total_count = Σ catalog` indices with replacement weighted by `catalog/Σcatalog` is exactly a multinomial(N, p) draw with N = Σ catalog and pₖ = catalogₖ / N. (Retrieved by fetching the raw `sig_fit_bootstrap.R`; the quoted two lines were returned.)
2. **Refit per replicate:** the function "uses the resampling data of original input and runs `sig_fit()` multiple times to estimate the exposure."
3. **Replicate count:** the documentation states "Bootstrap replicates >= 100 is recommended." (Retrieved from the sigminer-doc Chapter 4 page.)
4. **Primary reference:** sigminer cites Huang X, Wojtowicz D, Przytycka TM — "Detecting presence of mutational signatures in cancer with confidence" (Bioinformatics 2018) as the underlying method.

### Huang, Wojtowicz & Przytycka (2018) — Detecting presence of mutational signatures with confidence

**URL:** https://academic.oup.com/bioinformatics/article/34/2/330/4209996 ; PubMed: https://pubmed.ncbi.nlm.nih.gov/29028923/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed: *Bioinformatics* 34(2):330–337)

**Key Extracted Points:**

1. **Bootstrap principle:** the paper "proposed two complementary ways of measuring confidence and stability of decomposition results" and "employed bootstrap methods to evaluate the robustness of results obtained from mutation data decomposition." It is the primary method implemented by the `SignatureEstimation` R package (`bootstrapSigExposures`). (Retrieved via web search of the article landing page / PubMed record; the abstract text was returned.)

### Percentile bootstrap confidence interval (Efron 1979)

**URL:** https://math.mit.edu/~dav/05.dir/class24-prep-a.pdf (MIT 18.05 lecture notes summarizing Efron's percentile method)
**Accessed:** 2026-06-14
**Authority rank:** 1 (statistical method primary: Efron B., 1979, *Annals of Statistics* 7(1):1–26)

**Key Extracted Points:**

1. **Percentile interval definition:** "The percentile bootstrap is derived by using the 2.5 and the 97.5 percentiles of the bootstrap distribution as the 95% confidence interval." More generally, for confidence c the interval is the [100·½(1−c), 100·(1−½(1−c))] percentiles of the bootstrap distribution of the statistic. (Retrieved via web search for the Efron percentile method; the quoted definition was returned.)

### Hyndman & Fan (1996) — sample quantile type 7

**URL:** https://en.wikipedia.org/wiki/Quantile (Estimating quantiles from a sample → Type 7 / "linear interpolation of the modes")
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia, citing the primary Hyndman R.J. & Fan Y., 1996, *The American Statistician* 50(4):361–365)

**Key Extracted Points:**

1. **Type-7 quantile:** the default quantile estimator in R and NumPy places the p-quantile at 0-based rank h = p·(n−1) of the sorted sample and linearly interpolates: Q(p) = x₍⌊h⌋₎ + (h − ⌊h⌋)·(x₍⌊h⌋₊₁₎ − x₍⌊h⌋₎). This is the convention used to realize the percentile interval on the finite set of bootstrap exposures.

---

## Documented Corner Cases and Failure Modes

### From Senkin (2021) / sigminer

1. **Zero-mutation catalog (N = 0):** with no mutations there is nothing to resample; every resampled catalog is all-zero, so every replicate exposure is 0 and the interval collapses to [0, 0].
2. **Degenerate single-channel catalog:** a multinomial draw over a single non-zero channel is deterministic (all N mutations fall in that channel), so the resample equals the observed catalog exactly and every replicate exposure equals the point estimate.

### From bootstrap theory (Efron 1979)

1. **Single replicate (R = 1):** the percentile of a one-element distribution is that single value, so lower = upper = mean = the single replicate's exposure.

---

## Test Datasets

### Dataset: Type-7 percentile worked values (hand-derived from Hyndman & Fan 1996)

**Source:** Hyndman & Fan (1996) type-7 formula Q(p) = x₍⌊h⌋₎ + (h−⌊h⌋)(x₍⌊h⌋₊₁₎ − x₍⌊h⌋₎), h = p·(n−1).

| Sorted sample | p | h = p·(n−1) | Q(p) |
|---------------|------|-------------|------|
| [0,1,2,3,4] (n=5) | 0.5 | 2.0 | 2.0 |
| [0,1,2,3,4] (n=5) | 0.025 | 0.1 | 0 + 0.1·(1−0) = 0.1 |
| [0,1,2,3,4] (n=5) | 0.975 | 3.9 | 3 + 0.9·(4−3) = 3.9 |
| [2,4,6,8] (n=4) | 0.025 | 0.075 | 2 + 0.075·(4−2) = 2.15 |
| [2,4,6,8] (n=4) | 0.975 | 2.925 | 6 + 0.925·(8−6) = 7.85 |
| [2,4,6,8] (n=4) | 0.5 | 1.5 | 4 + 0.5·(6−4) = 5.0 |

### Dataset: Deterministic single-channel bootstrap (multinomial collapse)

**Source:** Senkin (2021) / sigminer multinomial resampling; degenerate single-outcome multinomial.

| Parameter | Value |
|-----------|-------|
| catalog | [10] (one channel, N = 10) |
| signatures | [[1.0]] (one signature, single channel) |
| Every resampled catalog | [10] (multinomial of 10 over one channel is deterministic) |
| NNLS exposure per replicate | 10 (fit of [10] onto [[1.0]]) |
| Point estimate / mean / lower / upper | 10 / 10 / 10 / 10 |

### Dataset: Zero-mutation catalog

**Source:** documented corner case (N = 0).

| Parameter | Value |
|-----------|-------|
| catalog | [0, 0, 0] |
| signatures | [[1,0,0],[0,1,0]] |
| Every replicate exposure | 0 |
| Interval (per signature) | point 0, mean 0, lower 0, upper 0 |

---

## Assumptions

1. **ASSUMPTION: Type-7 quantile convention** — the sources fix the percentile *method* ([2.5%, 97.5%] percentiles) but not the exact interpolation rule used to estimate a percentile from a finite sample. We adopt the type-7 (R/NumPy default; Hyndman & Fan 1996) linear-interpolation estimator. This is correctness-affecting (it sets the exact bound values) but is the documented default of the reference tooling (R `quantile`, used by sigminer), so it is the conventional, source-aligned choice rather than an invented value. Tests pin it with hand-derived type-7 values; the dominant degenerate/deterministic tests are convention-independent.
2. **ASSUMPTION: Fixed RNG seed (42)** — the bootstrap is randomized; sources do not prescribe a seed. A fixed default seed makes results reproducible (matching the repository's Phylogenetics bootstrap convention). Non-correctness-affecting for the deterministic test cases (single-channel collapse, N = 0) whose outcomes are seed-independent.

---

## Recommendations for Test Coverage

1. **MUST Test:** Percentile (type-7) returns the exact hand-derived interpolated values for the [0,1,2,3,4] and [2,4,6,8] samples at p = 0.025, 0.5, 0.975 (and p = 0, 1 → min/max). — Evidence: Hyndman & Fan (1996) type-7.
2. **MUST Test:** Single-channel catalog [10] with signature [[1.0]] yields, for any seed/replicate count, point = mean = lower = upper = 10 (multinomial collapse + NNLS). — Evidence: Senkin (2021) / sigminer multinomial resample; NNLS fit.
3. **MUST Test:** Zero-mutation catalog yields all-zero intervals (point/mean/lower/upper = 0) per signature. — Evidence: documented N = 0 corner case.
4. **MUST Test:** Lower ≤ point estimate region ≤ Upper and Lower ≤ Mean ≤ Upper, all exposures ≥ 0 (invariant), on a non-degenerate two-signature catalog with a fixed seed (reproducible). — Evidence: Efron (1979) percentile interval ordering; NNLS x ≥ 0 (Lawson & Hanson 1974).
5. **MUST Test:** Determinism — two calls with the same seed return identical intervals; the 97.5 bound ≥ the 2.5 bound. — Evidence: fixed-seed reproducibility.
6. **SHOULD Test:** One interval per signature, in signature order; point estimate equals `FitSignatures(observed).Exposures`. — Rationale: contract consistency with ONCO-SIG-002.
7. **MUST Test (failure modes):** null catalog/signatures → ArgumentNullException; ragged/empty signatures, catalog-length mismatch, negative count → ArgumentException; replicates < 1 or confidence ∉ (0,1) → ArgumentOutOfRangeException. — Evidence: input-validation contract mirrored from `FitSignatures`.

---

## References

1. Senkin S. (2021). MSA: reproducible mutational signature attribution with confidence based on simulations. *BMC Bioinformatics* 22:540. https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/
2. Huang X., Wojtowicz D., Przytycka T.M. (2018). Detecting presence of mutational signatures in cancer with confidence. *Bioinformatics* 34(2):330–337. https://academic.oup.com/bioinformatics/article/34/2/330/4209996
3. Wang S. et al. sigminer: `sig_fit_bootstrap`. https://shixiangwang.github.io/sigminer-doc/sigfit.html (source: https://raw.githubusercontent.com/ShixiangWang/sigminer/master/R/sig_fit_bootstrap.R)
4. Efron B. (1979). Bootstrap Methods: Another Look at the Jackknife. *Annals of Statistics* 7(1):1–26. https://doi.org/10.1214/aos/1176344552 (percentile method summarized at https://math.mit.edu/~dav/05.dir/class24-prep-a.pdf)
5. Hyndman R.J., Fan Y. (1996). Sample Quantiles in Statistical Packages. *The American Statistician* 50(4):361–365. https://doi.org/10.1080/00031305.1996.10473566 (type-7 default; https://en.wikipedia.org/wiki/Quantile)
6. Lawson C.L., Hanson R.J. (1974). Solving Least Squares Problems, Ch. 23 (NNLS). https://en.wikipedia.org/wiki/Non-negative_least_squares

---

## Change History

- **2026-06-14**: Initial documentation.
