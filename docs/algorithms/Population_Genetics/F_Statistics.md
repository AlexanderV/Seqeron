| Field | Value |
|-------|-------|
| Algorithm Group | Population Genetics |
| Test Unit ID | POP-FST-001 |
| Related Projects | Seqeron |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

This document covers the F-statistics suite $F_{ST}$, $F_{IS}$, and $F_{IT}$ as used by the repository's population-genetics code. These statistics describe population differentiation and inbreeding by comparing allele-frequency variance or observed versus expected heterozygosity across populations (Wright 1950; Wright 1965; Holsinger and Weir 2009). In this repository, `CalculateFst` implements a two-population Wright variance-based $F_{ST}$ calculation, while `CalculateFStatistics` returns heterozygosity-based $F_{IS}$, $F_{IT}$, and $F_{ST}$ components from per-locus summaries. `CalculatePairwiseFst` lifts the pairwise $F_{ST}$ calculation to a symmetric all-pairs matrix.

## 2. Scientific / Formal Basis

> A = F_ST, B = F_IS, C = F_IT

### 2.A F_ST

#### Domain Context

$F_{ST}$ measures population differentiation attributable to subpopulation structure. In the standard interpretation cited by the repository evidence, $F_{ST} = 0$ corresponds to no differentiation and $F_{ST} = 1$ corresponds to complete differentiation (Wright 1965; Wikipedia "Fixation index").

#### Core Model

Wright's variance-based definition is:

$$
F_{ST} = \frac{\sigma_S^2}{\bar{p}(1-\bar{p})}
$$

where $\bar{p}$ is the population-size-weighted mean allele frequency and $\sigma_S^2$ is the weighted variance among subpopulations (Wright 1965; Wikipedia "Fixation index"). For two populations with sizes $n_1, n_2$ and allele frequencies $p_1, p_2$:

- $\bar{p} = (n_1 p_1 + n_2 p_2) / (n_1 + n_2)$
- $\sigma_S^2 = \bigl(n_1 (p_1 - \bar{p})^2 + n_2 (p_2 - \bar{p})^2\bigr) / (n_1 + n_2)$

When $F_{ST}$ is expressed alongside $F_{IS}$ and $F_{IT}$, the heterozygosity form is:

$$
F_{ST} = 1 - \frac{H_S}{H_T}
$$

where $H_S$ is expected heterozygosity within subpopulations and $H_T$ is expected heterozygosity in the total population (Wikipedia "F-statistics").

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-FST-01 | $0 \leq F_{ST} \leq 1$ | The repository evidence and tests treat Wright's variance-based $F_{ST}$ in this range (Wikipedia "Fixation index"; POP-FST-001 evidence and tests). |
| INV-FST-02 | If the compared populations have identical allele frequencies at every compared locus, then $F_{ST} = 0$ | The among-population variance term $\sigma_S^2$ is zero when each compared locus has the same allele frequency in both populations (Wright 1965; POP-FST-001 tests). |
| INV-FST-03 | For equal-size two-population comparisons where every compared locus has $p_1 = 1$ and $p_2 = 0$, $F_{ST} = 1$ | Each locus contributes $\bar{p} = 0.5$, $\sigma_S^2 = 0.25$, and $\bar{p}(1-\bar{p}) = 0.25$, so the ratio is exactly $1$ (Wright 1965; POP-FST-001 tests). |

#### Comparison with Related Methods

| Aspect | Wright variance-based $F_{ST}$ | Weir-Cockerham $\theta$ |
|--------|-------------------------------|--------------------------|
| Definition basis | Ratio of among-subpopulation variance to total expected heterozygosity (Wright 1965; Wikipedia "Fixation index") | ANOVA variance components with finite-sample correction (Weir and Cockerham 1984; Holsinger and Weir 2009) |
| Statistical role | Direct population-parameter formulation from allele frequencies and population-size weights | Finite-sample estimator for sampled population data |

### 2.B F_IS

#### Domain Context

$F_{IS}$ measures the inbreeding coefficient of an individual relative to its subpopulation by comparing observed heterozygosity with the subpopulation expectation (Wikipedia "F-statistics").

#### Core Model

The heterozygosity-based definition is:

$$
F_{IS} = 1 - \frac{H_I}{H_S}
$$

where $H_I$ is observed heterozygosity and $H_S$ is expected heterozygosity within subpopulations (Wikipedia "F-statistics").

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-FIS-01 | $-1 \leq F_{IS} \leq 1$ | The repository tests assert this range and the evidence file documents that $F_{IS}$ can be negative under excess heterozygosity (Wikipedia "F-statistics"; POP-FST-001 evidence and tests). |
| INV-FIS-02 | $F_{IS}$ may be negative when observed heterozygosity exceeds the within-subpopulation expectation | The definition $1 - H_I / H_S$ becomes negative whenever $H_I > H_S$ (Wikipedia "F-statistics"; POP-FST-001 tests). |

### 2.C F_IT

#### Domain Context

$F_{IT}$ measures the inbreeding coefficient of an individual relative to the total population by comparing observed heterozygosity with the total-population expectation (Wikipedia "F-statistics").

#### Core Model

The heterozygosity-based definition is:

$$
F_{IT} = 1 - \frac{H_I}{H_T}
$$

where $H_I$ is observed heterozygosity and $H_T$ is expected heterozygosity in the total population (Wikipedia "F-statistics"). The partition identity relating the suite is:

$$
(1 - F_{IT}) = (1 - F_{IS})(1 - F_{ST})
$$

(Wright 1965; Wikipedia "F-statistics").

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-FIT-01 | $-1 \leq F_{IT} \leq 1$ | The POP-FST-001 tests assert this range for the heterozygosity-based implementation, consistent with the cited definition (Wikipedia "F-statistics"). |
| INV-FIT-02 | $(1 - F_{IT}) = (1 - F_{IS})(1 - F_{ST})$ | Wright's partition identity is algebraic for the heterozygosity definitions and is tested as an exact equality in POP-FST-001. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[CalculateFst] population1` | `IEnumerable<(double AlleleFreq, int SampleSize)>` | required | Per-locus allele frequency and sample size for the first population | Compared locus-by-locus with `population2`; empty input is accepted |
| `[CalculateFst] population2` | `IEnumerable<(double AlleleFreq, int SampleSize)>` | required | Per-locus allele frequency and sample size for the second population | Compared locus-by-locus with `population1`; empty input is accepted |
| `[CalculatePairwiseFst] populations` | `IEnumerable<(string PopulationId, IReadOnlyList<(double AlleleFreq, int SampleSize)> Variants)>` | required | Population IDs and per-locus allele-frequency lists for pairwise matrix construction | Matrix size equals the number of supplied populations |
| `[CalculateFStatistics] pop1Name` | `string` | required | Label copied into the returned `FStatistics` record | Preserved verbatim in the result |
| `[CalculateFStatistics] pop2Name` | `string` | required | Label copied into the returned `FStatistics` record | Preserved verbatim in the result |
| `[CalculateFStatistics] variantData` | `IEnumerable<(int HetObs1, int N1, int HetObs2, int N2, double AlleleFreq1, double AlleleFreq2)>` | required | Per-locus observed heterozygote counts, sample sizes, and allele frequencies for two populations | Empty input is accepted |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[CalculateFst] return value` | `double` | Pairwise Wright variance-based $F_{ST}$ for the compared populations |
| `[CalculatePairwiseFst] return value` | `double[,]` | Symmetric matrix of pairwise $F_{ST}$ values for the supplied populations |
| `[CalculateFStatistics] FStatistics.Fst` | `double` | Heterozygosity-based $F_{ST}$ component |
| `[CalculateFStatistics] FStatistics.Fis` | `double` | Heterozygosity-based $F_{IS}$ component |
| `[CalculateFStatistics] FStatistics.Fit` | `double` | Heterozygosity-based $F_{IT}$ component |
| `[CalculateFStatistics] FStatistics.Population1` | `string` | First population label copied from `pop1Name` |
| `[CalculateFStatistics] FStatistics.Population2` | `string` | Second population label copied from `pop2Name` |

### 3.3 Preconditions and Validation

The documented methods operate on typed tuple inputs and do not add explicit domain validation beyond what is visible in the source. `CalculateFst` returns `0` when either compared population sequence is empty, `CalculateFStatistics` returns an all-zero `FStatistics` record when `variantData` is empty, and `CalculatePairwiseFst` returns an empty matrix when no populations are supplied. The source does not perform explicit range checks on allele frequencies, sample sizes, or locus alignment; the repository-specific consequences of that choice are documented in Section 5.

## 4. Algorithm

### 4.A F_ST

#### High-Level Steps

1. For each compared locus, compute the population-size-weighted mean allele frequency $\bar{p}$.
2. Compute the weighted among-population variance $\sigma_S^2$ for that locus.
3. Compute the expected heterozygosity contribution $\bar{p}(1-\bar{p})$ for that locus.
4. Sum the variance terms and the heterozygosity terms across loci.
5. Return the ratio of the summed variance to the summed heterozygosity (Wright 1965; Wikipedia "Fixation index").

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Pairwise scalar $F_{ST}$ across $L$ compared loci | O(L) | O(1) | One weighted mean, variance, and heterozygosity calculation per compared locus |
| Pairwise $F_{ST}$ matrix across $k$ populations and $L$ compared loci per pair | O(k^2 x L) | O(k^2) | One scalar $F_{ST}$ calculation per unordered population pair |

### 4.B F_IS

#### High-Level Steps

1. Aggregate observed heterozygosity across the provided loci.
2. Aggregate expected heterozygosity within each subpopulation using $2p(1-p)$ for each population and locus.
3. Normalize both aggregates by the total sample count to obtain $H_I$ and $H_S$.
4. Return $1 - H_I/H_S$ (Wikipedia "F-statistics").

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Heterozygosity-based $F_{IS}$ across $L$ loci | O(L) | O(1) | One pass over the supplied locus summaries |

### 4.C F_IT

#### High-Level Steps

1. Aggregate observed heterozygosity across the provided loci.
2. Compute the total-population allele frequency $\bar{p}$ for each locus from the two population frequencies and sample sizes.
3. Aggregate total expected heterozygosity using $2\bar{p}(1-\bar{p})$.
4. Normalize by the total sample count to obtain $H_I$ and $H_T$.
5. Return $1 - H_I/H_T$ and preserve the partition identity with $F_{IS}$ and $F_{ST}$ (Wright 1965; Wikipedia "F-statistics").

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Heterozygosity-based $F_{IT}$ across $L$ loci | O(L) | O(1) | One pass over the supplied locus summaries |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)

- `PopulationGeneticsAnalyzer.CalculateFst(...)`: computes Wright's variance-based pairwise $F_{ST}$ from two per-locus allele-frequency sequences.
- `PopulationGeneticsAnalyzer.CalculatePairwiseFst(...)`: builds a symmetric pairwise $F_{ST}$ matrix by reusing `CalculateFst(...)`.
- `PopulationGeneticsAnalyzer.CalculateFStatistics(...)`: aggregates observed and expected heterozygosity and returns the `FStatistics` record.
- `PopulationGeneticsAnalyzer.FStatistics`: record type that carries `Fst`, `Fis`, `Fit`, `Population1`, and `Population2`.

### 5.2 Current Behavior

- `CalculateFst(...)` materializes both input sequences, returns `0` when either sequence is empty, iterates only across the first `Math.Min(pop1.Count, pop2.Count)` loci, and returns `0` when the accumulated heterozygosity denominator is `0`.
- `CalculatePairwiseFst(...)` materializes the population list, allocates an `n x n` matrix, computes only the upper triangle by calling `CalculateFst(...)`, mirrors each value into the lower triangle, and leaves the diagonal at the array default `0`.
- `CalculateFStatistics(...)` returns `new FStatistics(0, 0, 0, pop1Name, pop2Name)` for empty input and otherwise accumulates observed heterozygosity, within-population expected heterozygosity, and total expected heterozygosity before computing `Fis`, `Fit`, and `Fst`.
- None of the three methods performs explicit range validation on allele frequencies or sample sizes in the supplied tuples.

### 5.3 Conformance to Theory / Spec

#### 5.3.A F_ST

**Implemented (verbatim from the cited theory/spec):**

- Wright's variance-based pairwise $F_{ST}$ formula using population-size-weighted mean allele frequency and weighted among-population variance.
- Ratio-of-sums aggregation across compared loci in `CalculateFst(...)`.
- Symmetric pairwise matrix construction in `CalculatePairwiseFst(...)`, with `Fst(i, i) = 0` from the matrix initialization and mirrored off-diagonal entries.

**Intentionally simplified:**

- Weir-Cockerham $\theta$ and other finite-sample ANOVA-based corrections are not implemented; **consequence:** this repository's standalone $F_{ST}$ method returns Wright's direct variance-based value, which can differ from sample-corrected estimators discussed by Weir and Cockerham (1984).

**Not implemented:**

- (none)

#### 5.3.B F_IS

**Implemented (verbatim from the cited theory/spec):**

- Heterozygosity-based $F_{IS} = 1 - H_I/H_S$ in `CalculateFStatistics(...)`.
- Negative $F_{IS}$ values when observed heterozygosity exceeds the within-subpopulation expectation.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

#### 5.3.C F_IT

**Implemented (verbatim from the cited theory/spec):**

- Heterozygosity-based $F_{IT} = 1 - H_I/H_T$ in `CalculateFStatistics(...)`.
- Exact partition identity `(1 - Fit) = (1 - Fis)(1 - Fst)` for the heterozygosity-based computation used by `CalculateFStatistics(...)`.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `CalculateFst(...)` truncates to the shared prefix length of the two input sequences | Assumption | Extra loci in the longer input are ignored rather than rejected, so callers must align locus lists before comparing populations | accepted | Implemented as `for (int i = 0; i < Math.Min(pop1.Count, pop2.Count); i++)`; inherited by `CalculatePairwiseFst(...)` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `CalculateFst(...)` with either population sequence empty | Returns `0` | The method checks the materialized counts and returns `0` before any per-locus calculation |
| Both compared populations fixed for the same allele at every compared locus | Returns `0` | The accumulated heterozygosity denominator is `0`, and the implementation returns `0` instead of dividing by zero |
| Identical compared populations | $F_{ST} = 0$ | The weighted among-population variance is zero at every compared locus |
| Equal-size fixed opposite alleles at every compared locus | $F_{ST} = 1$ | Each locus contributes the same numerator and denominator value, so the ratio is exactly `1` |
| Excess heterozygosity in `CalculateFStatistics(...)` | `Fis` may be negative | The heterozygosity definition yields negative `Fis` when `H_I > H_S`, which is covered by the tests |

### 6.2 Limitations

The repository's standalone `CalculateFst(...)` method implements Wright's variance-based two-population $F_{ST}$ rather than Weir-Cockerham $\theta$ or other finite-sample estimators (Wright 1965; Weir and Cockerham 1984; Holsinger and Weir 2009). `CalculateFst(...)` and `CalculatePairwiseFst(...)` rely on enumeration order rather than explicit locus identifiers, so misaligned input sequences can silently compare different loci because only the shared prefix length is processed. The evidence file also records the standard note that $F_{ST}$ is not a metric distance, because it does not satisfy the triangle inequality (Wikipedia "Fixation index").

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PopulationGeneticsAnalyzer_FStatistics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_FStatistics_Tests.cs) - covers the exact-value checks, range checks, pairwise matrix symmetry, monomorphic cases, and negative `Fis` scenarios documented here.
- Test spec: [POP-FST-001.md](../../../tests/TestSpecs/POP-FST-001.md) - defines the POP-FST-001 contract and expected invariants for the suite.
- Evidence: [POP-FST-001-Evidence.md](../../../docs/Evidence/POP-FST-001-Evidence.md) - collects the literature and derivations used by the tests.
- Related algorithms: [Allele_Frequency.md](Allele_Frequency.md), [Hardy_Weinberg_Test.md](Hardy_Weinberg_Test.md), [Diversity_Statistics.md](Diversity_Statistics.md)

## 8. References

1. Wright, S. 1950. Genetical structure of populations. Nature 166. https://doi.org/10.1038/166247a0
2. Wright, S. 1965. The interpretation of population structure by F-statistics with special regard to systems of mating. Evolution 19. https://doi.org/10.2307/2406450
3. Weir, B. S., and C. C. Cockerham. 1984. Estimating F-statistics for the analysis of population structure. Evolution 38. https://doi.org/10.2307/2408641
4. Holsinger, K. E., and B. S. Weir. 2009. Genetics in geographically structured populations: defining, estimating and interpreting FST. Nature Reviews Genetics 10. https://doi.org/10.1038/nrg2611
5. Wikipedia contributors. 2026. Fixation index. Wikipedia. https://en.wikipedia.org/wiki/Fixation_index
6. Wikipedia contributors. 2026. F-statistics. Wikipedia. https://en.wikipedia.org/wiki/F-statistics
