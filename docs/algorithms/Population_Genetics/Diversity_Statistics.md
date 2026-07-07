# Diversity Statistics

| Field | Value |
|-------|-------|
| Algorithm Group | Population Genetics |
| Test Unit ID | POP-DIV-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Diversity statistics summarize how much sequence variation is present within a sample of aligned sequences.[1][2][3] This document covers nucleotide diversity (`pi`), Watterson's theta (`theta_W`), Tajima's D, and the two heterozygosity fields returned by the repository's `DiversityStatistics` record.[1][2][3][4] In the repository, `PopulationGeneticsAnalyzer.CalculateDiversityStatistics(...)` computes all of these deterministic summary statistics from one aligned sequence set.[6][7][8]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

These statistics quantify different aspects of within-population sequence variation and are commonly used to study mutation, neutrality, and population history.[1][2][3] Nucleotide diversity averages pairwise sequence differences, Watterson's theta rescales the number of segregating sites by the harmonic number `a_n`, and Tajima's D compares the two estimators to detect departures from the neutral model.[1][2][3]

### 2.2 Core Model

Nucleotide diversity is the average number of nucleotide differences per site between all unordered sequence pairs:[1][5]

$$
\pi = \frac{\sum_{i<j} d_{ij}}{\binom{n}{2} L}
$$

where `d_ij` is the number of differing positions between sequences `i` and `j`, `n` is the number of sequences, and `L` is sequence length.[1][5]

Watterson's theta is the per-site estimator based on the number of segregating sites `S`:[2][5]

$$
	heta_W = \frac{S}{a_n L}
$$

with

$$
a_n = \sum_{i=1}^{n-1} \frac{1}{i}
$$

The unnormalized Watterson estimate used inside Tajima's D is `S / a_1`.[2][3][5]

Tajima's D compares the average number of pairwise differences `\hat{k}` with the Watterson estimate:[3][5]

$$
D = \frac{\hat{k} - \frac{S}{a_1}}{\sqrt{e_1 S + e_2 S(S-1)}}
$$

where `\hat{k}` is the average number of pairwise differences per pair of sequences, not the per-site value `\pi`.[3][5][6][8] The variance terms use the standard constants

$$
a_1 = \sum_{i=1}^{n-1} \frac{1}{i}, \quad a_2 = \sum_{i=1}^{n-1} \frac{1}{i^2}
$$

$$
b_1 = \frac{n+1}{3(n-1)}, \quad b_2 = \frac{2(n^2+n+3)}{9n(n-1)}
$$

$$
c_1 = b_1 - \frac{1}{a_1}, \quad c_2 = b_2 - \frac{n+2}{a_1 n} + \frac{a_2}{a_1^2}
$$

$$
e_1 = \frac{c_1}{a_1}, \quad e_2 = \frac{c_2}{a_1^2 + a_2}
$$

Negative Tajima's D indicates an excess of rare alleles relative to the neutral expectation, whereas positive Tajima's D indicates a deficit of rare alleles.[3][5][6] The original document also noted the common rule of thumb that values beyond approximately `+/-2` are often treated as notable departures from neutrality.[5]

Expected heterozygosity, or basic gene diversity, is:[4][5]

$$
H_e = 1 - \sum_{i=1}^{m} f_i^2
$$

Nei's unbiased gene diversity rescales that quantity by `n/(n-1)`:[4]

$$
\hat{H} = \frac{n}{n-1}\left(1 - \sum_i p_i^2\right)
$$

For haploid sequence data, the cited current document and test specification treat Nei's unbiased per-site gene diversity as mathematically equivalent to nucleotide diversity `pi`.[4][6]

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `pi >= 0` and `theta_W >= 0` | Both are normalized counts of differences or segregating sites.[1][2][5] |
| INV-02 | `pi = 0` for identical sequences | Every pairwise difference count `d_ij` is zero.[1][6] |
| INV-03 | `theta_W = 0` when `S = 0` | The numerator in Watterson's estimator is the number of segregating sites.[2][6] |
| INV-04 | `D = 0` when `\hat{k} = S / a_1` | Tajima's numerator is the difference between those two estimators.[3][5][6] |
| INV-05 | `\hat{H} = \frac{n}{n-1} H_e` for haploid per-site allele frequencies when `n > 1` | This is Nei's unbiased correction.[4][6] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[CalculateNucleotideDiversity] sequences` | `IEnumerable<IReadOnlyList<char>>` | required | Collection of aligned sequences compared position by position | Meaningful use requires sequences of common length.[1][6][8] |
| `[CalculateWattersonTheta] segregatingSites` | `int` | required | Number of polymorphic positions `S` | `S = 0` returns `0`.[2][6][8] |
| `[CalculateWattersonTheta] sampleSize` | `int` | required | Number of sampled sequences `n` | `n < 2` returns `0`.[2][6][8] |
| `[CalculateWattersonTheta] sequenceLength` | `int` | required | Length `L` used for per-site normalization | `L <= 0` returns `0`.[2][6][8] |
| `[CalculateTajimasD] averagePairwiseDifferences` | `double` | required | `\hat{k}`, the average number of pairwise differences per sequence pair | This is the unnormalized quantity, not per-site `pi`.[3][5][6][8] |
| `[CalculateTajimasD] segregatingSites` | `int` | required | Number of segregating sites `S` | `S = 0` returns `0`.[3][6][8] |
| `[CalculateTajimasD] sampleSize` | `int` | required | Number of sampled sequences `n` | `n < 3` returns `0`.[3][6][8] |
| `[CalculateDiversityStatistics] sequences` | `IEnumerable<IReadOnlyList<char>>` | required | Sequence collection used to compute the full `DiversityStatistics` record | Fewer than two sequences return an all-zero statistics record except for `SampleSize`.[6][8] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[CalculateNucleotideDiversity] return value` | `double` | Per-site nucleotide diversity `pi`.[1][8] |
| `[CalculateWattersonTheta] return value` | `double` | Per-site Watterson estimator `theta_W`.[2][8] |
| `[CalculateTajimasD] return value` | `double` | Tajima's D statistic computed from `\hat{k}`, `S`, and `n`.[3][8] |
| `DiversityStatistics.NucleotideDiversity` | `double` | `pi`, the per-site average pairwise difference.[1][8] |
| `DiversityStatistics.WattersonTheta` | `double` | Per-site Watterson estimator.[2][8] |
| `DiversityStatistics.TajimasD` | `double` | Tajima's D computed from unnormalized `\hat{k}`.[3][8] |
| `DiversityStatistics.SegregratingSites` | `int` | Count of polymorphic positions; the field name retains the source-code typo for API compatibility.[6][8] |
| `DiversityStatistics.SampleSize` | `int` | Number of input sequences.[6][8] |
| `DiversityStatistics.HeterozygosityObserved` | `double` | Nei's unbiased per-site gene diversity.[4][6][8] |
| `DiversityStatistics.HeterozygosityExpected` | `double` | Basic per-site gene diversity `1 - sum p_i^2`.[4][5][8] |

### 3.3 Preconditions and Validation

`CalculateNucleotideDiversity` returns `0` when fewer than two sequences are supplied.[6][7][8] `CalculateWattersonTheta` returns `0` when `sampleSize < 2` or `sequenceLength <= 0`.[6][7][8] `CalculateTajimasD` returns `0` when `segregatingSites = 0`, when `sampleSize < 3`, or when the variance term in the denominator is non-positive.[6][8] `CalculateDiversityStatistics` returns an all-zero record for `0` or `1` input sequences while preserving the observed `SampleSize`.[6][7][8] The current implementation does not perform an explicit equal-length validation pass, so callers should supply aligned sequences with the same length.[8]

## 4. Algorithm

### 4.1 High-Level Steps

1. Materialize the sequence collection.
2. For `CalculateNucleotideDiversity`, compare every unordered pair of sequences position by position, count mismatches, and normalize by the number of pairs and the sequence length.
3. For `CalculateWattersonTheta`, compute the harmonic number `a_n` and divide `S` by `a_n * L`.
4. For `CalculateTajimasD`, compute `a_1`, `a_2`, the Tajima variance constants, the variance term, and then normalize `\hat{k} - S / a_1`.
5. For `CalculateDiversityStatistics`, count segregating sites, compute `pi`, convert `pi` to `\hat{k} = pi * L`, compute Tajima's D, and then compute the expected and bias-corrected heterozygosity summaries.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository's integrated path follows the test specification exactly: `CalculateDiversityStatistics` first computes per-site `pi`, then converts it to the unnormalized `\hat{k}` required by Tajima's D using `\hat{k} = pi * L`.[6][8] The two heterozygosity outputs are not duplicates: `HeterozygosityExpected` is the basic `1 - sum p_i^2` value, whereas `HeterozygosityObserved` applies Nei's `n/(n-1)` bias correction before averaging across sites.[4][6][8]

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateNucleotideDiversity` | `O(n^2 * m)` | `O(n)` | `n` = number of sequences, `m` = sequence length.[6][8] |
| `CalculateWattersonTheta` | `O(n)` | `O(1)` | Dominated by the harmonic-number loop from `1` to `n - 1`.[6][8] |
| `CalculateTajimasD` | `O(n)` | `O(1)` | Harmonic-number and variance-constant computation over sample size `n`.[6][8] |
| `CalculateDiversityStatistics` | `O(n^2 * m)` | `O(n)` plus output | Uses the pairwise-diversity routine and position-wise scans over the same sequence set.[6][8] |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)

- `PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(IEnumerable<IReadOnlyList<char>>)`: Computes per-site `pi` by explicit pairwise mismatch counting.
- `PopulationGeneticsAnalyzer.CalculateWattersonTheta(int, int, int)`: Computes per-site `theta_W` from `S`, `n`, and `L`.
- `PopulationGeneticsAnalyzer.CalculateTajimasD(double, int, int)`: Computes Tajima's D from unnormalized `\hat{k}`.
- `PopulationGeneticsAnalyzer.CalculateDiversityStatistics(IEnumerable<IReadOnlyList<char>>)`: Returns the full `DiversityStatistics` record.
- `PopulationGeneticsAnalyzer.CalculateObservedHeterozygosity(List<IReadOnlyList<char>>)`: Private helper for Nei's bias-corrected heterozygosity.
- `PopulationGeneticsAnalyzer.CalculateExpectedHeterozygosity(List<IReadOnlyList<char>>)`: Private helper for basic gene diversity.

### 5.2 Current Behavior

`CalculateDiversityStatistics` counts a segregating site by comparing every sequence at one position against the first sequence's character at that position.[8] The method computes Tajima's numerator from `kHat = pi * length`, matching the test specification and correcting the common confusion between per-site `pi` and unnormalized `\hat{k}`.[6][8] The record field is spelled `SegregratingSites` in source and tests, and that typo is part of the public API.[6][7][8] The heterozygosity fields are repository-specific naming for Nei's bias-corrected and uncorrected per-site gene-diversity calculations.[4][6][8]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Per-site nucleotide diversity `pi` as the average pairwise difference rate.[1][8]
- Per-site Watterson's theta `theta_W = S / (a_n L)`.[2][8]
- Tajima's D using the unnormalized `\hat{k}` numerator and the standard variance constants.[3][6][8]
- Expected heterozygosity and Nei's bias-corrected gene diversity averaged across sites.[4][8]

**Intentionally simplified:**

- Sequence comparison is a direct position-by-position character comparison with no ambiguity-code, gap, or missing-data model; **consequence:** callers must pre-normalize aligned sequence data before interpreting these statistics biologically.

**Not implemented:**

- Significance or p-value calculation for Tajima's D; **users should rely on:** no current alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Equal-length aligned sequences are assumed rather than validated explicitly | Assumption | Shorter sequences can produce invalid indexing or uninterpretable site-wise comparisons | accepted | The implementation uses the first sequence length and indexes every other sequence at the same positions.[8] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input | All diversity metrics are `0`, and `SampleSize = 0` | The integrated method special-cases fewer than two sequences.[6][7][8] |
| Single sequence | All diversity metrics are `0`, and `SampleSize = 1` | Pairwise and segregating-site statistics are undefined for `n < 2`.[6][7][8] |
| Identical sequences | `pi = 0`, `theta_W = 0`, `TajimasD = 0`, and `SegregratingSites = 0` | There are no pairwise differences or polymorphic positions.[1][2][6][7] |
| Two sequences that differ at every position | `pi = 1.0` | Every site contributes one mismatch in the only unordered pair.[1][6][7] |
| `sampleSize = 2` for `CalculateTajimasD` | Returns `0` | The repository requires `n >= 3` for Tajima's D.[3][6][8] |

### 6.2 Limitations

These routines assume aligned sequences and compare characters exactly as stored, without modeling ambiguous nucleotides, indels, recombination, or missing data.[1][2][3][8] Tajima's D is reported as a statistic only; the repository does not provide a significance test or p-value layer on top of it.[6][8]

## 7. Examples and Related Material

### 7.1 Worked Example

Using the five-sequence Tajima's D example preserved in the current document and the tests, the aligned sequences have length `20`, `4` segregating sites, and `20` total pairwise differences across `10` unordered pairs.[3][5][6][7] This gives

$$
\hat{k} = \frac{20}{10} = 2.0
$$

$$
\pi = \frac{2.0}{20} = 0.1
$$

$$
	heta_W = \frac{4}{2.0833 \times 20} \approx 0.096
$$

and

$$
D \approx 0.273
$$

which matches both the reference walk-through and the repository tests.[3][5][6][7]

## 8. References

1. Nei, M., and W.-H. Li. 1979. Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases. Proceedings of the National Academy of Sciences 76(10):5269-5273.
2. Watterson, G. A. 1975. On the number of segregating sites in genetical models without recombination. Theoretical Population Biology 7(2):256-276.
3. Tajima, F. 1989. Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics 123(3):585-595.
4. Nei, M. 1978. Estimation of average heterozygosity and genetic distance from a small number of individuals. Genetics 89(3):583-590.
5. Wikipedia contributors. Tajima's D. https://en.wikipedia.org/wiki/Tajima%27s_D
6. [POP-DIV-001.md](../../../tests/TestSpecs/POP-DIV-001.md)
7. [PopulationGeneticsAnalyzer_Diversity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_Diversity_Tests.cs)
8. [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)
