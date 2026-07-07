# Allele Frequency

| Field | Value |
|-------|-------|
| Algorithm Group | Population Genetics |
| Test Unit ID | POP-FREQ-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Allele frequency measures the relative abundance of an allele at one genetic locus within a population.[1][2] In the biallelic diploid case documented here, the frequencies are obtained by counting the allele copies contributed by homozygous and heterozygous genotypes.[1][4] This file also covers minor allele frequency (MAF), defined as the smaller of the two allele frequencies, and the repository helper that filters variants by an inclusive MAF interval.[2][5] The repository implementation provides deterministic algebraic summaries rather than a heuristic or probabilistic model.[5][6][7]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Allele frequency, also called gene frequency, is a basic population-genetics quantity used to describe genetic variation at a locus.[1][2] For a diploid biallelic locus, the two allele frequencies summarize how often each allele copy appears across all sampled chromosomes.[1] Minor allele frequency is often used to distinguish balanced polymorphisms from rare or nearly fixed variants in downstream analyses.[2]

### 2.2 Core Model

For genotype counts $n_{AA}$, $n_{Aa}$, and $n_{aa}$ in a diploid population, the allele frequencies are:[1][4]

$$
p = \frac{2n_{AA} + n_{Aa}}{2(n_{AA} + n_{Aa} + n_{aa})}
$$

$$
q = \frac{2n_{aa} + n_{Aa}}{2(n_{AA} + n_{Aa} + n_{aa})}
$$

The minor allele frequency is the smaller of the two allele frequencies:[2]

$$
MAF = \min(p, q) = \min(f, 1 - f)
$$

The original document also recorded the common population-genetics convention that variants with $MAF \ge 0.05$ are often treated as common and those with $MAF < 0.05$ as rarer variants.[2]

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | $p + q = 1$ for any non-empty biallelic diploid sample | Every sampled chromosome copy is counted as either the first or the second allele.[1][4] |
| INV-02 | $0 \le p \le 1$ and $0 \le q \le 1$ | Each frequency is a count divided by the total allele count.[1][2] |
| INV-03 | $0 \le MAF \le 0.5$ | `MAF` is defined as the smaller of two complementary frequencies.[2][5] |
| INV-04 | Major-allele copies plus minor-allele copies equals the total allele count | Diploid allele accounting counts two alleles per individual.[1][4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[CalculateAlleleFrequencies] homozygousMajor` | `int` | required | Count of major-major genotypes | Must be non-negative; negative input throws `ArgumentOutOfRangeException`.[5][6][7] |
| `[CalculateAlleleFrequencies] heterozygous` | `int` | required | Count of heterozygous genotypes | Must be non-negative; negative input throws `ArgumentOutOfRangeException`.[5][6][7] |
| `[CalculateAlleleFrequencies] homozygousMinor` | `int` | required | Count of minor-minor genotypes | Must be non-negative; negative input throws `ArgumentOutOfRangeException`.[5][6][7] |
| `[CalculateMAF] genotypes` | `IEnumerable<int>` | required | Diploid genotype codes | The documented convention is `0 = hom ref`, `1 = het`, `2 = hom alt`.[5][6][7] |
| `[FilterByMAF] variants` | `IEnumerable<Variant>` | required | Variants to test against the MAF interval | Each `Variant` must already carry an `AlleleFrequency` value used by the filter.[6][7] |
| `[FilterByMAF] minMAF` | `double` | `0.01` | Lower inclusive MAF threshold | Compared with `>=` in the current implementation.[7] |
| `[FilterByMAF] maxMAF` | `double` | `0.5` | Upper inclusive MAF threshold | Compared with `<=` in the current implementation.[7] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[CalculateAlleleFrequencies] MajorFreq` | `double` | Frequency of the allele counted by `2 * homozygousMajor + heterozygous`.[1][7] |
| `[CalculateAlleleFrequencies] MinorFreq` | `double` | Frequency of the allele counted by `2 * homozygousMinor + heterozygous`.[1][7] |
| `[CalculateMAF] return value` | `double` | Minor allele frequency in the range `[0, 0.5]` for the documented genotype encoding.[2][5][7] |
| `[FilterByMAF] return value` | `IEnumerable<Variant>` | Lazy sequence of variants whose computed `MAF = min(AlleleFrequency, 1 - AlleleFrequency)` lies inside `[minMAF, maxMAF]`.[5][7] |

### 3.3 Preconditions and Validation

`CalculateAlleleFrequencies` rejects negative genotype counts with `ArgumentOutOfRangeException` and returns `(0, 0)` when the total allele count is zero.[5][6][7] `CalculateMAF` returns `0` for an empty genotype collection.[5][6][7] `FilterByMAF` uses inclusive lower and upper bounds, preserves the input order because it yields matching items as it iterates, and does not recompute allele frequency from genotype calls.[5][6][7] The repository code does not perform additional validation on genotype codes passed to `CalculateMAF` or on `Variant.AlleleFrequency` values passed to `FilterByMAF`.[7]

## 4. Algorithm

### 4.1 High-Level Steps

1. For `CalculateAlleleFrequencies`, validate that all three genotype counts are non-negative.
2. Compute the total number of allele copies as `2 * (homozygousMajor + heterozygous + homozygousMinor)`.
3. Return `(0, 0)` when that total is zero; otherwise divide the counted major and minor allele copies by the total.
4. For `CalculateMAF`, materialize the genotype sequence, sum alternate-allele copies from the documented `0/1/2` encoding, convert that sum to an alternate-allele frequency, and return `min(f, 1 - f)`.
5. For `FilterByMAF`, compute `min(AlleleFrequency, 1 - AlleleFrequency)` for each input `Variant` and yield only the records inside the inclusive threshold interval.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository uses the genotype coding documented in both the current source and the tests: `0 = homozygous reference`, `1 = heterozygous`, and `2 = homozygous alternate`.[5][6][7] `FilterByMAF` does not inspect genotype calls; it treats `Variant.AlleleFrequency` as the source value and converts that stored allele frequency to `MAF` with `min(f, 1 - f)`.[7]

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateAlleleFrequencies` | `O(1)` | `O(1)` | Constant-time arithmetic on three counts.[5][7] |
| `CalculateMAF` | `O(n)` | `O(n)` | Materializes the genotype sequence with `ToList()` before summing genotype codes.[5][7] |
| `FilterByMAF` | `O(n)` | `O(1)` incremental | Iterator over the input sequence; output size depends on the number of retained variants.[5][7] |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)

- `PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(int, int, int)`: Validates non-negative counts and returns a `(MajorFreq, MinorFreq)` tuple.
- `PopulationGeneticsAnalyzer.CalculateMAF(IEnumerable<int>)`: Computes `MAF` from genotype-coded alternate-allele counts.
- `PopulationGeneticsAnalyzer.FilterByMAF(IEnumerable<PopulationGeneticsAnalyzer.Variant>, double, double)`: Lazily filters `Variant` records by an inclusive MAF interval.
- `PopulationGeneticsAnalyzer.Variant`: Stores the `AlleleFrequency` value consumed by `FilterByMAF`.

### 5.2 Current Behavior

The repository returns allele frequencies in `MajorFreq` / `MinorFreq` order rather than using reference/alternate naming.[7] `CalculateMAF` materializes the input sequence and sums the raw genotype integers, so the method assumes the documented `0/1/2` encoding rather than validating each element.[7] `FilterByMAF` is implemented with `yield return`, which preserves input order and evaluates items lazily as callers enumerate the result.[5][7]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Biallelic diploid allele-frequency calculation from genotype counts using `2n_hom + n_het` allele accounting.[1][4][7]
- Minor allele frequency as `min(f, 1 - f)`.[2][5][7]

**Intentionally simplified:**

- `CalculateMAF` accepts pre-encoded diploid genotype integers instead of parsing genotype strings or phased haplotypes; **consequence:** callers must normalize genotype representation before using the helper.
- `FilterByMAF` operates on the stored `Variant.AlleleFrequency` field only; **consequence:** genotype-level allele counting is outside the filtering API.

**Not implemented:**

- Multi-allelic allele-frequency and MAF calculation; **users should rely on:** no current alternative.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Zero genotype counts in `CalculateAlleleFrequencies` | Returns `(0, 0)` | The implementation special-cases zero total allele count.[5][7] |
| All homozygous major | Returns `(1.0, 0.0)` | All allele copies belong to the first allele.[5][6] |
| All homozygous minor | Returns `(0.0, 1.0)` | All allele copies belong to the second allele.[5][6] |
| All heterozygous | Returns `(0.5, 0.5)` | Each individual contributes one copy of each allele.[5][6] |
| Empty genotype list in `CalculateMAF` | Returns `0` | The implementation returns zero for empty input.[5][6][7] |
| Monomorphic genotype list | Returns `0` for `MAF` | A fixed locus has no minor allele frequency.[2][5][6] |
| Empty variant list in `FilterByMAF` | Returns an empty sequence | The iterator simply yields no results.[5][6][7] |

### 6.2 Limitations

The documented formulas and the repository helpers are restricted to diploid biallelic data.[1][2][7] `CalculateMAF` assumes `0/1/2` genotype encoding, and `FilterByMAF` assumes each `Variant` already carries a meaningful allele-frequency estimate.[6][7]

## 7. Examples and Related Material

### 7.1 Worked Example

Using the four-o'clock plant example from the original document and the cited genotype-frequency reference, a population with `49` `AA`, `42` `Aa`, and `9` `aa` individuals has `200` total allele copies, `140` copies of `A`, and `60` copies of `a`.[4][6]

$$
p(A) = \frac{2 \times 49 + 42}{200} = \frac{140}{200} = 0.70
$$

$$
q(a) = \frac{2 \times 9 + 42}{200} = \frac{60}{200} = 0.30
$$

So the minor allele frequency is `0.30`.[2][4]

### 7.2 Applications and Use Cases

Minor allele frequency is commonly used to distinguish more common from rarer variants, to filter loci before downstream association analyses, and to summarize per-locus diversity in population datasets.[2][4]

## 8. References

1. Gillespie, J. H. 2004. Population Genetics: A Concise Guide. 2nd ed. Johns Hopkins University Press.
2. Wikipedia contributors. Allele frequency. https://en.wikipedia.org/wiki/Allele_frequency
3. Wikipedia contributors. Minor allele frequency. https://en.wikipedia.org/wiki/Minor_allele_frequency
4. Wikipedia contributors. Genotype frequency. https://en.wikipedia.org/wiki/Genotype_frequency
5. [POP-FREQ-001.md](../../../tests/TestSpecs/POP-FREQ-001.md)
6. [PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs)
7. [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)
