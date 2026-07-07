# Linkage Disequilibrium

| Field | Value |
|-------|-------|
| Algorithm Group | Population Genetics |
| Test Unit ID | POP-LD-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Linkage disequilibrium (LD) measures the non-random association of alleles at different loci within a population.[1][2][5] This file covers the classical disequilibrium coefficient `D`, Lewontin's normalized `D'`, the squared correlation `r^2`, and the repository's haplotype-block helper.[1][2][3][5] In the repository, `PopulationGeneticsAnalyzer.CalculateLD(...)` reports `DPrime` and `RSquared` for one variant pair, and `FindHaplotypeBlocks(...)` groups consecutive variants whose adjacent-pair LD exceeds a threshold.[6][7][8][9] The block caller is intentionally simplified relative to the full Gabriel et al. procedure.[3][6][8][9]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Two loci are in linkage disequilibrium when allele combinations occur together more or less often than expected under random association.[1][5] Haplotype blocks are genomic regions with high internal LD and comparatively limited historical recombination, so neighboring variants tend to be inherited together.[3][6]

### 2.2 Core Model

The classical disequilibrium coefficient is:[1][5]

$$
D = p_{AB} - p_A p_B
$$

where `p_AB` is the haplotype frequency of the `AB` combination and `p_A`, `p_B` are the marginal allele frequencies.[1][5]

Lewontin's normalized coefficient is:[1][5]

$$
D' = \frac{D}{D_{max}}
$$

In practice, many summaries report the magnitude `|D'|`; the repository's public `DPrime` field exposes that non-negative form, so observable outputs lie in `[0, 1]`.[6][9]

with

$$
D_{max} = \begin{cases}
\min(p_A p_B, q_A q_B), & D < 0 \\
\min(p_A q_B, q_A p_B), & D \ge 0
\end{cases}
$$

and `q_A = 1 - p_A`, `q_B = 1 - p_B`.[1][5]

The squared correlation measure is:[2][5]

$$
r^2 = \frac{D^2}{p_A q_A p_B q_B}
$$

with `0 <= r^2 <= 1`.[2][5] The current document also recorded the common practical convention that `r^2 >= 0.8` is often treated as strong LD.[5][6]

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `D = 0` under random association | `p_AB` then equals `p_A p_B` by definition.[1][5] |
| INV-02 | `0 <= r^2 <= 1` | `r^2` is a squared correlation measure.[2][5][6] |
| INV-03 | `0 <= |D'| <= 1` | `D'` normalizes `D` by its theoretical maximum magnitude.[1][5][6] |
| INV-04 | A haplotype block must contain at least two variants | One marker cannot define a multi-locus association block.[3][6] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[CalculateLD] variant1Id` | `string` | required | Identifier for the first variant | Preserved in the output record.[7][9] |
| `[CalculateLD] variant2Id` | `string` | required | Identifier for the second variant | Preserved in the output record.[7][9] |
| `[CalculateLD] genotypes` | `IEnumerable<(int Geno1, int Geno2)>` | required | Paired genotype values for the two loci | The repository documents `0 = homozygous major`, `1 = heterozygous`, `2 = homozygous minor`.[6][7][9] |
| `[CalculateLD] distance` | `int` | required | Physical distance between the two loci | Copied into the returned record as `Distance`.[7][9] |
| `[FindHaplotypeBlocks] variants` | `IEnumerable<(string VariantId, int Position, IReadOnlyList<int> Genotypes)>` | required | Variant metadata plus genotype vectors | Fewer than two variants produce no blocks.[6][7][9] |
| `[FindHaplotypeBlocks] ldThreshold` | `double` | `0.7` | Minimum adjacent-pair `r^2` required to keep extending a block | The repository default is `0.7`.[6][7][9] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `LinkageDisequilibrium.Variant1` | `string` | Identifier of the first input variant.[7][9] |
| `LinkageDisequilibrium.Variant2` | `string` | Identifier of the second input variant.[7][9] |
| `LinkageDisequilibrium.DPrime` | `double` | Normalized non-negative LD summary returned by the repository.[7][9] |
| `LinkageDisequilibrium.RSquared` | `double` | Squared correlation summary for the input genotype pair sequence.[2][7][9] |
| `LinkageDisequilibrium.Distance` | `double` | Distance copied from the input `distance` value.[7][9] |
| `HaplotypeBlock.Start` | `int` | Start position of one detected block.[7][9] |
| `HaplotypeBlock.End` | `int` | End position of one detected block.[7][9] |
| `HaplotypeBlock.Variants` | `IReadOnlyList<string>` | Variant IDs included in the block.[7][9] |
| `HaplotypeBlock.Haplotypes` | `IReadOnlyList<(string Haplotype, double Frequency)>` | Haplotype list placeholder returned by the record type; the current block finder leaves it empty.[9] |

### 3.3 Preconditions and Validation

`CalculateLD` returns `DPrime = 0` and `RSquared = 0` for an empty genotype-pair sequence while preserving the variant IDs and distance.[7][9] Monomorphic loci are protected from division-by-zero in the `r^2` computation by returning `0` when either genotype variance is zero.[7][9] `FindHaplotypeBlocks` sorts variants by position before analysis and returns no blocks when fewer than two variants are available.[6][7][9] The current implementation assumes the documented `0/1/2` genotype encoding and does not perform additional validation of genotype values.[9]

## 4. Algorithm

### 4.1 High-Level Steps

1. For `CalculateLD`, materialize the paired genotype sequence.
2. Estimate locus-wise allele frequencies from the genotype values.
3. Compute the LD summary statistics for that variant pair and return them with the input IDs and distance.
4. For `FindHaplotypeBlocks`, sort the input variants by genomic position.
5. Compute LD for each adjacent variant pair and extend the current block while the adjacent-pair `RSquared` value stays at or above `ldThreshold`.
6. Emit only blocks containing at least two variants.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository's public contract uses genotype values as diploid allele-count encodings (`0`, `1`, `2`) and a default block-extension threshold of `0.7`.[6][7][9] Block membership is therefore based on runs of adjacent pairs whose returned `RSquared` meets that threshold, not on a separate haplotype-frequency table.[7][9]

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateLD` | `O(n)` | `O(n)` | Materializes `n` paired genotypes into a list before computing means, covariance, and variances.[7][9] |
| `FindHaplotypeBlocks` | `O(v log v + v * g)` | `O(v + g)` | `v` = number of variants and `g` = genotype-vector length per adjacent comparison; the method sorts variants and then zips each adjacent genotype pair once.[9] |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)

- `PopulationGeneticsAnalyzer.CalculateLD(string, string, IEnumerable<(int Geno1, int Geno2)>, int)`: Computes `DPrime` and `RSquared` for one variant pair.
- `PopulationGeneticsAnalyzer.FindHaplotypeBlocks(IEnumerable<(string VariantId, int Position, IReadOnlyList<int> Genotypes)>, double)`: Detects simplified haplotype blocks from adjacent-pair LD.

### 5.2 Current Behavior

`CalculateLD` computes `RSquared` as the squared Pearson correlation of genotype values `0`, `1`, and `2`, and it estimates `D` from diploid genotype covariance using `D = Cov(X1, X2) / 2`.[6][9] `DPrime` is returned as a non-negative quantity because the implementation uses `abs(D) / DMax` and clamps the result to `1.0`.[7][9] `FindHaplotypeBlocks` sorts input variants by position, considers only adjacent pairs, and builds blocks from contiguous runs whose adjacent `RSquared` values satisfy the threshold.[7][9] The returned `HaplotypeBlock.Haplotypes` lists are always empty, and genotype vectors of unequal length are truncated to their shared prefix because the method pairs them with `Zip(...)`.[9]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Pairwise LD summaries based on `r^2` and normalized `D'`.[1][2][6][9]
- Haplotype blocks emitted only when at least two variants remain in a high-LD run.[3][6][7][9]

**Intentionally simplified:**

- `D` is inferred from diploid genotype covariance rather than from phased haplotype counts; **consequence:** the API works on unphased genotype vectors but does not expose haplotype frequencies.
- Block detection uses only adjacent-pair `RSquared >= ldThreshold`; **consequence:** non-adjacent LD evidence and the full Gabriel confidence-interval procedure are not part of the result.
- `HaplotypeBlock.Haplotypes` is left empty; **consequence:** callers receive block membership and boundaries, not haplotype composition.

**Not implemented:**

- Phased haplotype estimation and explicit haplotype-frequency reporting; **users should rely on:** no current alternative.
- Full Gabriel confidence-interval block calling; **users should rely on:** no current alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Diploid genotype covariance is used as the practical proxy for haplotype disequilibrium | Deviation | `DPrime` and `RSquared` are available from unphased genotype vectors, but explicit haplotype frequencies are not | accepted | Confirmed in the source and evidence documents.[6][9] |
| 2 | Haplotype blocks are defined from adjacent-pair thresholding only | Deviation | Long-range evidence and full Gabriel confidence intervals do not affect block boundaries | accepted | The method calls `CalculateLD` only for consecutive position-sorted variants.[7][9] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty genotype sequence | Returns `RSquared = 0` and `DPrime = 0` | The implementation special-cases empty input.[7][9] |
| Monomorphic locus | Returns `RSquared = 0` | The method suppresses division by zero when either genotype variance is zero.[7][9] |
| Perfect LD with identical genotype vectors | Returns `RSquared = 1` and `DPrime = 1` | Identical vectors have perfect correlation in the test suite and source rationale.[6][7] |
| Balanced independent design | Returns `RSquared = 0` and `DPrime = 0` | The tests use a covariance-zero construction to represent no LD.[6][7] |
| Single variant in block detection | Returns no blocks | A block requires at least two variants.[3][6][7] |

### 6.2 Limitations

The repository does not infer phased haplotypes, does not report haplotype frequencies, and does not implement the full Gabriel block-definition procedure.[3][6][9] Block detection only considers adjacent pairs, so non-adjacent LD relationships do not influence block boundaries.[7][9]

## 8. References

1. Lewontin, R. C. 1964. The interaction of selection and linkage. Genetics 49(1):49-67.
2. Hill, W. G., and A. Robertson. 1968. Linkage disequilibrium in finite populations. Theoretical and Applied Genetics 38(6):226-231.
3. Gabriel, S. B., et al. 2002. The Structure of Haplotype Blocks in the Human Genome. Science 296(5576):2225-2229.
4. Wikipedia contributors. Linkage disequilibrium. https://en.wikipedia.org/wiki/Linkage_disequilibrium
5. Wikipedia contributors. Haplotype block. https://en.wikipedia.org/wiki/Haplotype_block
6. [POP-LD-001.md](../../../docs/Evidence/POP-LD-001.md)
7. [POP-LD-001.md](../../../tests/TestSpecs/POP-LD-001.md)
8. [PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs)
9. [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)
