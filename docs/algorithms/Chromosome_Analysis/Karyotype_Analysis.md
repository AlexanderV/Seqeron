# Karyotype Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Chromosome Analysis |
| Test Unit ID | CHROM-KARYO-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Karyotype analysis summarizes chromosome composition by counting chromosomes, separating autosomes from sex chromosomes, and reporting whole-chromosome count abnormalities. In this repository, `AnalyzeKaryotype` operates on an explicit list of chromosome tuples, while `DetectPloidy` estimates ploidy from normalized read depth values. The terminology for monosomy, trisomy, tetrasomy, and pentasomy follows standard cytogenetic usage.[1][3] The implementation is intentionally compact: it is useful for chromosome-set summaries and depth-based ploidy heuristics, but it does not model full cytogenetic morphology, banding, or sex-chromosome abnormality logic.[1][4]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A karyotype is the complete chromosome complement of a cell, including chromosome number and the organization of the chromosome set.[1] The existing repository documentation also highlights standard ploidy terminology.[2]

| State | Meaning |
|-------|---------|
| Haploid (`n`) | One chromosome set |
| Diploid (`2n`) | Two chromosome sets |
| Polyploid (`>2n`) | More than two chromosome sets |

Aneuploidy refers to gain or loss of individual chromosomes rather than whole chromosome sets.[3] Examples already documented in the repository include Down syndrome (trisomy 21) and Turner syndrome (monosomy X).[3]

### 2.2 Core Model

`AnalyzeKaryotype` uses absolute chromosome-copy counts per chromosome group and labels any non-expected autosome count with standard cytogenetic terminology. `DetectPloidy` uses the ratio between the median normalized depth and an expected diploid depth to estimate the most likely ploidy level, with a confidence score that decreases as the observed ratio moves away from an integer copy-state boundary.[2][3]

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Autosome copies use a shared base name with optional numeric suffixes such as `chr1_1`, `chr1_2`. | Homologous chromosomes can be split into separate groups and abnormalities can be under- or over-counted. |
| ASM-02 | `normalizedDepths` are on a scale where diploid material is represented by `expectedDiploidDepth`. | Ploidy estimates shift because the median depth ratio no longer reflects chromosome-set copy number. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `TotalChromosomes == AutosomeCount + SexChromosomes.Count`. | `AnalyzeKaryotype` partitions the input into autosomes and sex chromosomes before computing the summary. |
| INV-02 | `TotalGenomeSize` equals the sum of input chromosome lengths, and `MeanChromosomeLength` is derived from that total when chromosomes are present. | The implementation computes both directly from the materialized input list. |
| INV-03 | `HasAneuploidy` is true exactly when the abnormality list is non-empty. | The flag is set only when a chromosome group count differs from `expectedPloidyLevel`. |
| INV-04 | `DetectPloidy` returns `PloidyLevel` in `[1, 8]` and `Confidence` in `[0, 1]` for non-empty input, and `(2, 0)` for empty input. | The implementation clamps ploidy to `[1, 8]`, clamps confidence at zero, and special-cases empty input. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[AnalyzeKaryotype] chromosomes` | `IEnumerable<(string Name, long Length, bool IsSexChromosome)>` | required | Chromosome descriptors used to build a karyotype summary. | Empty input returns a zeroed `Karyotype`. Autosomes are grouped by base name after removing an underscore followed by a numeric suffix. |
| `[AnalyzeKaryotype] expectedPloidyLevel` | `int` | `2` | Expected autosome copy count used to decide whether a group is abnormal. | Copied into the returned `Karyotype.PloidyLevel` for non-empty input. |
| `[DetectPloidy] normalizedDepths` | `IEnumerable<double>` | required | Normalized depth values used to estimate ploidy. | Empty input returns `(2, 0)`. |
| `[DetectPloidy] expectedDiploidDepth` | `double` | `1.0` | Baseline depth corresponding to diploid material. | Used as the denominator for the median-depth ratio. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Karyotype.TotalChromosomes` | `int` | Number of chromosome tuples provided to `AnalyzeKaryotype`. |
| `Karyotype.AutosomeCount` | `int` | Number of chromosomes marked as non-sex chromosomes. |
| `Karyotype.SexChromosomes` | `IReadOnlyList<string>` | Names of chromosomes flagged as sex chromosomes. |
| `Karyotype.TotalGenomeSize` | `long` | Sum of all chromosome lengths. |
| `Karyotype.MeanChromosomeLength` | `double` | Average chromosome length across the input set. |
| `Karyotype.PloidyLevel` | `int` | Expected ploidy used for the analysis, or `0` for empty input. |
| `Karyotype.HasAneuploidy` | `bool` | Whether any autosome group count differs from `expectedPloidyLevel`. |
| `Karyotype.Abnormalities` | `IReadOnlyList<string>` | Per-group abnormality labels such as `Trisomy chr21`. |
| `DetectPloidy.PloidyLevel` | `int` | Median-depth-based ploidy estimate. |
| `DetectPloidy.Confidence` | `double` | Distance-from-integer confidence score for that estimate. |

### 3.3 Preconditions and Validation

`AnalyzeKaryotype` materializes its input and returns an empty `Karyotype` when no chromosomes are provided. The helper that derives a chromosome base name only strips an underscore followed by a numeric suffix, so names outside that convention are treated literally. `DetectPloidy` materializes and sorts the input values, computes the median of the list, rounds `ratio * 2` to the nearest ploidy, clamps the result to `[1, 8]`, and returns `(2, 0)` when no depths are supplied.

## 4. Algorithm

### 4.1 High-Level Steps

1. Materialize the chromosome list or normalized-depth sequence.
2. For `AnalyzeKaryotype`, partition chromosomes into sex chromosomes and autosomes.
3. Group autosomes by base chromosome name and mark any group whose count differs from `expectedPloidyLevel` with an absolute cytogenetic term.
4. Compute total chromosome count, total genome size, and mean chromosome length, then return a `Karyotype` summary.
5. For `DetectPloidy`, sort the depth values and compute the true median.
6. Convert the median-depth ratio into a rounded ploidy estimate, clamp it to `[1, 8]`, compute confidence from the distance to the nearest integer copy state, and return the result.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository implementation uses the following operational rules for ploidy estimation:

```text
medianDepth = median(sortedDepths)
ratio = medianDepth / expectedDiploidDepth
ploidy = round(ratio * 2)
ploidy = clamp(ploidy, 1, 8)
fractionalPart = |ratio * 2 - ploidy|
confidence = max(0, 1 - fractionalPart * 2)
```

Autosome group counts are labeled with the following absolute terms: `0 -> Nullisomy`, `1 -> Monosomy`, `2 -> Disomy`, `3 -> Trisomy`, `4 -> Tetrasomy`, `5 -> Pentasomy`, and `Polysomy (N copies)` for higher counts.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `AnalyzeKaryotype` | `O(n)` | `O(n)` | `n` is the number of chromosome tuples. |
| `DetectPloidy` | `O(n log n)` | `O(n)` | Sorting the depth values dominates the median calculation. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ChromosomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs)

- `ChromosomeAnalyzer.AnalyzeKaryotype(...)`: builds a `Karyotype` summary from chromosome tuples.
- `ChromosomeAnalyzer.DetectPloidy(...)`: estimates ploidy from normalized depth values.

### 5.2 Current Behavior

`AnalyzeKaryotype` only evaluates autosome groups when building the abnormality list; sex chromosomes are preserved in `SexChromosomes` but are not classified as abnormal or normal. The method strips only underscore-plus-number suffixes when collapsing names to a base chromosome identifier. `DetectPloidy` returns `(2, 0)` for empty input, uses the exact median of the sorted depth list, and clamps extreme estimates into the range `[1, 8]`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Absolute cytogenetic terminology such as `Monosomy`, `Trisomy`, `Tetrasomy`, and `Pentasomy` is used for chromosome counts.[3]
- Diploid-like, haploid-like, and tetraploid-like depth ratios map to `2`, `1`, and `4` respectively in the ploidy estimator.[2][3]

**Intentionally simplified:**

- `AnalyzeKaryotype` treats chromosome abnormality detection as a count-comparison problem over autosome groups; **consequence:** morphology, centromere position, and banding-pattern information are outside the result contract.
- `DetectPloidy` uses only the median normalized depth and not a full ploidy model; **consequence:** local mixture structure and copy-number heterogeneity are collapsed into a single global estimate.

**Not implemented:**

- Sex-chromosome abnormality calling in `AnalyzeKaryotype`; **users should rely on:** `DetectAneuploidy` for depth-based chromosome-state summaries, with awareness of its own documented limitations.
- Banding-pattern or structural-karyotype analysis; **users should rely on:** no current alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `AnalyzeKaryotype` excludes sex chromosomes from abnormality grouping. | Deviation | Sex-chromosome abnormalities are not surfaced in `Abnormalities`. | accepted | Directly confirmed from the source code and covered by the separation invariant in [CHROM-KARYO-001.md](../../../tests/TestSpecs/CHROM-KARYO-001.md). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty chromosome list | Returns a zeroed `Karyotype` with no abnormalities. | Source code special-cases empty input. |
| Empty depth list | Returns `(2, 0)`. | The implementation treats diploid as the default with zero confidence when no data are available. |
| `ratio = 0.5` | Returns ploidy `1` with high confidence. | This is the documented haploid boundary case. |
| `ratio = 1.0` | Returns ploidy `2` with high confidence. | This is the documented diploid boundary case. |
| `ratio = 2.0` | Returns ploidy `4` with high confidence. | This is the documented tetraploid boundary case. |
| Extreme depth ratios | Ploidy is clamped into `[1, 8]`. | The implementation bounds the estimate to a fixed range. |

### 6.2 Limitations

The repository implementation provides count-based karyotype summaries and a depth-based ploidy heuristic, not a full cytogenetic analysis. `AnalyzeKaryotype` depends on a chromosome-naming convention to group homologs and does not classify sex-chromosome abnormalities. `DetectPloidy` reduces the entire depth profile to a single median-derived estimate, so it does not distinguish mixed ploidy states or regional copy-number structure.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ChromosomeAnalyzer_Karyotype_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_Karyotype_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [CHROM-KARYO-001.md](../../../tests/TestSpecs/CHROM-KARYO-001.md)
- Related algorithms: [Aneuploidy_Detection.md](Aneuploidy_Detection.md)
- Related algorithms: [Centromere_Analysis.md](Centromere_Analysis.md)
- Related algorithms: [Telomere_Analysis.md](Telomere_Analysis.md)

## 8. References

1. Wikipedia contributors. 2026. Karyotype. Wikipedia. https://en.wikipedia.org/wiki/Karyotype
2. Wikipedia contributors. 2026. Ploidy. Wikipedia. https://en.wikipedia.org/wiki/Ploidy
3. Wikipedia contributors. 2026. Aneuploidy. Wikipedia. https://en.wikipedia.org/wiki/Aneuploidy
4. Tjio JH, Levan A. 1956. The chromosome number of man. Hereditas. N/A
