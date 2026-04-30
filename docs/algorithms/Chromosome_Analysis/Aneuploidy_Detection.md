# Aneuploidy Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Chromosome Analysis |
| Test Unit ID | CHROM-ANEU-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Aneuploidy detection estimates abnormal chromosome copy number from sequencing depth and classifies whole-chromosome gains or losses. In this repository, `DetectAneuploidy` aggregates depth into genomic bins, converts the depth ratio relative to a diploid baseline into a log2 ratio and integer copy number, and emits per-bin `CopyNumberState` values. `IdentifyWholeChromosomeAneuploidy` then summarizes those states per chromosome and reports the dominant non-disomic class when it meets a minimum-fraction threshold. The implementation is a heuristic, read-depth-based approximation rather than a full clinical copy-number pipeline. Standard cytogenetic terms such as nullisomy, monosomy, trisomy, tetrasomy, and pentasomy follow the definitions summarized in the repository sources.[1][3]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Aneuploidy is the presence of an abnormal number of chromosomes in a cell; in humans the normal diploid complement is 46 chromosomes, and aneuploidy refers to gains or losses of individual chromosomes rather than whole chromosome sets.[1][3]

Standard copy-number terminology used by this document is summarized below.[1]

| Term | Definition | Copy Number |
|------|------------|-------------|
| Nullisomy | Complete absence of a chromosome pair | 0 |
| Monosomy | Single copy instead of a pair | 1 |
| Disomy | Normal two-copy diploid state | 2 |
| Trisomy | Three copies | 3 |
| Tetrasomy | Four copies | 4 |
| Pentasomy | Five copies | 5 |

Examples of human aneuploidies compatible with live birth that appear in the existing repository documentation are listed below.[1]

| Condition | Karyotype | Description |
|-----------|-----------|-------------|
| Down syndrome | 47,XX,+21 or 47,XY,+21 | Trisomy 21 |
| Edwards syndrome | 47,XX,+18 or 47,XY,+18 | Trisomy 18 |
| Patau syndrome | 47,XX,+13 or 47,XY,+13 | Trisomy 13 |
| Turner syndrome | 45,X | Monosomy X |
| Klinefelter syndrome | 47,XXY | Extra X chromosome |

### 2.2 Core Model

The documented biological model is that sequencing depth is proportional to chromosome copy number when a diploid baseline is available. Under that proportional model, a one-copy region is expected to have about half the read depth of a two-copy region, a three-copy region about 1.5 times the depth, and a four-copy region about twice the depth.[2]

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | `medianDepth` is a valid diploid baseline for the sample being analyzed. | Copy-number estimates are systematically shifted because all ratios are computed against the wrong baseline. |
| ASM-02 | Read depth is approximately proportional to copy number over the bins being compared. | Integer copy-number states inferred from depth ratios can be misleading even when the arithmetic is internally consistent. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Reported bin copy numbers are integers in the range `[0, 10]`. | The implementation rounds the inferred copy number and clamps the result to that interval before emitting `CopyNumberState`. |
| INV-02 | Reported confidence scores are in the range `[0, 1]`. | Confidence is computed as `1 - min(1, |expected - observed|)`, which cannot exceed that range. |
| INV-03 | Whole-chromosome classifications are emitted only for dominant non-disomic states whose fraction is at least `minFraction`. | `IdentifyWholeChromosomeAneuploidy` groups states by chromosome and copy number and only reports the dominant state when it meets the threshold and `CopyNumber != 2`. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[DetectAneuploidy] depthData` | `IEnumerable<(string Chromosome, int Position, double Depth)>` | required | Depth observations used to infer copy number. | Empty input produces no output. Positions are grouped by `Position / binSize`, and chromosome names are processed independently. |
| `[DetectAneuploidy] medianDepth` | `double` | required | Expected diploid depth baseline. | Must be `> 0`; otherwise the method yields no states. |
| `[DetectAneuploidy] binSize` | `int` | `1000000` | Bin width used to aggregate depth values. | Start and end coordinates are emitted as `binIndex * binSize` and `(binIndex + 1) * binSize - 1`. |
| `[IdentifyWholeChromosomeAneuploidy] copyNumberStates` | `IEnumerable<CopyNumberState>` | required | Bin-level copy-number states to summarize per chromosome. | Empty input produces no output. |
| `[IdentifyWholeChromosomeAneuploidy] minFraction` | `double` | `0.8` | Minimum dominant-state fraction required for whole-chromosome classification. | Only dominant non-disomic states at or above this threshold are reported. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CopyNumberState.Chromosome` | `string` | Chromosome name copied from the grouped input. |
| `CopyNumberState.Start` | `int` | Inclusive bin start coordinate computed from the integer bin index and `binSize`. |
| `CopyNumberState.End` | `int` | Inclusive bin end coordinate computed as `(binIndex + 1) * binSize - 1`. |
| `CopyNumberState.CopyNumber` | `int` | Rounded, clamped copy-number estimate for the bin. |
| `CopyNumberState.LogRatio` | `double` | `log2(meanDepth / medianDepth)` for the bin. |
| `CopyNumberState.Confidence` | `double` | Heuristic agreement score between the observed depth ratio and the rounded copy-number state. |
| `WholeChromosome.Chromosome` | `string` | Chromosome whose dominant bin-level state satisfied `minFraction`. |
| `WholeChromosome.CopyNumber` | `int` | Dominant copy number for that chromosome. |
| `WholeChromosome.Type` | `string` | Named class (`Nullisomy`, `Monosomy`, `Trisomy`, `Tetrasomy`, `Pentasomy`) or the fallback string `Copy number = N`. |

### 3.3 Preconditions and Validation

`DetectAneuploidy` materializes the input depth sequence and returns no results when the input is empty or when `medianDepth <= 0`. Depth values are averaged inside integer bins defined by `Position / binSize`; the method does not perform sex-chromosome normalization, ploidy normalization, or additional validation beyond those checks. `IdentifyWholeChromosomeAneuploidy` materializes its input, groups states by chromosome, and returns no output for chromosomes whose dominant copy number is disomic or whose dominant fraction is below `minFraction`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Materialize `depthData` and stop immediately if the sequence is empty or `medianDepth <= 0`.
2. Group observations by chromosome, then group each chromosome into bins using integer division by `binSize`.
3. For each bin, average the depth values, compute a log2 ratio against `medianDepth`, convert that ratio to an integer copy number, clamp the result to `[0, 10]`, and compute a confidence score.
4. Emit one `CopyNumberState` per chromosome/bin pair in ascending bin order within each chromosome.
5. Group `CopyNumberState` values by chromosome and then by `CopyNumber`.
6. For each chromosome, emit a whole-chromosome classification only when the dominant non-disomic copy number occupies at least `minFraction` of the bins.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository implementation uses the following operational rules for each binned depth measurement:

```text
logRatio = log2(meanDepth / medianDepth)
copyNumber = round((2 ^ logRatio) * 2)
expected = copyNumber / 2.0
observed = 2 ^ logRatio
confidence = 1.0 - min(1.0, |expected - observed|)
```

The depth-ratio interpretation documented in the original file is summarized below.

| Observed/Expected Ratio | Log2 Ratio | Copy Number |
|------------------------|------------|-------------|
| 0.0 | `-∞` | 0 |
| 0.5 | `-1.0` | 1 |
| 1.0 | `0.0` | 2 |
| 1.5 | `+0.58` | 3 |
| 2.0 | `+1.0` | 4 |
| 2.5 | `+1.32` | 5 |

Whole-chromosome labels are assigned with the following mapping from the dominant bin-level copy number: `0 -> Nullisomy`, `1 -> Monosomy`, `3 -> Trisomy`, `4 -> Tetrasomy`, `5 -> Pentasomy`, and `N -> Copy number = N` for other non-disomic values.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DetectAneuploidy` | `O(n)` | `O(n)` | `n` is the number of input depth observations. Existing repository documentation characterizes the method as linear in the input size. |
| `IdentifyWholeChromosomeAneuploidy` | `O(n)` | `O(n)` | `n` is the number of input `CopyNumberState` values. Existing repository documentation characterizes the method as linear in the input size. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ChromosomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs)

- `ChromosomeAnalyzer.DetectAneuploidy(...)`: bins depth values by chromosome and position, computes `CopyNumberState` outputs, and yields them in chromosome/bin order.
- `ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(...)`: summarizes per-bin states per chromosome and emits dominant non-disomic whole-chromosome calls.

### 5.2 Current Behavior

The implementation groups observations by chromosome before binning positions with integer division by `binSize`. Within each bin it uses the arithmetic mean depth, computes `logRatio` with `Math.Log2(meanDepth / medianDepth)`, rounds the implied copy number, and clamps it to `[0, 10]`. Output coordinates are the full bin bounds rather than the min/max positions seen in the bin. Whole-chromosome classification ignores dominant copy number `2`, uses `0.8` as the default dominance threshold, and falls back to the string `Copy number = N` for non-disomic copy numbers other than `0`, `1`, `3`, `4`, and `5`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Standard cytogenetic terminology for `Nullisomy`, `Monosomy`, `Disomy`, `Trisomy`, `Tetrasomy`, and `Pentasomy` is used for the corresponding copy-number states.[1]
- The documented ratio examples `0.5x -> 1`, `1.0x -> 2`, `1.5x -> 3`, and `2.0x -> 4` are realized by the repository’s depth-to-copy-number conversion.[2]

**Intentionally simplified:**

- Copy number is inferred from mean read depth in fixed bins against a single median baseline; **consequence:** the method behaves as a coarse read-depth heuristic rather than a full segmentation or bias-corrected CNV caller.
- Whole-chromosome classification is based on a dominant-fraction threshold (`minFraction`) over per-bin states; **consequence:** mixed states below the threshold are suppressed instead of being reported as mosaic or partial events.
- Confidence is a heuristic agreement score between the rounded copy-number state and the observed depth ratio; **consequence:** the value is bounded and useful for ordering or filtering, but it is not a calibrated probability.

**Not implemented:**

- Sex-chromosome normalization or special-casing; **users should rely on:** no current alternative.
- Partial monosomy/trisomy labels beyond the underlying per-bin `CopyNumberState` output; **users should rely on:** inspection of the per-bin states, with no current alternative for named partial-event calls.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Sex chromosomes are treated the same as autosomes. | Deviation | A normal male X or Y copy count can look monosomic under the current heuristic. | accepted | Documented in the original file and in [CHROM-ANEU-001.md](../../../tests/TestSpecs/CHROM-ANEU-001.md). |
| 2 | Regional copy-number changes are emitted per bin but not classified as partial chromosome events. | Deviation | Users must interpret bin-level states directly for partial aneuploidy. | accepted | The original file documents partial aneuploidy as detected at bin level rather than named at chromosome level. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input depth data | `DetectAneuploidy` yields no states. | The implementation returns immediately when the materialized input list is empty. |
| `medianDepth = 0` or `medianDepth < 0` | `DetectAneuploidy` yields no states. | The implementation rejects non-positive baselines to avoid invalid depth-ratio calculations. |
| Very high depth ratios | Copy number is capped at `10`. | The implementation clamps the rounded copy-number estimate to `[0, 10]`. |
| Very low depth ratios | Copy number is floored at `0`. | The implementation clamps the rounded copy-number estimate to `[0, 10]`. |
| Mixed copy numbers on one chromosome | Whole-chromosome output depends on whether the dominant state reaches `minFraction`. | `IdentifyWholeChromosomeAneuploidy` reports only dominant non-disomic states at or above the threshold. |
| Single depth observation in one bin | A single `CopyNumberState` can still be produced. | The implementation averages whatever values are present in the bin and does not require multiple measurements. |

### 6.2 Limitations

This repository implementation is limited to a fixed-bin, read-depth heuristic. It does not perform sex-chromosome normalization, purity or ploidy correction, GC-bias correction, or segmentation into partial chromosome events. Regional copy-number changes are visible in the emitted per-bin states, but whole-chromosome classification remains a dominant-state summary.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ChromosomeAnalyzer_Aneuploidy_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Aneuploidy_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`
- Test specification: [CHROM-ANEU-001.md](../../../tests/TestSpecs/CHROM-ANEU-001.md)
- Related algorithms: [Karyotype_Analysis.md](Karyotype_Analysis.md)

## 8. References

1. Wikipedia contributors. 2026. Aneuploidy. Wikipedia. https://en.wikipedia.org/wiki/Aneuploidy
2. Wikipedia contributors. 2026. Copy number variation. Wikipedia. https://en.wikipedia.org/wiki/Copy_number_variation
3. Griffiths AJF, Miller JH, Suzuki DT, Lewontin RC, Gelbart WM. 2000. An Introduction to Genetic Analysis. 7th ed. W. H. Freeman. N/A
