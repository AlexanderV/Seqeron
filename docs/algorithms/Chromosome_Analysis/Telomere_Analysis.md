# Telomere Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Chromosome Analysis |
| Test Unit ID | CHROM-TELO-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Telomere analysis detects telomeric repeat tracts at chromosome ends and provides a separate helper for converting qPCR T/S ratios into approximate telomere length. In this repository, `AnalyzeTelomeres` is a port of Heng Li's `seqtk telo` maximal-scoring-segment scan:[6] it scans the 5' and 3' ends of a supplied sequence for motif-rotation hits, reports end-specific tract lengths and repeat purities, and flags critically short telomeres based on configurable thresholds. `EstimateTelomereLengthFromTSRatio` applies the proportional T/S-ratio relationship described in the cited qPCR literature.[3] The sequence-based detection is heuristic and end-focused, so it is appropriate for approximate repeat-end assessment rather than complete telomere biology inference.[1][2]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A telomere is a repetitive DNA structure at a chromosome end that protects chromosome termini from degradation and fusion.[1] The existing repository documentation records the following background facts.[1][2]

| Feature | Value |
|---------|-------|
| Vertebrate canonical repeat | `TTAGGG` |
| Reverse complement at the 5' end | `CCCTAA` |
| Repeat unit length | 6 bp |
| Human telomere length at birth | 5,000–15,000 bp |

The documented chromosome-end orientation is:

| End | Repeat Sequence | Direction |
|-----|-----------------|-----------|
| 5' end | `CCCTAA` | Toward the chromosome interior |
| 3' end | `TTAGGG` | Toward the chromosome terminus |

The original document also records several organism-specific telomeric repeat patterns.[1][4]

| Organism | Repeat | Notes |
|----------|--------|-------|
| Vertebrates | `TTAGGG` | Conserved across vertebrates |
| Arabidopsis | `TTTAGGG` | 7-bp repeat |
| Tetrahymena | `TTGGGG` | Discovery organism for a canonical telomere repeat |
| S. cerevisiae | Variable | Irregular repeat pattern |

The existing repository documentation further states a normal human range of 5,000–15,000 bp, a critical threshold around 3,000 bp, and an association between short telomeres and aging or disease risk.[5]

### 2.2 Core Model

The sequence-based model is that telomeric repeats occur at chromosome ends. Following `seqtk telo`,[6] every k-mer (k = motif length) equal to any rotation of the expected end motif is a hit; scanning inward from each terminus, each scored position adds `+1` for a hit and `-penalty` otherwise, the tract ends at the position of maximal cumulative score, and the scan stops once the score falls more than `maxDrop` below the maximum (X-drop). Using rotations makes detection independent of the repeat phase at the terminus. For qPCR data, the repository uses the proportional T/S-ratio model described by Cawthon (2002):

```text
estimatedLength = referenceLength * (tsRatio / referenceRatio)
```

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Canonical telomeric repeats are concentrated within `searchLength` bases of the sequence ends. | End scanning can underestimate or miss telomeres that extend deeper into the sequence or are absent from the provided ends. |
| ASM-02 | The supplied repeat unit matches the organism being analyzed. | End-specific length and purity estimates become unreliable because the comparison uses the wrong motif. |
| ASM-03 | The measured T/S ratio is proportional to average telomere length for the assay context. | `EstimateTelomereLengthFromTSRatio` returns a scaled value that may not correspond to actual telomere length. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 <= TelomereLength5Prime, TelomereLength3Prime <= min(searchLength, n)`. | Lengths are `maxPos + 1` (5') / `n - maxPos` (3') inside the scanned window, `0` when the maximal score is `<= 0`. |
| INV-02 | `0 <= RepeatPurity <= 1`; for a non-empty tract `RepeatPurity > 1/(1+penalty)`. | Purity = hits / scored positions in the tract; the tract's score `hits - penalty*misses` is `> 0`. |
| INV-03 | `Has5PrimeTelomere` and `Has3PrimeTelomere` are derived from the measured lengths and `minTelomereLength`. | The public method compares each measured length to the threshold after scanning. |
| INV-04 | For non-negative `tsRatio`, the T/S helper returns a non-negative result and follows `referenceLength * tsRatio / referenceRatio`. | The implementation is a direct proportional calculation. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[AnalyzeTelomeres] chromosomeName` | `string` | required | Chromosome identifier copied into the result. | Preserved verbatim in `TelomereResult.Chromosome`. |
| `[AnalyzeTelomeres] sequence` | `string` | required | DNA sequence whose ends are scanned for telomeric repeats. | Empty or `null` input returns no telomeres and `IsCriticallyShort = true`. |
| `[AnalyzeTelomeres] telomereRepeat` | `string` | `"TTAGGG"` | Repeat unit used for the 3' scan. | The 5' scan uses its reverse complement; must be non-empty A/C/G/T (case-insensitive), else `ArgumentException`. |
| `[AnalyzeTelomeres] searchLength` | `int` | `10000` | Maximum distance from each chromosome end to inspect. | Effective scan length is `min(searchLength, sequence.Length)` on each end. |
| `[AnalyzeTelomeres] minTelomereLength` | `int` | `500` | Minimum measured repeat length required for `Has*Telomere` to be true. | Applied independently to the 5' and 3' ends. |
| `[AnalyzeTelomeres] criticalLength` | `int` | `3000` | Threshold used for the `IsCriticallyShort` flag when a telomere is detected. | Applied only after end-specific presence/absence is determined. |
| `[AnalyzeTelomeres] penalty` | `int` | `1` | Score for a non-hit position is `-penalty` (seqtk `-p`). | Sign ignored, as in seqtk. |
| `[AnalyzeTelomeres] maxDrop` | `int` | `2000` | X-drop: stop scanning when the score falls more than this below the maximum (seqtk `-d`). | |
| `[EstimateTelomereLengthFromTSRatio] tsRatio` | `double` | required | Observed T/S ratio. | Interpreted as proportional to average telomere length. |
| `[EstimateTelomereLengthFromTSRatio] referenceRatio` | `double` | `1.0` | Reference-sample T/S ratio. | Used as the denominator in the proportional formula. |
| `[EstimateTelomereLengthFromTSRatio] referenceLength` | `double` | `7000` | Reference-sample telomere length in base pairs. | Used as the scale factor in the proportional formula. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `TelomereResult.Chromosome` | `string` | Chromosome name passed to the method. |
| `TelomereResult.Has5PrimeTelomere` | `bool` | Whether the measured 5' repeat tract reaches `minTelomereLength`. |
| `TelomereResult.TelomereLength5Prime` | `int` | Measured 5' tract length in bases. |
| `TelomereResult.Has3PrimeTelomere` | `bool` | Whether the measured 3' repeat tract reaches `minTelomereLength`. |
| `TelomereResult.TelomereLength3Prime` | `int` | Measured 3' tract length in bases. |
| `TelomereResult.RepeatPurity5Prime` | `double` | Fraction of matching repeat bases in the counted 5' tract. |
| `TelomereResult.RepeatPurity3Prime` | `double` | Fraction of matching repeat bases in the counted 3' tract. |
| `TelomereResult.IsCriticallyShort` | `bool` | Flag indicating critically short detected telomeres, with a special-case `true` result for empty input. |
| `EstimatedLength` | `double` | Base-pair estimate returned by `EstimateTelomereLengthFromTSRatio`. |

### 3.3 Preconditions and Validation

`AnalyzeTelomeres` uppercases the input sequence and repeat motif, computes the reverse complement of the repeat for the 5' end, and scans only the configured end windows. Any base other than A/C/G/T resets the rolling k-mer (no hit for the next k−1 positions) without stopping the scan. When the sequence is shorter than the repeat length, both end lengths remain zero. `EstimateTelomereLengthFromTSRatio` performs a direct proportional calculation and does not impose additional validation beyond the numeric inputs supplied by the caller.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return a no-telomere result with `IsCriticallyShort = true` for empty or `null` input.
2. Uppercase the sequence and repeat unit and compute the reverse complement of the repeat.
3. Build the hit sets: all rotations of the reverse-complement motif (5') and of the motif (3').
4. 5' end: for i = 0.. within the window, the k-mer ending at i is a hit if it is in the 5' set; from i ≥ k add `+1`/`-penalty`; track the maximum score and its position; stop when `max − score > maxDrop`. Length = `maxPos + 1` if the maximum is > 0.
5. 3' end: the same from i = n−1 downward (k-mer starting at i, scoring from n − i ≥ k), never entering an accepted 5' tract; length = `n − maxPos`.
6. Purity = hits / scored positions within each tract; `Has*Telomere` = length ≥ `minTelomereLength` (and > 0).
7. Mark the result as critically short when a detected telomere is shorter than `criticalLength`.
8. For qPCR data, compute `referenceLength * tsRatio / referenceRatio`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The documented sequence-orientation rules are:

| End | Expected Repeat |
|-----|-----------------|
| 5' end | Reverse complement of `telomereRepeat` |
| 3' end | `telomereRepeat` |

The qPCR helper uses the proportional T/S-ratio formula from Section 2.2.[3]

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `AnalyzeTelomeres` | `O(n)` | `O(1)` | `n = min(searchLength, sequence.Length)` for the scanned end windows in the original documentation. |
| `EstimateTelomereLengthFromTSRatio` | `O(1)` | `O(1)` | Direct arithmetic. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ChromosomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs)

- `ChromosomeAnalyzer.AnalyzeTelomeres(...)`: scans both ends for repeat tracts and returns `TelomereResult`.
- `ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(...)`: converts a T/S ratio to a base-pair estimate.

### 5.2 Current Behavior

The implementation is a line-by-line port of `stk_telo` (seqtk 1.5-r133), including its scoring offsets (5' positions scored from i ≥ k, so a single 5' unit yields no tract, while one 3' unit yields 6). Lengths are base-resolution and may include a terminal partial unit. Seqeron extensions: the `searchLength` window, presence gated by tract length (`minTelomereLength`) instead of seqtk's min score (default 300), and the purity statistic. For non-empty input, `IsCriticallyShort` becomes true only when a detected telomere is present and its measured length is below `criticalLength`; a non-empty sequence with no detected telomere yields `IsCriticallyShort = false`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Vertebrate-style `TTAGGG`/`CCCTAA` orientation at chromosome ends is modeled explicitly.[1][2]
- The T/S-ratio helper uses the proportional relationship `referenceLength * (tsRatio / referenceRatio)`.[3]

**Intentionally simplified:**

- End detection uses exact motif-rotation k-mer hits (seqtk telo); variant repeats (e.g. TCAGGG, TGAGGG) count as misses; **consequence:** tracts rich in variant repeats score lower or may be truncated.
- Only the first and last `searchLength` bases are examined; **consequence:** end-distal or truncated sequence context can hide telomeric sequence outside the scanned windows.
- The default repeat and threshold values are configurable but not species-specific; **consequence:** non-vertebrate use cases require callers to supply an appropriate repeat motif and interpretability remains heuristic.

**Not implemented:**

- Detection of interstitial telomeric repeats away from chromosome ends; **users should rely on:** no current alternative.
- Assay calibration, population normalization, or clinical interpretation beyond the proportional T/S estimate; **users should rely on:** no current alternative.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns no telomeres and `IsCriticallyShort = true`. | The public method special-cases empty input. |
| Sequence shorter than the repeat unit | Returns zero-length telomeres. | The internal repeat matcher exits immediately when the region is shorter than the repeat. |
| No telomeric repeats | Returns zero lengths and `Has*Telomere = false`. | No scanned window meets the similarity threshold. |
| Divergent repeats | Sporadic non-motif units lower purity (e.g. every 10th unit TTAGGA over 200 units: length 1194, purity 1075/1189); a tract made only of a non-motif hexamer (TTAGGA×200) is not detected. | Values from `seqtk telo`. |
| Terminal partial repeat | Detected with base-resolution length (A×1000+(TTAGGG)×200+`TTAG` → 1204). | Motif rotations are phase-independent. |
| Custom repeat motif | Uses the supplied motif and its reverse complement. | The repeat string is a public parameter. |

### 6.2 Limitations

The repository implementation is an end-focused repeat scanner. It does not model telomere-associated proteins, interstitial telomeric repeats, assay-specific calibration, or organism-specific defaults beyond the caller-supplied repeat motif. The helper for T/S ratios is a direct proportional conversion and does not by itself validate assay quality or biological interpretation.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ChromosomeAnalyzer_Telomere_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_Telomere_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [CHROM-TELO-001.md](../../../tests/TestSpecs/CHROM-TELO-001.md)
- Related algorithms: [Centromere_Analysis.md](Centromere_Analysis.md)

## 8. References

1. Wikipedia contributors. 2026. Telomere. Wikipedia. https://en.wikipedia.org/wiki/Telomere
2. Meyne J, Ratliff RL, Moyzis RK. 1989. Conservation of the human telomere sequence (TTAGGG)n among vertebrates. Proceedings of the National Academy of Sciences. N/A
3. Cawthon RM. 2002. Telomere measurement by quantitative PCR. Nucleic Acids Research. doi:10.1093/nar/30.10.e47
4. Blackburn EH, Gall JG. 1978. A tandemly repeated sequence at the termini of the extrachromosomal ribosomal RNA genes in Tetrahymena. Journal of Molecular Biology. N/A
5. Rossiello et al. 2022. N/A. Nature Cell Biology. N/A
6. Li H. seqtk 1.5-r133, `stk_telo` (seqtk.c). https://github.com/lh3/seqtk (source: https://raw.githubusercontent.com/lh3/seqtk/master/seqtk.c)
