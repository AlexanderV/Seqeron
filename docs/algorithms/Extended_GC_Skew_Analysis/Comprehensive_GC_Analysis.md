# Comprehensive GC Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Composition / Extended GC-Skew Analysis |
| Test Unit ID | SEQ-GC-ANALYSIS-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Comprehensive GC analysis summarizes the nucleotide composition of a DNA sequence in a single pass: the overall GC content (percentage of G+C bases), the overall GC skew and AT skew (strand-asymmetry indicators), the GC-skew and GC-content profiles over sliding windows, and the compositional variability (population variance) of those windowed profiles. It is a specification-driven aggregation of exact, closed-form composition statistics — there is no estimation or heuristic; every reported number is determined exactly by the base counts. It is used for at-a-glance characterization of a sequence's compositional structure and as a precursor to replication-origin analysis [1][4].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In double-stranded DNA the two strands replicate asymmetrically; the leading and lagging strands accumulate different mutational biases, producing measurable over- or under-representation of G relative to C (and A relative to T) along a chromosome [1][4]. The aggregate strength of these biases, together with the bulk G+C fraction, characterizes a genome's compositional landscape [3].

### 2.2 Core Model

For a sequence with base counts G, C, A, T:

- **GC content (percentage):** `GC% = (G + C) / (A + T + G + C) × 100` [3].
- **GC skew:** `GC skew = (G − C) / (G + C)`, range [−1, +1]; defined as 0 when G + C = 0 [1][2][5].
- **AT skew:** `AT skew = (A − T) / (A + T)`, range [−1, +1]; defined as 0 when A + T = 0 [6].
- **Windowed profiles:** GC skew and GC content evaluated on each sliding window of length `w` advancing by `step` ("multiple windows along the sequence") [5].
- **Compositional variability:** the **population** variance of the per-window values, `σ² = Σ(xᵢ − μ)² / N`, where N is the number of windows [7].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | OverallGcSkew ∈ [−1, +1] | numerator |G−C| ≤ denominator G+C [2] |
| INV-02 | OverallAtSkew ∈ [−1, +1] | |A−T| ≤ A+T [6] |
| INV-03 | OverallGcContent ∈ [0, 100] | G+C ≤ total bases, then ×100 [3] |
| INV-04 | GcContentVariance ≥ 0 and GcSkewVariance ≥ 0 | a sum of squares divided by N ≥ 0 [7] |
| INV-05 | Number of windows = ⌊(n−w)/step⌋ + 1 when n ≥ w, else 0 | only windows fully inside the sequence are emitted [5] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `DnaSequence` or `string` | required | DNA sequence; case-insensitive | only A/C/G/T counted, other symbols ignored |
| windowSize | `int` | 1000 | sliding-window length for profiles | ≥ 1 (validated by the windowed cores) |
| stepSize | `int` | 100 | step between window starts | ≥ 1 |

### 3.2 Output / Return Value

`GcAnalysisResult` (record):

| Field | Type | Description |
|-------|------|-------------|
| OverallGcContent | double | GC% over the whole sequence (0–100) |
| OverallGcSkew | double | (G−C)/(G+C) over the whole sequence |
| OverallAtSkew | double | (A−T)/(A+T) over the whole sequence |
| GcContentVariance | double | population variance of windowed GC% |
| GcSkewVariance | double | population variance of windowed GC skew |
| WindowedGcSkew | IReadOnlyList&lt;GcSkewPoint&gt; | per-window GC skew with positions |
| WindowedGcContent | IReadOnlyList&lt;GcContentPoint&gt; | per-window GC% with positions |
| SequenceLength | int | length of the analyzed sequence |

### 3.3 Preconditions and Validation

A null `DnaSequence` throws `ArgumentNullException`. A null/empty string returns a zero result with empty windowed lists and `SequenceLength = 0`. Counting is case-insensitive (uppercased internally); only A/C/G/T affect the metrics — ambiguous/other symbols are ignored in numerators and denominators (matching Biopython `GC_skew`, which ignores ambiguous bases) [5]. Indexing of window positions is 0-based; `WindowStart`/`WindowEnd` are inclusive and `Position` is the window midpoint `start + windowSize/2`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input to uppercase (string overload) or read `DnaSequence.Sequence`.
2. Compute windowed GC-skew and GC-content profiles over all full windows.
3. Compute overall GC content, GC skew, AT skew across the whole sequence.
4. Compute the population variance of the windowed GC-content and GC-skew values (0 if there are no windows).
5. Assemble the `GcAnalysisResult`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| AnalyzeGcContent | O(n + W·w) | O(W) | n = sequence length, W = number of windows, w = windowSize; each window is recounted independently. Overall scalar metrics are O(n). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GcSkewCalculator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs)

- `GcSkewCalculator.AnalyzeGcContent(DnaSequence, windowSize, stepSize)`: canonical entry point; validates non-null and delegates to the core.
- `GcSkewCalculator.AnalyzeGcContent(string, windowSize, stepSize)`: string overload (API parity with sibling methods); zero result for null/empty.

### 5.2 Current Behavior

Each window's counts are recomputed independently (no incremental sliding accumulator); this keeps the code identical to the per-window cores already used for `CalculateWindowedGcSkew`. The variances are population variances (÷N), consistent with treating the emitted windows as the complete population. The unit performs counting/aggregation only — no substring search or pattern matching — so the repository suffix tree is not applicable and was not used.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- GC content `(G+C)/(A+T+G+C)×100` [3].
- GC skew `(G−C)/(G+C)` with G+C=0 → 0 [1][2][5].
- AT skew `(A−T)/(A+T)` with A+T=0 → 0 [6].
- Population variance `Σ(xᵢ−μ)²/N` of windowed values [7].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Cumulative GC-skew diagram and explicit origin/terminus calls are out of scope here; **users should rely on:** `GcSkewCalculator.PredictReplicationOrigin` (SEQ-REPLICATION-001) and `CalculateCumulativeGcSkew`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | GC content reported as percentage, not Biopython [0,1] fraction | Assumption | value differs from `gc_fraction` by ×100 | accepted | repository/Brock convention [3]; locked by tests |
| 2 | "Variability" = population variance (÷N), not sample (÷N−1) | Assumption | smaller than Bessel-corrected variance | accepted | windows are the full population [7] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null `DnaSequence` | `ArgumentNullException` | contract parity with sibling methods |
| null/empty string | zero result, empty windows, length 0 | string-overload contract |
| sequence shorter than window | empty windowed lists, variances 0, scalars still computed | only full windows are emitted [5] |
| no G/C bases | OverallGcSkew = 0, GcContent = 0 | zero-division → 0 [5]; numerator 0 [3] |
| pure-G / pure-C sequence | OverallGcSkew = +1 / −1, GcContent = 100 | skew bounds [2] |

### 6.2 Limitations

Windows are recomputed per step (no incremental optimization); for very large windows this is O(W·w). Only A/C/G/T are counted; degenerate IUPAC codes do not contribute. The aggregation does not itself locate replication origins — use the dedicated origin predictor.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** `GGGCCAT` — G=3, C=2, A=1, T=1, length 7.
GC content = (3+2)/7×100 = 71.42857142857143; GC skew = (3−2)/(3+2) = 0.2; AT skew = (1−1)/(1+1) = 0.0.
For `GGCC` with window 2, step 2: windows `GG` (skew +1, GC% 100) and `CC` (skew −1, GC% 100); GcSkewVariance = ((1−0)²+(−1−0)²)/2 = 1.0; GcContentVariance = 0.0 [1][3][7].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GcSkewCalculator_AnalyzeGcContent_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculator_AnalyzeGcContent_Tests.cs) — covers INV-01..INV-05
- Evidence: [SEQ-GC-ANALYSIS-001-Evidence.md](../../../docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md)
- Related algorithms: [AT_Skew](./AT_Skew.md)

## 8. References

1. Lobry JR. 1996. Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution* 13(5):660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
2. Grigoriev A. 1998. Analyzing genomes with cumulative skew diagrams. *Nucleic Acids Research* 26(10):2286–2290. https://doi.org/10.1093/nar/26.10.2286
3. Madigan MT, Martinko JM. 2003. *Brock Biology of Microorganisms*, 10th ed. Pearson-Prentice Hall. (GC% formula via https://en.wikipedia.org/wiki/GC-content)
4. Wikipedia. GC skew. https://en.wikipedia.org/wiki/GC_skew
5. Biopython contributors. Bio.SeqUtils package v1.84 (`GC_skew`, `gc_fraction`). https://biopython.org/docs/1.84/api/Bio.SeqUtils.html
6. Charneski CA, Honti F, Bryant JM, Hurst LD, Feil EJ. 2011. Atypical AT skew in Firmicute genomes. *PLoS Genetics* 7(9):e1002283. https://doi.org/10.1371/journal.pgen.1002283
7. Population Variance — Definition, Formula, Examples. Cuemath. https://www.cuemath.com/data/population-variance/
