# GC Skew

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-GCSKEW-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-09-28 |

## 1. Overview

GC skew measures strand-specific asymmetry between guanine and cytosine counts and is commonly used to study replication-associated composition bias. In this repository, the documented surface covers whole-sequence skew, sliding-window skew, cumulative skew, and a heuristic origin/terminus prediction based on cumulative-skew extrema. The core skew formula is exact, while the replication-boundary predictor adds implementation-specific heuristics.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

GC skew analysis is commonly used to identify replication origins and termini in bacterial and archaeal genomes. The original document notes that leading strands often show positive GC skew, lagging strands often show negative GC skew, and the skew typically changes sign near replication boundaries. It also attributes this asymmetry to strand-specific mutational pressures during replication. Sources: Lobry (1996), Grigoriev (1998), Tillier & Collins (2000), Wikipedia (GC skew).

### 2.2 Core Model

GC skew is defined as:

$$
GC\ skew = \frac{G - C}{G + C}
$$

where `G` and `C` are counts of guanine and cytosine in the analyzed region. The cumulative form used for boundary detection is:

$$
Cumulative\ GC\ skew(n) = \sum_{i=1}^{n} GC\ skew(window_i)
$$

i.e. "a sum of (G−C)/(G+C) in adjacent windows from an arbitrary start to a given point in a sequence" (Grigoriev 1998, abstract). Windows are adjacent and non-overlapping (step = window size).

The original document interprets the global minimum of cumulative skew as the replication origin and the global maximum as the terminus for typical circular bacterial chromosomes.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `-1 <= GC skew <= 1` | The numerator is bounded by the denominator in absolute value |
| INV-02 | Empty input or windows with no `G` or `C` bases yield `0` | The implementation guards against division by zero |
| INV-03 | Windowed positions are reported at `WindowStart + WindowSize / 2` | That coordinate formula is explicit in source |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | Sequence to analyze | Null `DnaSequence` input throws `ArgumentNullException`; empty string yields `0` or no points |
| `[CalculateWindowedGcSkew/CalculateCumulativeGcSkew DnaSequence] windowSize` | `int` | `1000` | Sliding-window or cumulative window length | Must be `>= 1` |
| `[CalculateWindowedGcSkew DnaSequence] stepSize` | `int` | `100` | Step size for windowed GC skew | Must be `>= 1` |
| `[string] windowSize` | `int` | `1000` | Sliding-window or cumulative window length | Must be `>= 1` (validated eagerly, same as the typed overloads) |
| `[string] stepSize` | `int` | `100` | Step size for windowed GC skew | Must be `>= 1` (validated eagerly) |

`PredictReplicationOrigin(...)` (SEQ-REPLICATION-001) takes no window parameter; `AnalyzeGcContent(...)` is covered by SEQ-GC-ANALYSIS-001.

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `skew` | `double` | Overall GC skew for the sequence |
| `GcSkewPoint` | record | Window center, skew value, window start, and window end |
| `CumulativeGcSkewPoint` | record | Window center, window skew, and cumulative skew |
| `ReplicationOriginPrediction` | record | Predicted origin and terminus positions, their cumulative-skew values, and a significance flag |

### 3.3 Preconditions and Validation

`CalculateGcSkew(...)` returns `0` for empty string input and throws `ArgumentNullException` for null `DnaSequence` input. Both the typed and the raw-string `CalculateWindowedGcSkew(...)` / `CalculateCumulativeGcSkew(...)` overloads throw `ArgumentOutOfRangeException` when `windowSize < 1` or `stepSize < 1`; validation is eager (at call time, not at first enumeration). Before the 2026-09 review the raw-string overloads did not validate, and `stepSize = 0` (windowed) / `windowSize = 0` (cumulative) produced non-terminating enumerations. Null or empty strings yield `0` / an empty sequence of points.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count `G` and `C` bases in the full sequence or current window.
2. Compute `(G - C) / (G + C)` and return `0` when the denominator is zero.
3. For sliding-window analysis, emit the skew at each window center.
4. For cumulative skew, sum each window's skew value across the traversal.
5. For origin prediction (SEQ-REPLICATION-001), build the per-nucleotide cumulative skew (G = +1, C = −1) and choose its first global minimum as the origin and first global maximum as the terminus.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Related metric documented alongside GC skew:

$$
AT\ skew = \frac{A - T}{A + T}
$$

The same source file also provides `CalculateAtSkew(...)` helpers and a combined GC-content analysis surface.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateGcSkew` | `O(n)` | `O(1)` | Single pass over the sequence |
| `CalculateWindowedGcSkew` | `O(n·w/s)` | `O(w)` streaming | Recounts each window (w = window, s = step); `O(n)` when `s = w` |
| `CalculateCumulativeGcSkew` | `O(n)` | `O(w)` streaming | Adjacent non-overlapping windows |
| `PredictReplicationOrigin` | `O(n)` | `O(1)` | Single pass over the per-nucleotide cumulative skew |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GcSkewCalculator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs)

- `GcSkewCalculator.CalculateGcSkew(...)`: Computes whole-sequence GC skew.
- `GcSkewCalculator.CalculateWindowedGcSkew(...)`: Computes windowed skew with configurable step size.
- `GcSkewCalculator.CalculateCumulativeGcSkew(...)`: Produces cumulative skew points.
- `GcSkewCalculator.PredictReplicationOrigin(...)`: Predicts origin and terminus positions from cumulative skew.

### 5.2 Current Behavior

Windowed GC skew reports positions at the center of each analyzed window. Cumulative GC skew uses non-overlapping windows because the source sets `stepSize = windowSize` inside the cumulative routine. Only complete windows are reported (window starts `0, s, 2s, …` while `start + w ≤ n`); a trailing partial window is dropped. This matches SkewIT `gcskew.py` (Lu & Salzberg 2020), whereas Biopython `Bio.SeqUtils.GC_skew` appends the partial tail window (e.g. `GC_skew("GGGGCCCCGG", 4)` = `[1.0, -1.0, 1.0]` vs Seqeron `[1.0, -1.0]`); on the complete windows the values agree exactly. Counting is case-insensitive; only `G`/`C` are counted (ambiguity codes such as `S` are ignored, as in Biopython). `PredictReplicationOrigin(...)` works on the per-nucleotide cumulative skew and flags significance when max > min (see SEQ-REPLICATION-001). The same class also provides `CalculateAtSkew(...)` and a combined `AnalyzeGcContent(...)` helper.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Standard GC skew calculation.
- Sliding-window and cumulative-skew analysis.
- Origin/terminus prediction from cumulative-skew extrema.

**Intentionally simplified:**

- Origin prediction assumes the cumulative-skew minimum and maximum map directly to ori/ter; **consequence:** more complex replication architectures are not modeled.
- Significance is only "non-zero amplitude" (max > min); **consequence:** no statistical test is applied.

**Not implemented:**

- Correction for horizontal gene transfer, inversions, or other genome-history effects; **users should rely on:** downstream comparative analysis when those effects matter.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns `0` or yields no points | Explicit source guard |
| No `G` or `C` bases | Returns `0` | Division-by-zero protection |
| `windowSize < 1` or `stepSize < 1` (typed or raw-string windowed/cumulative overloads) | Throws `ArgumentOutOfRangeException` at call time | Eager guard clause |
| Sequence shorter than the window / trailing partial window | Not reported | Complete-window convention (SkewIT); Biopython appends the partial window |

### 6.2 Limitations

Origin and terminus prediction assume a single circular chromosome with bidirectional replication, use a simple non-zero-amplitude significance flag, and do not account for genome rearrangements or horizontal transfer that may distort the skew profile.

## 8. References

1. Lobry, J.R. (1996). "Asymmetric substitution patterns in the two DNA strands of bacteria." *Molecular Biology and Evolution*, 13(5):660-665.
2. Grigoriev, A. (1998). "Analyzing genomes with cumulative skew diagrams." *Nucleic Acids Research*, 26(10):2286-2290.
3. Tillier, E.R. & Collins, R.A. (2000). "The contributions of replication orientation, gene direction, and signal sequences to base-composition asymmetries in bacterial genomes." *Journal of Molecular Evolution*, 50:249-257.
4. Wikipedia contributors. "GC skew." *Wikipedia, The Free Encyclopedia*.
5. Biopython 1.88, `Bio.SeqUtils.GC_skew` / `xGC_skew` (reference implementation, numerically cross-checked).
6. Lu, J. & Salzberg, S.L. (2020). SkewIT: The Skew Index Test for large-scale GC Skew analysis of bacterial genomes. *PLoS Comput Biol* 16(12):e1008439; `gcskew.py` (github.com/jenniferlu717/SkewIT).
