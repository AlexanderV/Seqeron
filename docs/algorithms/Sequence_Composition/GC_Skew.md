# GC Skew

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-GCSKEW-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

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
| `[PredictReplicationOrigin/AnalyzeGcContent DnaSequence] windowSize` | `int` | `1000` | Window length used by higher-level helpers | Positive values are expected, but the current helpers do not validate them before calling the core loops |
| `[AnalyzeGcContent DnaSequence] stepSize` | `int` | `100` | Step size used by the comprehensive analysis helper | Positive values are expected, but the current helper does not validate them before calling the core loops |
| `[string] windowSize` | `int` | `1000` | Sliding-window or cumulative window length | Nonpositive values are not validated and are currently unsupported |
| `[string] stepSize` | `int` | `100` | Step size for windowed GC skew | Nonpositive values are not validated and are currently unsupported |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `skew` | `double` | Overall GC skew for the sequence |
| `GcSkewPoint` | record | Window center, skew value, window start, and window end |
| `CumulativeGcSkewPoint` | record | Window center, window skew, and cumulative skew |
| `ReplicationOriginPrediction` | record | Predicted origin and terminus positions, their cumulative-skew values, and a significance flag |

### 3.3 Preconditions and Validation

`CalculateGcSkew(...)` returns `0` for empty string input and throws `ArgumentNullException` for null `DnaSequence` input. The typed `CalculateWindowedGcSkew(...)` and `CalculateCumulativeGcSkew(...)` overloads throw `ArgumentOutOfRangeException` when `windowSize < 1` or `stepSize < 1`. The raw-string overloads only guard empty input and otherwise delegate directly to the core loops, so nonpositive sizes are currently unsupported rather than consistently validated; in particular, `stepSize = 0` in `CalculateWindowedGcSkew(string, ...)` and `windowSize = 0` in `CalculateCumulativeGcSkew(string, ...)` can produce non-terminating enumerations. The higher-level typed helpers `PredictReplicationOrigin(...)` and `AnalyzeGcContent(...)` also call the unvalidated core routines directly, so they likewise assume positive `windowSize` and `stepSize` without enforcing them. `PredictReplicationOrigin(...)` returns a zeroed prediction when no cumulative points are available.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count `G` and `C` bases in the full sequence or current window.
2. Compute `(G - C) / (G + C)` and return `0` when the denominator is zero.
3. For sliding-window analysis, emit the skew at each window center.
4. For cumulative skew, sum each window's skew value across the traversal.
5. For origin prediction, choose the cumulative-skew minimum as the origin and the maximum as the terminus.

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
| `CalculateWindowedGcSkew` | `O(n)` | `O(1)` streaming | Emits one point per valid window |
| `CalculateCumulativeGcSkew` | `O(n)` | `O(1)` streaming | Uses fixed-size windows |
| `PredictReplicationOrigin` | `O(n)` | `O(w)` | Materializes cumulative points to find extrema |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GcSkewCalculator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs)

- `GcSkewCalculator.CalculateGcSkew(...)`: Computes whole-sequence GC skew.
- `GcSkewCalculator.CalculateWindowedGcSkew(...)`: Computes windowed skew with configurable step size.
- `GcSkewCalculator.CalculateCumulativeGcSkew(...)`: Produces cumulative skew points.
- `GcSkewCalculator.PredictReplicationOrigin(...)`: Predicts origin and terminus positions from cumulative skew.

### 5.2 Current Behavior

Windowed GC skew reports positions at the center of each analyzed window. Cumulative GC skew uses non-overlapping windows because the source sets `stepSize = windowSize` inside the cumulative routine. `PredictReplicationOrigin(...)` flags the result as significant only when the cumulative-skew amplitude exceeds `0.01 × pointCount`. The raw-string windowed and cumulative overloads uppercase input but do not validate nonpositive window or step values before entering the shared loops. The higher-level typed helpers `PredictReplicationOrigin(...)` and `AnalyzeGcContent(...)` also bypass the validated typed windowed/cumulative entry points and call the same core loops directly. The same class also provides `CalculateAtSkew(...)` and a combined `AnalyzeGcContent(...)` helper.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Standard GC skew calculation.
- Sliding-window and cumulative-skew analysis.
- Origin/terminus prediction from cumulative-skew extrema.

**Intentionally simplified:**

- Origin prediction assumes the cumulative-skew minimum and maximum map directly to ori/ter; **consequence:** more complex replication architectures are not modeled.
- Significance is determined by a fixed amplitude heuristic; **consequence:** predictions depend on the source-specific threshold rather than a statistical test.

**Not implemented:**

- Correction for horizontal gene transfer, inversions, or other genome-history effects; **users should rely on:** downstream comparative analysis when those effects matter.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns `0` or yields no points | Explicit source guard |
| No `G` or `C` bases | Returns `0` | Division-by-zero protection |
| Invalid window or step size on validated typed windowed/cumulative overloads | Throws `ArgumentOutOfRangeException` | Explicit validation exists only in those typed overloads |
| Nonpositive `windowSize`/`stepSize` in `PredictReplicationOrigin(...)` or `AnalyzeGcContent(...)` | Unsupported; the helpers do not validate before calling the core loops | These methods bypass the validated typed windowed/cumulative entry points |
| Nonpositive window or step size on raw-string windowed/cumulative overloads | Unsupported; zero values can fail to terminate | Those overloads delegate directly to the core loops without parameter validation |

### 6.2 Limitations

Origin and terminus prediction assume a single circular chromosome with bidirectional replication, use a heuristic significance threshold, and do not account for genome rearrangements or horizontal transfer that may distort the skew profile. The raw-string windowed and cumulative overloads also do not consistently validate nonpositive `windowSize` or `stepSize` values.

## 8. References

1. Lobry, J.R. (1996). "Asymmetric substitution patterns in the two DNA strands of bacteria." *Molecular Biology and Evolution*, 13(5):660-665.
2. Grigoriev, A. (1998). "Analyzing genomes with cumulative skew diagrams." *Nucleic Acids Research*, 26(10):2286-2290.
3. Tillier, E.R. & Collins, R.A. (2000). "The contributions of replication orientation, gene direction, and signal sequences to base-composition asymmetries in bacterial genomes." *Journal of Molecular Evolution*, 50:249-257.
4. Wikipedia contributors. "GC skew." *Wikipedia, The Free Encyclopedia*.
