# GC Content Profile

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-GC-PROFILE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

The GC content profile reports the local guanine+cytosine content along a nucleotide
sequence by evaluating GC content on a sliding window of fixed width that advances by a
fixed step. Each window's value is the GC-content percentage `(G + C) / (A + T + G + C) × 100`
[1]. It is an exact, specification-driven statistic — no estimation or heuristic — used to
visualise compositional structure (GC-rich/poor regions, isochores, GC islands) and as a
precursor to replication-origin and gene-finding analyses [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

GC content is the proportion of a nucleic-acid sequence composed of guanine and cytosine.
Because G·C base pairs are joined by three hydrogen bonds versus two for A·T, GC content
correlates with duplex stability and varies systematically across genomes and along
chromosomes [1]. Evaluating it on a sliding window turns the single scalar into a
positional profile of local composition.

### 2.2 Core Model

For a window of bases, GC content as a percentage is [1]:

```
GC% = (G + C) / (A + T + G + C) × 100
```

where `G`, `C`, `A`, `T` are counts of the respective standard bases in the window; for RNA,
`U` is counted in the same place as `T` (a non-GC base) [2]. The denominator is the count of
standard bases A+T+G+C only; ambiguous symbols such as `N` are excluded — this is the
denominator used by Biopython `gc_fraction` under its default `ambiguous="remove"` mode,
where `gc_fraction("ACTGN") = 0.50` (N removed from the length) [2]. Biopython's
`gc_fraction` returns the same quantity as a fraction in [0, 1]; multiplying by 100 yields
the percentage used here, consistent with Biopython `GC123`, which returns "percentages
between 0 and 100" [2].

The profile applies this formula to each window at offsets `0, step, 2·step, …` up to
`length − windowSize`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each window value lies in [0, 100] | numerator G+C ≤ denominator A+T+G+C, then ×100 [1] |
| INV-02 | Denominator counts only A/T/U/G/C; N and other symbols are excluded | Biopython default `remove`: `gc_fraction("ACTGN")=0.50` [2] |
| INV-03 | Window count = ⌊(n − w)/step⌋ + 1 for w ≤ n, else 0; offsets are 0, step, 2·step, … | sliding-window enumeration [1] |
| INV-04 | U is treated as a non-GC base equivalent to T | Biopython RNA `gc_fraction("GGAUCUUCGGAUCU")=0.50` [2] |
| INV-05 | A window with no A/T/U/G/C base yields 0 | repository zero-division convention (matches SEQ-GC-ANALYSIS-001) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | Nucleotide sequence (DNA or RNA) | case-insensitive; A/C/G/T/U standard, other symbols excluded from denominator |
| windowSize | int | 100 | Window width W in bases | windows produced only when W ≤ length |
| stepSize | int | 1 | Window advance in bases | ≥ 1 for forward progress |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | IEnumerable&lt;double&gt; | GC% (0–100) per window, in order of offset 0, step, 2·step, … |

### 3.3 Preconditions and Validation

Null or empty `sequence`, or `windowSize` greater than the sequence length, yields an empty
profile (no windows). Input is case-folded to uppercase before counting (lowercase accepted).
T and U are both non-GC bases. The denominator counts only the standard bases A/T/U/G/C in
the window; ambiguous symbols (e.g. `N`) are excluded from the denominator. A window whose
bases are all non-standard yields 0 (zero-division convention). No exceptions are thrown for
these input classes.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return empty if the sequence is null/empty or `windowSize` exceeds its length.
2. Case-fold the sequence to uppercase.
3. For each offset `i = 0, step, 2·step, …` while `i ≤ length − windowSize`:
   1. Count `gc` (G or C) and `total` (A, T, U, G, or C) over the window `[i, i+windowSize)`.
   2. Yield `gc / total × 100` if `total > 0`, else `0`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateGcContentProfile | O(W · windowSize) | O(1) streaming | W = number of windows = ⌊(n − windowSize)/step⌋ + 1; each window is recounted independently. For step ≥ windowSize windows are disjoint and total work is O(n). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateGcContentProfile(string, int, int)`: streams the GC% of each
  sliding window via deferred `yield return`.

### 5.2 Current Behavior

The method is a single short pass per window; it is not a substring-search/matching
operation, so the repository suffix tree is **not** applicable (no occurrence enumeration —
it counts G/C/A/T/U characters in place). Counting is per-window and independent; values are
produced lazily (deferred execution), so callers materialise with `ToList()`/`ToArray()`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `GC% = (G + C) / (A + T + G + C) × 100` per window [1].
- Denominator excludes ambiguous symbols (N), matching Biopython `gc_fraction` default
  `ambiguous="remove"` [2].
- U treated as a non-GC base equivalent to T [2].

**Intentionally simplified:**

- (none).

**Not implemented:**

- Ambiguity-weighted GC counting (Biopython `ambiguous="weighted"`/`"ignore"` modes);
  **users should rely on:** the default `remove` convention only — degenerate IUPAC codes
  beyond N are not weighted.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | GC content reported as percentage (0–100), not Biopython's [0,1] fraction | Deviation | value differs from `gc_fraction` by ×100 | accepted | Wikipedia ×100 form [1]; repository convention (SEQ-GC-ANALYSIS-001); locked by tests |
| 2 | Empty-window (no standard base) → 0 | Assumption | degenerate all-N window | accepted | zero-division convention, consistent with SEQ-GC-ANALYSIS-001 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty sequence | empty profile | no windows (§3.3) |
| windowSize > length | empty profile | no full window (INV-03) |
| windowSize == length | one window = whole-sequence GC% | INV-03 |
| window with N | N excluded from denominator | INV-02 [2] |
| all-N window | 0 | INV-05 (zero-division convention) |
| lowercase input | same as uppercase | case-folded counting (§3.3) |

### 6.2 Limitations

Each window is recounted from scratch (no incremental sliding sum), so cost scales with
`windowSize` per overlapping window. Only N-style exclusion is modelled; other IUPAC
degenerate codes are treated as non-standard (excluded from the denominator) rather than
fractionally weighted. The profile reports composition only — it does not classify GC
islands or isochores.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
// GGGAAATGCC, window 4, step 3 → windows GGGA, AAAT, TGCC
var profile = SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", windowSize: 4, stepSize: 3).ToList();
// profile == [75.0, 0.0, 75.0]
```

**Numerical / biological walk-through:**

`GGGAAATGCC` has windows at offsets 0, 3, 6. `GGGA`: G+C = 3, A+T+G+C = 4 → 3/4×100 = 75.0.
`AAAT`: G+C = 0 → 0.0. `TGCC` (T,G,C,C): G+C = 3, total = 4 → 75.0. Tests pin the exact
offsets and values from the formula [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateGcContentProfile_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateGcContentProfile_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [SEQ-GC-PROFILE-001-Evidence.md](../../../docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md)
- Related algorithms: [Entropy_Profile](../Statistics/Entropy_Profile.md); [Comprehensive_GC_Analysis](../Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md)

## 8. References

1. Wikipedia. 2026. GC-content. https://en.wikipedia.org/wiki/GC-content (accessed 2026-06-14).
2. Cock P.J.A. et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25(11):1422–1423. https://doi.org/10.1093/bioinformatics/btp163 ; `Bio.SeqUtils.gc_fraction` source: https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py (accessed 2026-06-14).
