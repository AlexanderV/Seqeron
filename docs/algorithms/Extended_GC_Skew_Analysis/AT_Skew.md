# AT Skew

| Field | Value |
|-------|-------|
| Algorithm Group | Extended GC Skew Analysis (Composition) |
| Test Unit ID | SEQ-ATSKEW-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

AT skew measures the compositional asymmetry between adenine and thymine within a single DNA strand as `(A − T) / (A + T)` [1][2]. Together with GC skew it is used to detect strand-specific mutational and selective bias and, by its cumulative form, to locate replication origins and termini in bacterial genomes [1][3]. This unit computes the single scalar AT skew of a whole sequence (or window). It is an exact arithmetic statistic, not a heuristic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In a strand at intrastrand equilibrium ("second parity rule") A ≈ T and C ≈ G. Lobry (1996) showed that the two strands of bacterial DNA depart from this equifrequency, the leading strand being relatively enriched in G over C and T over A [1]. AT skew quantifies the A-vs-T part of that departure for one given strand [1][2].

### 2.2 Core Model

For a sequence with adenine count `A` and thymine count `T`:

```
AT skew = (A − T) / (A + T)
```

verbatim from Charneski et al. (2011) [2] and corroborated by the Lobry (1996) primary as cited in [3]. Only A and T contribute; G, C and any other symbol are not counted [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Result ∈ [−1, +1] | numerator magnitude `|A − T| ≤ A + T` = denominator [3] |
| INV-02 | A = T (and A + T > 0) ⇒ result = 0 | numerator A − T = 0 [2] |
| INV-03 | A + T = 0 ⇒ result = 0 (no exception/NaN) | zero-denominator convention from Biopython `GC_skew` [4] |
| INV-04 | Result is case-insensitive | counting ignores case [4]; input normalized to upper case |
| INV-05 | Non-A/T symbols do not change the value | only A and T are counted [4] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | AT skew | GC skew |
|--------|---------|---------|
| Bases compared | A vs T | G vs C [1] |
| Formula | (A − T)/(A + T) [2] | (G − C)/(G + C) [3] |
| Primary use | strand A/T asymmetry; complements GC skew | replication origin/terminus location [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `string` | required | DNA sequence; only A/T counted | case-insensitive; null/empty ⇒ 0 |
| sequence | `DnaSequence` | required | DNA sequence value object (already upper-cased) | non-null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `double` | AT skew in [−1, +1]; 0 when A + T = 0 |

### 3.3 Preconditions and Validation

- String overload: `null` or empty ⇒ returns `0`. Input is upper-cased (`ToUpperInvariant`) before counting, so case does not matter.
- `DnaSequence` overload: `null` ⇒ `ArgumentNullException`; the value object is constructed already normalized to upper case.
- Accepted alphabet: A and T are counted; every other symbol (G, C, N, gaps, IUPAC ambiguity) is ignored and affects neither numerator nor denominator [4]. No T↔U normalization is performed (DNA input expected).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input (null/empty handling per overload).
2. Count occurrences of `A` and `T` (case-insensitive).
3. If `A + T = 0`, return 0; otherwise return `(A − T) / (A + T)`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateAtSkew | O(n) | O(1) | two linear passes counting A and T; n = sequence length |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GcSkewCalculator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs)

- `GcSkewCalculator.CalculateAtSkew(string)`: canonical entry; upper-cases input, counts A/T, returns the skew.
- `GcSkewCalculator.CalculateAtSkew(DnaSequence)`: delegate; forwards the normalized `DnaSequence.Sequence` to the same core.

### 5.2 Current Behavior

The two public overloads share a private `CalculateAtSkewCore`. Counting uses `string.Count` over the characters 'A' and 'T'. The suffix tree was **not** evaluated/used: AT skew is a single linear count of two base symbols, not a substring-search or occurrence-enumeration problem, so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `(A − T) / (A + T)` exactly as defined in Charneski et al. (2011) [2] and Lobry (1996) [1][3].
- Range [−1, +1] (INV-01) [3].
- Zero-denominator ⇒ 0 (INV-03), matching the Biopython `GC_skew` ZeroDivisionError → 0.0 convention [4].
- Only A/T counted; other symbols ignored (INV-05) [4].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Windowed / cumulative AT skew profiles and replication-origin location from AT skew: out of scope for this unit; **users should rely on** `CalculateWindowedGcSkew` / `CalculateCumulativeGcSkew` / `PredictReplicationOrigin` (GC-skew based) in the same class.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Case-insensitivity and "ignore non-A/T symbols" for the AT-skew analog | Assumption | Determines counts for lowercase/ambiguous input | accepted | Taken from the directly analogous Biopython `GC_skew` [4]; formula itself fully sourced [1][2] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` / empty string | 0 | documented validation |
| `null` DnaSequence | `ArgumentNullException` | documented validation |
| no A and no T (e.g. "GGCC") | 0 | A + T = 0 (INV-03) [4] |
| pure A | +1.0 | T = 0 (INV-01) [3] |
| pure T | −1.0 | A = 0 (INV-01) [3] |
| lowercase input | same as upper case | case-insensitive (INV-04) [4] |

### 6.2 Limitations

Computes a single global statistic; it does not localize asymmetry along the sequence (use the windowed/cumulative GC-skew methods for that). No T↔U conversion: RNA input ("U") would be ignored rather than treated as T.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
double skew = GcSkewCalculator.CalculateAtSkew("AAAT"); // 0.5
```

**Numerical walk-through:** `"AAATGGGCCC"` → A = 3, T = 1, G/C ignored ⇒ (3 − 1)/(3 + 1) = 0.5.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GcSkewCalculator_CalculateAtSkew_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculator_CalculateAtSkew_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [SEQ-ATSKEW-001-Evidence.md](../../../docs/Evidence/SEQ-ATSKEW-001-Evidence.md)

## 8. References

1. Lobry, J. R. 1996. Asymmetric substitution patterns in the two DNA strands of bacteria. Molecular Biology and Evolution 13(5):660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
2. Charneski, C. A., Honti, F., Bryant, J. M., Hurst, L. D., Feil, E. J. 2011. Atypical AT Skew in Firmicute Genomes Results from Selection and Not from Mutation. PLoS Genetics 7(9):e1002283. https://doi.org/10.1371/journal.pgen.1002283
3. Wikipedia contributors. GC skew. https://en.wikipedia.org/wiki/GC_skew (accessed 2026-06-14).
4. Biopython project. Bio.SeqUtils (GC_skew). https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/__init__.py (accessed 2026-06-14).
