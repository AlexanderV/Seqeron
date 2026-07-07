# Sequence Summary

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-SUMMARY-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

`SummarizeNucleotideSequence` produces a single `SequenceSummary` record that bundles the
most common descriptive statistics of a DNA/RNA sequence: length, GC content, Shannon
entropy, linguistic complexity, melting temperature, and a per-symbol nucleotide
composition. It performs no new computation of its own; it is a pure aggregation that
delegates each field to the corresponding already-defined canonical metric method. It is
exact in the sense that each field reproduces, bit-for-bit, the value the underlying metric
method returns on the same input.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Sequence "summary" / "stats" records are a standard convenience in sequence-analysis
toolkits: rather than calling several functions, a caller gets one struct with the
headline descriptors. Each descriptor has its own formal basis (composition counting,
Shannon information entropy, Trifonov linguistic complexity, oligo melting temperature);
the summary's only formal obligation is field-wise consistency with those metrics.

### 2.2 Core Model

For an input sequence `S`, the summary fields are defined as:

- **Length** = `|S|` (raw character count).
- **GcContent** = GC fraction = (#G + #C) / (#counted bases), case-insensitive, 0 for empty input [1].
- **Entropy** = Shannon entropy `H = − Σ p·log₂ p` over the per-symbol frequencies, in bits [2].
- **Complexity** = linguistic complexity, a vocabulary-usage measure (observed vs possible words) combined across word sizes, in the range (0,1) [3].
- **MeltingTemperature** = Wallace rule `2(A+T) + 4(G+C)` for short oligos, otherwise the GC/Marmur-Doty formula `64.9 + 41·(GC − 16.4)/N` [4][5].
- **Composition** = the counts of A, T, G, C, U, N [1].

The summary selects the Wallace branch when `|S| < 14` and the GC branch otherwise [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `summary.Length == S.Length` | field copies `CalculateNucleotideComposition(S).Length` |
| INV-02 | `summary.GcContent == CalculateNucleotideComposition(S).GcContent` | field is read directly from the composition record [1] |
| INV-03 | `summary.Entropy == CalculateShannonEntropy(S)` | field is the return of that method [2] |
| INV-04 | `summary.Complexity == CalculateLinguisticComplexity(S)` | field is the return of that method [3] |
| INV-05 | `summary.MeltingTemperature == CalculateMeltingTemperature(S, S.Length < 14)` | field is the return of that method with that flag [4][5] |
| INV-06 | Composition dict A,T,G,C,U,N counts equal `CalculateNucleotideComposition(S)` counts | dict is built directly from those counts [1] |
| INV-07 | 0 ≤ GcContent ≤ 1 and 0 ≤ Complexity < 1 (DNA fragments) | fraction and vocabulary-usage bounds [1][3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | DNA/RNA sequence | case-insensitive; null/empty allowed (degenerate result) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Length` | `int` | character count of the input |
| `GcContent` | `double` | GC fraction in [0,1] |
| `Entropy` | `double` | Shannon entropy in bits |
| `Complexity` | `double` | linguistic complexity in (0,1) |
| `MeltingTemperature` | `double` | Tm in °C |
| `Composition` | `IReadOnlyDictionary<char,int>` | counts for A,T,G,C,U,N |

### 3.3 Preconditions and Validation

Null or empty input returns a degenerate summary (Length 0, GcContent 0, Entropy 0,
Complexity 0, MeltingTemperature 0, all composition counts 0); no exception is thrown,
matching the empty-sequence handling of each per-metric method [1]. Input is
case-insensitive (each per-metric method uppercases internally). T and U are counted as
distinct symbols; N is counted; other characters are excluded from GC/entropy as defined
by the per-metric methods.

## 4. Algorithm

### 4.1 High-Level Steps

1. Compute `comp = CalculateNucleotideComposition(sequence)`.
2. Compute `entropy = CalculateShannonEntropy(sequence)`.
3. Compute `complexity = CalculateLinguisticComplexity(sequence)`.
4. Compute `tm = CalculateMeltingTemperature(sequence, useWallaceRule: sequence.Length < 14)`.
5. Build the composition dictionary from the composition counts and assemble the record.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `SummarizeNucleotideSequence` | O(n) | O(n) | dominated by linguistic complexity (k-mer sets up to k=6); other metrics are O(n) single passes |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.SummarizeNucleotideSequence(string)`: aggregates the metrics into a `SequenceSummary`.
- `SequenceStatistics.SequenceSummary`: the returned `readonly record struct`.

### 5.2 Current Behavior

The method is a thin aggregator: it calls the four canonical metric methods and copies
their results into the record, plus a 6-entry composition dictionary (A,T,G,C,U,N) built
from the composition counts. No string searching/matching is performed by the summary
itself, so the repository suffix tree is **not** applicable here (N/A — this is a counting
and arithmetic aggregation, not occurrence enumeration). The underlying
`CalculateLinguisticComplexity` builds short k-mer sets directly; it does not query the
suffix tree, and the summary does not change that.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- GcContent as GC fraction over counted bases [1]; Shannon entropy in bits H = −Σ p·log₂ p [2]; melting temperature Wallace / GC-Marmur-Doty branch selection [4][5]; composition counts A,T,G,C,U,N [1].
- Field-wise equality: every summary field equals its canonical per-metric method's value on the same input (the aggregation contract).

**Intentionally simplified:**

- Complexity: the aggregated `CalculateLinguisticComplexity` computes the **mean** of per-word-size vocabulary-usage ratios rather than the **product** `C = U₁U₂…Uw` defined by Trifonov [3]; **consequence:** the Complexity value differs from a strict Trifonov product. This is a property of the linguistic-complexity method (its own unit), not of the aggregation; the summary faithfully reports that method's value.

**Not implemented:**

- (none) — the summary adds no metric beyond what the per-metric methods provide.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Tm threshold length<14 | Assumption | selects Wallace vs GC formula | accepted | sibling SEQ-TM-001 convention (`ThermoConstants.WallaceMaxLength`); summary tested for equality with `CalculateMeltingTemperature` |
| 2 | Complexity = mean (not Trifonov product) | Deviation | Complexity value differs from strict Trifonov | accepted | belongs to the linguistic-complexity method; out of scope for the aggregation |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| empty / null sequence | degenerate summary (all zero counts/metrics) | per-metric empty handling [1] |
| lowercase input | identical summary to uppercase | per-metric methods uppercase internally |
| RNA input (U) | U counted; GC/entropy include U as a symbol | composition counts U [1] |
| length exactly 14 | GC/Marmur-Doty branch (14 is not < 14) | threshold is strict `<` |

### 6.2 Limitations

The summary inherits every limitation of its component metrics (e.g., the mean-based
linguistic complexity, the simple Tm formulas without salt/nearest-neighbor correction).
It is descriptive only and does not validate that the input is biologically meaningful.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var s = SequenceStatistics.SummarizeNucleotideSequence("ATGCATGC");
// s.Length == 8, s.GcContent == 0.5, s.Entropy == 2.0 (log2 4),
// s.MeltingTemperature == 24.0 (Wallace: 2*(A+T)+4*(G+C) = 2*4+4*4)
```

**Numerical walk-through:** for "ATGCATGC" the composition is A=2,T=2,G=2,C=2 →
GcContent = 4/8 = 0.5; four equally frequent symbols → H = log₂ 4 = 2.0 bits; length 8 < 14
→ Wallace Tm = 2·(2+2) + 4·(2+2) = 8 + 16 = 24.0 °C.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_SummarizeNucleotideSequence_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_SummarizeNucleotideSequence_Tests.cs) — covers `INV-01`..`INV-07`
- Evidence: [SEQ-SUMMARY-001-Evidence.md](../../../docs/Evidence/SEQ-SUMMARY-001-Evidence.md)
- Related algorithms: [Entropy_Profile](../Statistics/Entropy_Profile.md), [Melting_Temperature](../Statistics/Melting_Temperature.md)

## 8. References

1. Cock, P. J. A. et al. 2009. Biopython (`Bio.SeqUtils.gc_fraction`). Bioinformatics 25(11):1422–1423. https://doi.org/10.1093/bioinformatics/btp163 (source: https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py).
2. Shannon, C. E. 1948. A Mathematical Theory of Communication. Bell System Technical Journal 27(3):379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x (formula/units: https://en.wikipedia.org/wiki/Entropy_(information_theory)).
3. Trifonov, E. N. 1990. Linguistic sequence complexity (vocabulary usage). https://en.wikipedia.org/wiki/Linguistic_sequence_complexity.
4. Biopython `Bio.SeqUtils.MeltingTemp` (`Tm_Wallace`, `Tm_GC`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py.
5. Marmur, J., Doty, P. 1962. Determination of the base composition of DNA from its thermal denaturation temperature. J Mol Biol 5:109–118.
