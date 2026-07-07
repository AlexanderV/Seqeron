# Sequence Composition Statistics

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-STATS-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

This unit computes single-pass composition statistics for a DNA/RNA sequence: per-base counts (A, T, G, C, U, N, other), GC content, AT content, and the strand-asymmetry measures GC skew and AT skew. It also exposes a protein-composition counterpart and a nucleotide summary aggregator. The computations are exact, deterministic, specification-driven formulas (no heuristics) [1][2][3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

GC content is the fraction of a sequence that is guanine or cytosine; it correlates with thermal stability and varies by organism and genomic region. Strand compositional asymmetry — measured as GC skew and AT skew — was first reported by Lobry (1996), who found "a departure from intrastrand equifrequency between A and T or between C and G," and is used to locate replication origins/termini [1][3].

### 2.2 Core Model

For counts of each base in the sequence:

- **GC content** = (G + C) / (A + T + G + C + U) — Biopython `gc_fraction` computes G+C over the standard-alphabet length, returning a float in [0, 1] [2].
- **AT content** = (A + T + U) / (A + T + G + C + U).
- **GC skew** = (G − C) / (G + C); 0 when G + C = 0 [2][3].
- **AT skew** = (A − T) / (A + T); 0 when A + T = 0 [3].

Positive GC skew indicates G-richness, negative indicates C-richness [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ GcContent ≤ 1, 0 ≤ AtContent ≤ 1 | Ratio of non-negative subset count to total [2] |
| INV-02 | −1 ≤ GcSkew ≤ 1, −1 ≤ AtSkew ≤ 1 | (x−y)/(x+y) with x,y ≥ 0 lies in [−1,1] [3] |
| INV-03 | CountA+CountT+CountG+CountC+CountU+CountN+CountOther = Length | Counts partition every character |
| INV-04 | GcSkew = 0 if G+C = 0; AtSkew = 0 if A+T = 0 | Zero-denominator guard [2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | DNA/RNA sequence | null/empty allowed; case-insensitive; standard alphabet {A,T,G,C,U}, others → N/Other |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Length | int | Total characters in the input |
| CountA/T/G/C/U/N/Other | int | Per-base counts; non-ACGTU letters and symbols counted as N or Other |
| GcContent / AtContent | double | Fractions in [0,1] |
| GcSkew / AtSkew | double | Strand asymmetry in [−1,1] |

### 3.3 Preconditions and Validation

Null or empty input returns an all-zero `NucleotideComposition`. Input is upper-cased (`ToUpperInvariant`) before counting, so the result is case-insensitive. T and U are counted separately; degenerate IUPAC codes (S, W, R, …) are not in the standard alphabet and are counted as `Other`. No exceptions are thrown for any string input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return all-zero composition if the input is null or empty.
2. Single pass over the upper-cased sequence, incrementing A/T/G/C/U/N/Other counts.
3. total = A+T+G+C+U; compute GcContent, AtContent.
4. Compute GcSkew = (G−C)/(G+C) and AtSkew = (A−T)/(A+T), each guarded against a zero denominator.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateNucleotideComposition | O(n) | O(1) | one pass; `ToUpperInvariant` allocates one O(n) string |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateNucleotideComposition(string)`: counts + GC/AT content + GC/AT skew.
- `SequenceStatistics.SummarizeNucleotideSequence(string)`: delegates to the above plus entropy/complexity/Tm.
- `SequenceStatistics.CalculateAminoAcidComposition(string)`: protein residue counts/ratios (MW/pI/hydrophobicity belong to SEQ-MW/PI/HYDRO units).

### 5.2 Current Behavior

This is not a search/matching operation, so the repository suffix tree is not applicable — composition is a single linear scan. T and U are tracked as separate counts so both DNA and RNA inputs are handled. AT content includes U in the numerator (A+T+U), but AT skew uses the DNA-specific (A−T)/(A+T) formula without U, matching the Lobry/Wikipedia definition [1][3].

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- GC content = (G+C)/length over the standard alphabet [2].
- GC skew = (G−C)/(G+C) with 0 returned when G+C = 0 [2][3].
- AT skew = (A−T)/(A+T) with 0 returned when A+T = 0 [3].
- Case-insensitive counting [2].

**Intentionally simplified:**

- Degenerate IUPAC symbols (S, W, R, Y, …): counted as `Other` rather than partially toward GC/AT as Biopython's `gc_fraction` does for S/W; **consequence:** for sequences containing degenerate codes the GC content differs slightly from Biopython. For the standard {A,T,G,C,U} alphabet the results are identical.

**Not implemented:**

- Windowed / cumulative GC skew and origin prediction; **users should rely on:** the SEQ-GCSKEW-001 unit ([GC_Skew.md](GC_Skew.md)).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty | all-zero composition, GcContent 0 | Biopython `gc_fraction` empty → 0 [2] |
| no G/C | GcSkew 0 | zero-denominator guard [2] |
| no A/T | AtSkew 0 | zero-denominator guard |
| N/degenerate | counted as N/Other; excluded from GC/AT totals | standard-alphabet scope |
| mixed case | identical to upper-case result | `ToUpperInvariant` [2] |

### 6.2 Limitations

Does not model degenerate IUPAC ambiguity codes in GC/AT totals. Windowed skew and replication-boundary heuristics are out of scope (separate unit).

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
var c = SequenceStatistics.CalculateNucleotideComposition("GGGC");
// c.GcContent == 1.0, c.GcSkew == (3-1)/4 == 0.5
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateNucleotideComposition_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_CalculateNucleotideComposition_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [SEQ-STATS-001-Evidence.md](../../../docs/Evidence/SEQ-STATS-001-Evidence.md)
- Related algorithms: [GC_Skew](GC_Skew.md)

## 8. References

1. Lobry, J. R. 1996. Asymmetric substitution patterns in the two DNA strands of bacteria. Molecular Biology and Evolution 13(5):660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
2. Cock, P. J. A. et al. Biopython, Bio.SeqUtils (gc_fraction, GC_skew). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
3. Wikipedia contributors. GC skew. https://en.wikipedia.org/wiki/GC_skew
