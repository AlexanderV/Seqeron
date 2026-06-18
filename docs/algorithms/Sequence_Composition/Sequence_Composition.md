# Sequence Composition

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition / Statistics |
| Test Unit ID | SEQ-COMPOSITION-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Sequence composition counts the constituent symbols of a biological sequence and derives the
standard summary statistics from those counts. For nucleotide sequences this is the per-base
tally (A, T, G, C, U, N, other), GC/AT content and GC/AT compositional skew; for protein
sequences it is the per-residue tally over the 20 standard amino acids. The computation is exact
and specification-driven (no heuristics): it is a single linear pass over the sequence followed by
closed-form ratios [1][2][4].

> **Consolidation note.** This unit's two methods are identical to those already delivered under
> Test Unit **SEQ-STATS-001** (`SequenceStatistics.CalculateNucleotideComposition`,
> `CalculateAminoAcidComposition`). SEQ-COMPOSITION-001 is a duplicate Registry entry; it is
> resolved by consolidation, not re-implementation. See the SEQ-STATS-001 algorithm doc
> ([Sequence_Composition_Statistics.md](./Sequence_Composition_Statistics.md)) for the full
> treatment of GC content and skew.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Base/residue composition is the most elementary descriptor of a sequence. GC content correlates
with genome stability and taxonomy; GC and AT skew measure strand compositional asymmetry, first
reported by Lobry (1996) for bacterial genomes and used to locate replication origins [3].

### 2.2 Core Model

For a nucleotide sequence over the canonical alphabet {A, C, G, T, U} [4], with counts
`nX` of each base and `total = nA + nT + nG + nC + nU`:

- GC content = (nG + nC) / total [1]
- AT content = (nA + nT + nU) / total
- GC skew = (nG − nC) / (nG + nC) [2][3]
- AT skew = (nA − nT) / (nA + nT) [2]

For a protein sequence, composition is the multiset of single-letter amino-acid codes over the
20 standard residues {A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y} [4]; `Length` is the number of
counted residues.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ GcContent ≤ 1, 0 ≤ AtContent ≤ 1 | ratio of a subset count to the canonical-base total [1] |
| INV-02 | −1 ≤ GcSkew ≤ 1, −1 ≤ AtSkew ≤ 1 | difference over sum of non-negative counts [2] |
| INV-03 | CountA+T+G+C+U+N+Other = Length | each character is classified into exactly one bucket |
| INV-04 | amino-acid `Length` = Σ Counts | every counted residue increments exactly one entry [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | DNA/RNA or protein sequence | case-insensitive; null/empty allowed |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| NucleotideComposition | record struct | Length, per-base counts, GcContent, AtContent, GcSkew, AtSkew |
| AminoAcidComposition | record struct | Length, residue Counts, plus derived physico-chemical fields (covered by SEQ-MW/PI/HYDRO units) |

### 3.3 Preconditions and Validation

Null or empty input returns an all-zero composition (no exception) [1]. Counting is
case-insensitive (input is upper-cased) [1]. Non-canonical nucleotide letters are routed to
`CountN` (for `N`) or `CountOther`; they do not contribute to GC/AT totals.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return the zero composition for null/empty input.
2. Single pass: upper-case and classify each character into A/T/G/C/U/N/other.
3. Compute `total`, GC/AT content and GC/AT skew with zero-denominator guards.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateNucleotideComposition | O(n) | O(1) | single pass; fixed-size counters |
| CalculateAminoAcidComposition | O(n) | O(σ) | σ = distinct residues ≤ 26 |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateNucleotideComposition(string)`: per-base counts + GC/AT content + GC/AT skew.
- `SequenceStatistics.CalculateAminoAcidComposition(string)`: per-residue counts + Length.

### 5.2 Current Behavior

No search/matching is involved (a single linear scan of counters), so the repository suffix tree
is not applicable here. GC/AT skew use zero-denominator guards returning 0.0, matching Biopython's
`ZeroDivisionError` handling [1].

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- GC content = (G+C)/total, GC skew = (G−C)/(G+C), AT skew = (A−T)/(A+T) [1][2][3].
- Case-insensitive counting [1]; empty/null → zero composition [1].
- Amino-acid composition over the 20 IUPAC single-letter codes [4].

**Intentionally simplified:**

- Degenerate IUPAC codes: Biopython counts `S` toward GC and `W` toward the denominator; this
  implementation routes all non-{A,T,G,C,U} letters to `CountN`/`CountOther`. **Consequence:**
  results agree exactly with Biopython over the canonical alphabet; they differ only for sequences
  containing degenerate symbols.

**Not implemented:**

- Weighted ambiguous-base GC (`ambiguous="weighted"`); **users should rely on:** Biopython
  `gc_fraction` for degenerate-heavy sequences.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| empty / null | all-zero composition | Biopython empty handling [1] |
| no G/C | GcSkew = 0 | zero-denominator guard [1] |
| no A/T | AtSkew = 0 | zero-denominator guard [1] |
| lowercase / mixed case | same as uppercase | case-insensitive counting [1] |

### 6.2 Limitations

Degenerate IUPAC codes are not folded into GC/AT totals (see 5.3). Physico-chemical fields of the
amino-acid composition (molecular weight, pI, hydrophobicity) are owned by the SEQ-MW/PI/HYDRO
units, not this one.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
var comp = SequenceStatistics.CalculateNucleotideComposition("ATGC");
// comp.GcContent == 0.5, comp.GcSkew == 0.0, comp.Length == 4
```

For `GGGC`: GcContent = 4/4 = 1.0, GcSkew = (3−1)/4 = 0.5 [2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateNucleotideComposition_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateNucleotideComposition_Tests.cs) — covers `INV-01`–`INV-04` (shared canonical fixture with SEQ-STATS-001)
- Evidence: [SEQ-COMPOSITION-001-Evidence.md](../../../docs/Evidence/SEQ-COMPOSITION-001-Evidence.md)
- Related algorithms: [Sequence_Composition_Statistics.md](./Sequence_Composition_Statistics.md) (SEQ-STATS-001 — same methods)

## 8. References

1. Cock, P. J. A. et al. Biopython, `Bio.SeqUtils` (`gc_fraction`, `GC_skew`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py (accessed 2026-06-14)
2. Wikipedia contributors. GC skew. https://en.wikipedia.org/wiki/GC_skew (accessed 2026-06-14)
3. Lobry, J. R. 1996. Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution* 13(5):660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
4. IUPAC nucleotide and amino-acid single-letter codes. https://www.bioinformatics.org/sms/iupac.html (accessed 2026-06-14)
