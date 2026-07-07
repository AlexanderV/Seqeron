# Consensus Sequence from a Multiple Alignment

| Field | Value |
|-------|-------|
| Algorithm Group | Matching / Motif analysis |
| Test Unit ID | MOTIF-CONS-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Given a set of equal-length aligned DNA sequences, this algorithm computes the consensus string: at each column it emits the most frequently occurring nucleotide [1][2]. It is the classical column-wise majority consensus used to summarise a multiple alignment as a single representative sequence. The computation is exact and deterministic; ties between equally-frequent bases are resolved in a fixed alphabetical order [4]. It differs from the IUPAC-degenerate consensus (`GenerateConsensus`), which encodes ambiguous columns with IUPAC codes rather than picking a single base.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A multiple sequence alignment places homologous sequences in a common coordinate frame so that each column contains positionally-corresponding residues. Summarising the alignment as one "typical" sequence — the consensus — is a standard step in motif description and primer/probe design [1].

### 2.2 Core Model

For aligned strings of length *n*, build a 4×*n* profile matrix *P* where *P*[b, j] is the number of times base *b* ∈ {A, C, G, T} occurs in column *j* [2]. The consensus *c* is defined position-wise: the *j*th symbol of *c* is "the symbol having the maximum value in the *j*-th column of the profile matrix" [2], i.e. the most frequent residue at that position [1]. When several symbols share the maximum count there may be more than one valid consensus [2]; this implementation fixes the choice to the alphabetically-earliest tied base [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output length = common input length | one symbol emitted per column [2] |
| INV-02 | Each output symbol is a base attaining the maximum count in its column | definition of consensus [2] |
| INV-03 | Ties resolve to the alphabetically-earliest base (A<C<G<T) | fixed scan order over the alphabet [4] |
| INV-04 | Identical inputs ⇒ output equals that sequence | every column is unanimous [1][2] |
| INV-05 | Deterministic for any valid input | INV-03 removes the only ambiguity [4] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Most-frequent consensus (this) | IUPAC-degenerate consensus (`GenerateConsensus`) |
|--------|-------------------------------|--------------------------------------------------|
| Ambiguous column | single most-frequent base | IUPAC ambiguity code (e.g. R for A/G) |
| Output alphabet | A, C, G, T | A, C, G, T + IUPAC codes |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| alignedSequences | `IEnumerable<string>` | required | Aligned DNA sequences | Equal length; alphabet {A,C,G,T}; case-insensitive |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `string` | Consensus; one uppercase base per alignment column. Empty string for an empty collection |

### 3.3 Preconditions and Validation

Null collection → `ArgumentNullException`. Empty collection → `""`. Sequences of unequal length → `ArgumentException`. Any character outside {A,C,G,T} (after uppercasing) → `ArgumentException`. Input is uppercased before processing (case-insensitive); indexing is 0-based.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input; uppercase all sequences.
2. For each column, count occurrences of A, C, G, T (the profile column) [2].
3. Select the base with the maximum count, scanning A→C→G→T so ties resolve alphabetically [4].
4. Append the selected base to the consensus.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

Alphabet/order table: `{'A','C','G','T'}` — also the tie-break order [4]. No scoring matrix or plurality threshold is applied (contrast EMBOSS `cons`, which gates output on a weighted plurality value defaulting to half the total sequence weight [3]).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Consensus | O(n × m) | O(n) | n = column count, m = sequence count; O(1) extra per column for the 4-element profile |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.CreateConsensusFromAlignment(IEnumerable<string>)`: column-wise most-frequent consensus with alphabetical tie-break.

### 5.2 Current Behavior

Sequences are uppercased via `ToUpperInvariant`. Per-column counts use a 4-element array in alphabetical order; the maximum is found with a strict `>` comparison while iterating in alphabetical order, so the first (alphabetically-earliest) maximum wins on a tie. **Search reuse:** the repository suffix tree was evaluated and is N/A — this is a column-wise tally over aligned positions, not a substring/occurrence search, so no pattern matching is involved.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Most-frequent residue per column = profile-matrix column maximum [1][2].
- Deterministic alphabetical tie-break (A<C<G<T) [4].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Weighted plurality threshold and no-consensus 'n'/'x' output from EMBOSS `cons` [3]; **users should rely on:** the IUPAC-degenerate `MotifFinder.GenerateConsensus` for ambiguity-aware consensus, or an external EMBOSS run for threshold-gated consensus.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty collection | `""` | nothing to summarise; mirrors `GenerateConsensus` |
| Single sequence | returns it (uppercased) | each column's only base is its maximum [2] |
| Tie column (A,G) | alphabetically-earliest (A) | tie-break rule [4] |
| Lowercase input | normalised to uppercase | case-insensitive contract |

### 6.2 Limitations

Operates on the DNA alphabet {A,C,G,T} only (no IUPAC ambiguity input, no gaps, no protein). Requires pre-aligned equal-length input; it does not perform alignment. No confidence/plurality threshold is applied — every column yields a base.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var aligned = new[]
{
    "ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC",
    "TTGGAACT", "ATGCCATT", "ATGGCACT"
};
string consensus = MotifFinder.CreateConsensusFromAlignment(aligned); // "ATGCAACT"
```

**Numerical walk-through:** For the Rosalind CONS sample [2], the profile is A=`5 1 0 0 5 5 0 0`, C=`0 0 1 4 2 0 6 1`, G=`1 1 6 3 0 1 0 0`, T=`1 5 0 0 0 1 1 6`. Taking the column maxima (5→A, 5→T, 6→G, 4→C, 5→A, 5→A, 6→C, 6→T) gives `ATGCAACT`.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MotifFinder_CreateConsensusFromAlignment_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_CreateConsensusFromAlignment_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [MOTIF-CONS-001-Evidence.md](../../../docs/Evidence/MOTIF-CONS-001-Evidence.md)

## 8. References

1. Wikipedia contributors. 2026. Consensus sequence. Wikipedia. https://en.wikipedia.org/wiki/Consensus_sequence
2. Rosalind. Consensus and Profile (CONS). https://rosalind.info/problems/cons/
3. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2 (program docs: https://www.bioinformatics.nl/cgi-bin/emboss/help/cons)
4. Los Alamos HIV Sequence Database. Advanced Consensus Maker — explanation. https://hfv.lanl.gov/content/sequence/CONSENSUS/AdvConExplain.html
