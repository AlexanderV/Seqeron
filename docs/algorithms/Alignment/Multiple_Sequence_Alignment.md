# Multiple Sequence Alignment (MSA)

## Overview

Multiple Sequence Alignment (MSA) is the process of aligning three or more biological sequences (DNA, RNA, or protein) by inserting gaps to maximize similarity across all sequences simultaneously.

**Implementation:** `SequenceAligner.MultipleAlign()`
**Complexity:** O(k² × m) where k = number of sequences, m = average sequence length
**Algorithm Variant:** Anchor-based Star Alignment (progressive alignment with suffix tree anchors)

---

## Theory

### Problem Definition

Given m sequences S₁, S₂, ..., Sₘ of varying lengths, an MSA transforms them into sequences S'₁, S'₂, ..., S'ₘ of equal length L by inserting gap characters (`-`), subject to:
- L ≥ max{|Sᵢ|} for all i
- No column consists entirely of gaps (Wikipedia MSA)
- Removing gaps from S'ᵢ recovers Sᵢ (Wikipedia MSA)

### Computational Complexity

Finding the optimal MSA is **NP-complete** (Wang & Jiang, 1994). For n sequences of length L, exact dynamic programming requires O(Lⁿ) time and space, making it infeasible for more than a few sequences.

### Progressive Alignment

The most widely used heuristic approach, developed by Feng & Doolittle (1987):

1. **Pairwise alignment phase:** Compute distances between all sequence pairs
2. **Guide tree construction:** Build tree (e.g., UPGMA, neighbor-joining) from distances
3. **Progressive alignment phase:** Align sequences/profiles following guide tree order

### Star Alignment

A simplified variant of progressive alignment:
- Select one sequence as the "center" (the sequence with highest total similarity to all others)
- Align all other sequences to the center using pairwise alignment
- Combine alignments into a multiple alignment

*Note: Star alignment is a well-known algorithmic concept in the MSA literature, but is NOT described on the Wikipedia Clustal page. ClustalW uses guide-tree-based progressive alignment, not center-star.*

**Advantages:** Simple, fast, O(k) pairwise alignments
**Disadvantages:** Result depends heavily on center selection; errors in center alignment propagate

---

## Implementation Details

### Algorithm Steps

```
1. Input validation: null check, empty check
2. Trivial case: single sequence returns unchanged
3. Anchor-based star alignment:
   a. Select center sequence via k-mer (4-mer) cosine similarity
      - For each pair (i,j): compute k-mer frequency vectors, cosine similarity
      - Center = sequence with highest total similarity to all others
      - Note: ClustalV uses k-tuple match counts (Wikipedia Clustal);
        cosine similarity on frequency vectors is an implementation design choice
   b. Build suffix tree on center sequence (implementation-specific)
   c. For each other sequence (parallel):
      - Find exact-match anchors via suffix tree
      - Align gaps between anchors using Needleman-Wunsch
      - Store pairwise alignment result
   d. Gap reconciliation:
      - Merge gap columns from independent pairwise alignments
      - Project all sequences into merged coordinate space
   e. Pad all sequences to equal length
4. Consensus generation (Wikipedia Consensus sequence: "most frequent residues... at each position"):
   a. For each column position:
      - Count occurrences of A, C, G, T, and '-'
      - Select character with highest count (gaps participate — implementation choice)
      - On tie between gap and nucleotide, prefer nucleotide (implementation choice)
   b. Return majority-voted consensus
5. Sum-of-pairs scoring (Wikipedia MSA: "sum of all of the pairs of characters at each position"):
   a. For each column, for each pair (i < j):
      - Both match: +Match
      - Both differ: +Mismatch
      - Gap-nucleotide: +GapExtend
      - Gap-gap: 0 (standard bioinformatics convention, not stated in Wikipedia)
   b. TotalScore = sum across all C(k,2) pairs and all columns
6. Return MultipleAlignmentResult
```

### Data Structures

```csharp
public sealed record MultipleAlignmentResult(
    string[] AlignedSequences,  // All sequences aligned to equal length
    string Consensus,           // Majority-voted consensus sequence (gaps included)
    int TotalScore)             // Sum-of-pairs score across all C(k,2) pairs
{
    public static MultipleAlignmentResult Empty =>
        new(Array.Empty<string>(), "", 0);
}
```

### Scoring

The implementation uses **column-based sum-of-pairs (SP) scoring** (Wikipedia MSA):
- TotalScore = Σ_{columns} Σ_{i<j} score(S'ᵢ[col], S'ⱼ[col])
- Sums scores across ALL C(k,2) sequence pairs, not just center-to-other
- Gap-gap pairs contribute 0 (neutral) — standard bioinformatics convention, not stated in Wikipedia
- Gap-nucleotide pairs contribute GapExtend penalty

---

## Invariants

| Invariant | Description | Source |
|-----------|-------------|--------|
| Equal length | All aligned sequences have the same length | Wikipedia MSA |
| No all-gap columns | No column consists entirely of gaps | Wikipedia MSA |
| Count preservation | Output contains same number of sequences as input | Wikipedia MSA |
| Sequence recovery | Removing gaps from aligned sequence recovers original | Wikipedia MSA |
| Valid characters | Aligned sequences contain only {A, C, G, T, -} | Wikipedia MSA |
| Consensus validity | Consensus characters from {A, C, G, T, -} | Wikipedia MSA |
| SP score | TotalScore equals column-based SP across all C(k,2) pairs | Wikipedia MSA |

---

## Implementation Design Choices (not from external sources)

The following decisions are implementation-specific and NOT derived from external sources:

| Choice | Rationale |
|--------|----------|
| 4-mer cosine similarity for center selection | ClustalV uses k-tuple match counts; cosine on frequency vectors is our variant |
| Suffix tree anchors for pairwise alignment | Implementation optimization; not from any MSA reference |
| Gaps participate in consensus vote | Reasonable default; Wikipedia Consensus sequence says "most frequent residues" without specifying gap handling |
| Nucleotide preferred over gap on tie | Conservative choice: prefer informative characters over gap |
| Gap-gap = 0 in SP scoring | Standard bioinformatics convention; not explicitly stated in Wikipedia |
| Star alignment (center-based) | Well-known MSA simplification; not described on Wikipedia Clustal page |

---

## Edge Cases

| Input | Behavior |
|-------|----------|
| `null` | Throws `ArgumentNullException` |
| Empty collection | Returns `MultipleAlignmentResult.Empty` |
| Single sequence | Returns alignment with that sequence; consensus = sequence; TotalScore = 0 |
| Two sequences | Equivalent to pairwise alignment; SP = one pair score |
| Identical sequences | Perfect alignment; consensus = any input sequence |
| Different lengths | Shorter sequences padded with gaps |

---

## Limitations

1. **Center dependence:** Result depends on which sequence is selected as center (star alignment property)
2. **No guide tree:** Simplified star alignment; does not construct phylogenetic guide tree (ClustalW uses UPGMA/neighbor-joining guide trees per Wikipedia Clustal)
3. **Error propagation:** Errors in early alignments propagate to final result (Wikipedia MSA)
4. **Suboptimal for divergent sequences:** Star alignment is less accurate than sophisticated methods (T-Coffee, MUSCLE) for evolutionarily distant sequences

---

## References

1. Wikipedia. "Multiple sequence alignment." https://en.wikipedia.org/wiki/Multiple_sequence_alignment
2. Wikipedia. "Clustal." https://en.wikipedia.org/wiki/Clustal
3. Wikipedia. "Consensus sequence." https://en.wikipedia.org/wiki/Consensus_sequence
4. Feng DF, Doolittle RF (1987). "Progressive sequence alignment." J Mol Evol. 25(4):351-360.
5. Wang L, Jiang T (1994). "On the complexity of multiple sequence alignment." J Comput Biol. 1(4):337-348.
6. Thompson JD, Higgins DG, Gibson TJ (1994). "CLUSTAL W." Nucleic Acids Res. 22(22):4673-80.
