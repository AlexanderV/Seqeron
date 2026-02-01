# Multiple Sequence Alignment (MSA)

## Overview

Multiple Sequence Alignment (MSA) is the process of aligning three or more biological sequences (DNA, RNA, or protein) by inserting gaps to maximize similarity across all sequences simultaneously.

**Implementation:** `SequenceAligner.MultipleAlign()`  
**Complexity:** O(n² × m) where n = number of sequences, m = average sequence length  
**Algorithm Variant:** Star Alignment (simplified progressive alignment)

---

## Theory

### Problem Definition

Given m sequences S₁, S₂, ..., Sₘ of varying lengths, an MSA transforms them into sequences S'₁, S'₂, ..., S'ₘ of equal length L by inserting gap characters (`-`), subject to:
- L ≥ max{|Sᵢ|} for all i
- No column consists entirely of gaps
- Removing gaps from S'ᵢ recovers Sᵢ

### Computational Complexity

Finding the optimal MSA is **NP-complete** (Wang & Jiang, 1994). For n sequences of length L, exact dynamic programming requires O(Lⁿ) time and space, making it infeasible for more than a few sequences.

### Progressive Alignment

The most widely used heuristic approach, developed by Feng & Doolittle (1987):

1. **Pairwise alignment phase:** Compute distances between all sequence pairs
2. **Guide tree construction:** Build tree (e.g., UPGMA, neighbor-joining) from distances
3. **Progressive alignment phase:** Align sequences/profiles following guide tree order

### Star Alignment

A simplified variant of progressive alignment:
- Select one sequence as the "center" (reference)
- Align all other sequences to the center using pairwise alignment
- Combine alignments into a multiple alignment

**Advantages:** Simple, fast, O(n) pairwise alignments  
**Disadvantages:** Result depends heavily on center selection; errors in center alignment propagate

---

## Implementation Details

### Algorithm Steps

```
1. Input validation: null check, empty check
2. Trivial case: single sequence returns unchanged
3. Star alignment:
   a. Use first sequence as reference (center)
   b. For each other sequence:
      - Perform global alignment with reference
      - Store aligned sequence
   c. Pad all sequences to equal length
4. Consensus generation:
   a. For each column position:
      - Count occurrences of A, C, G, T (excluding gaps)
      - Select most frequent nucleotide
   b. Return majority-voted consensus
5. Return MultipleAlignmentResult
```

### Data Structures

```csharp
public sealed record MultipleAlignmentResult(
    string[] AlignedSequences,  // All sequences aligned to equal length
    string Consensus,           // Majority-voted consensus sequence  
    int TotalScore)             // Sum of pairwise alignment scores
{
    public static MultipleAlignmentResult Empty => 
        new(Array.Empty<string>(), "", 0);
}
```

### Scoring

The implementation uses **sum-of-pairs scoring**:
- TotalScore = Σ Score(reference, sequence_i) for all i > 0

---

## Invariants

| Invariant | Description |
|-----------|-------------|
| Equal length | All aligned sequences have the same length |
| Count preservation | Output contains same number of sequences as input |
| Sequence recovery | Removing gaps from aligned sequence recovers original |
| Valid characters | Aligned sequences contain only {A, C, G, T, -} |
| Consensus validity | Consensus characters from {A, C, G, T, -} |

---

## Edge Cases

| Input | Behavior |
|-------|----------|
| `null` | Throws `ArgumentNullException` |
| Empty collection | Returns `MultipleAlignmentResult.Empty` |
| Single sequence | Returns alignment with that sequence; consensus = sequence |
| Two sequences | Equivalent to global alignment between them |
| Identical sequences | Perfect alignment; consensus = any input sequence |
| Different lengths | Shorter sequences padded with gaps |

---

## Limitations

1. **Order dependence:** Using first sequence as reference; different orderings may produce different results
2. **No guide tree:** Does not construct phylogenetic guide tree like ClustalW
3. **Error propagation:** Errors in early alignments propagate to final result
4. **Suboptimal for divergent sequences:** Star alignment is less accurate than sophisticated methods (T-Coffee, MUSCLE) for evolutionarily distant sequences

---

## References

1. Wikipedia. "Multiple sequence alignment." https://en.wikipedia.org/wiki/Multiple_sequence_alignment
2. Wikipedia. "Clustal." https://en.wikipedia.org/wiki/Clustal
3. Feng DF, Doolittle RF (1987). "Progressive sequence alignment." J Mol Evol. 25(4):351-360.
4. Wang L, Jiang T (1994). "On the complexity of multiple sequence alignment." J Comput Biol. 1(4):337-348.
5. Thompson JD, Higgins DG, Gibson TJ (1994). "CLUSTAL W." Nucleic Acids Res. 22(22):4673-80.
