# Semi-Global Alignment (Fitting / Query-in-Reference)

## Documented Theory (Authoritative Sources)

### Purpose
Semi-global alignment (also called "glocal" or "ends-free" alignment) is a **hybrid** of global and local alignment that searches for the best possible partial alignment of two sequences. It permits alignment to start and/or end at any position in one or both sequences, avoiding the penalties for unaligned ends.

The implementation uses the **fitting alignment** variant (also called "query-in-reference"), which aligns the entire query globally against the best-matching substring of the reference.

### Definition

From Wikipedia (Sequence alignment):
> "Hybrid methods, known as semi-global or 'glocal' (short for global-local) methods, search for the best possible partial alignment of the two sequences (in other words, a combination of one or both starts and one or both ends is stated to be aligned)."

> "Another case where semi-global alignment is useful is when one sequence is short (for example a gene sequence) and the other is very long (for example a chromosome sequence). In that case, the short sequence should be globally (fully) aligned but only a local (partial) alignment is desired for the long sequence."

From Rosalind (SIMS — Finding a Motif with Modifications):
> A fitting alignment is "an alignment of a substring of s against all of t." The aim is to find a substring s′ of s that maximizes an alignment score with respect to t.

From Rosalind (SMGB — Semiglobal Alignment glossary):
> "A semiglobal alignment of strings s and t is an alignment in which any gaps appearing as prefixes or suffixes of s and t do not contribute to the score of the alignment."

Note: The implementation specifically implements the **fitting alignment** (query-in-reference) variant, which is a subset of the broader glocal family.

### Algorithm Variants

| Variant | Free End Gaps | Typical Use Case | Rosalind Problem |
|---------|---------------|------------------|------------------|
| Query-in-reference (fitting) | Start and end of reference | Short read mapping, primer alignment | SIMS |
| Overlap | End of seq1, start of seq2 | Sequence assembly | OAP |
| Full semi-global | All four ends free | General substring finding | SMGB |

### Relationship to Needleman–Wunsch

The fitting alignment is a modification of the Needleman–Wunsch (global) algorithm with modified initialization and traceback:

| Aspect | Global (NW) | Fitting (Query-in-Reference) |
|--------|-------------|------------------------------|
| First row $F_{0,j}$ | $d \cdot j$ (gap penalties) | $0$ (free start gaps in reference) |
| First column $F_{i,0}$ | $d \cdot i$ (gap penalties) | $d \cdot i$ (query fully aligned) |
| Recurrence | $\max(\text{diag, up, left})$ | Same (no zero floor) |
| Traceback start | $F_{m,n}$ (bottom-right) | $\max_j F_{m,j}$ (max in last row) |
| End-gap handling | All penalized | Free in reference |

Sources: Wikipedia (Needleman–Wunsch algorithm), Rosalind (SIMS, SMGB).

### Scoring Model

The scoring uses the same recurrence as Needleman–Wunsch with linear gap cost:

$$F_{i,j} = \max\left(F_{i-1,j-1} + s(a_i, b_j),\; F_{i-1,j} + d,\; F_{i,j-1} + d\right)$$

Where:
- $s(a,b)$ is the substitution score (match/mismatch)
- $d$ is the linear gap penalty (typically negative)

The key difference from Smith–Waterman (local alignment) is the absence of a zero floor — the score can be negative.

Source: Wikipedia (Needleman–Wunsch algorithm), "Advanced presentation of algorithm" section.

### Complexity

For query length $m$ and reference length $n$:
- **Time**: $O(m \cdot n)$
- **Space**: $O(m \cdot n)$ for the scoring matrix

Source: Wikipedia (Needleman–Wunsch algorithm), "Complexity" section.

---

## Implementation Notes (Seqeron.Genomics)

**Implementation Location:** [Seqeron.Genomics.Alignment/SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

### Methods
- `SequenceAligner.SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` — public entry point
- `SemiGlobalAlignCore(string, string, ScoringMatrix)` — core DP implementation

### Implementation Details

1. **Matrix Initialization**:
   - First row: All zeros (free start gaps in reference)
   - First column: Accumulating gap penalties via `GapExtend` (linear cost; query fully aligned)

2. **Recurrence**: Standard NW with linear gap cost (`GapExtend` only; `GapOpen` is unused).

3. **Traceback Start**: Maximum score in last row $\max_j F_{m,j}$ (not bottom-right cell).

4. **Trailing Gap Handling**: Remaining reference bases after the alignment endpoint are appended as gaps in aligned seq1 (no penalty).

5. **Score**: Retrieved from $F_{m, \text{maxJ}}$ — the optimal fitting score.

### Key Invariants

| ID | Invariant | Description |
|----|-----------|-------------|
| INV-1 | AlignmentType = SemiGlobal | Type marker set correctly |
| INV-2 | Equal alignment length | `len(aligned1) == len(aligned2)` |
| INV-3 | Query fully represented | `RemoveGaps(aligned1) == query` |
| INV-4 | Reference is substring | `RemoveGaps(aligned2)` is a substring of original reference |
| INV-5 | Score = $\max_j F_{m,j}$ | Score equals the maximum in the last row of the DP matrix |

---

## References

1. Wikipedia. "Sequence alignment." https://en.wikipedia.org/wiki/Sequence_alignment
2. Wikipedia. "Needleman–Wunsch algorithm." https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
3. Rosalind. "Finding a Motif with Modifications (SIMS)." https://rosalind.info/problems/sims/
4. Rosalind. "Semiglobal Alignment (SMGB)." https://rosalind.info/problems/smgb/
5. Brudno M et al. (2003). "Glocal alignment: finding rearrangements during alignment." Bioinformatics 19 Suppl 1:i54-62.
