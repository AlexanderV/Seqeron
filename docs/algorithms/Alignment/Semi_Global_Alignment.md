# Semi-Global Alignment (Ends-Free / Glocal)

## Documented Theory (Authoritative Sources)

### Purpose
Semi-global alignment (also called "glocal" or "ends-free" alignment) is a **hybrid** of global and local alignment that searches for the best possible partial alignment of two sequences. It permits alignment to start and/or end at any position in one or both sequences, avoiding the penalties for unaligned ends. (Wikipedia: Sequence alignment)

### Definition
From Wikipedia (Sequence alignment):
> "Hybrid methods, known as semi-global or 'glocal' (short for global-local) methods, search for the best possible partial alignment of the two sequences (in other words, a combination of one or both starts and one or both ends is stated to be aligned)."

### Primary Use Cases

1. **Overlap alignment**: When the downstream part of one sequence overlaps with the upstream part of another. Neither global nor local alignment is ideal:
   - Global alignment forces extension beyond the overlap region
   - Local alignment might not fully cover the overlap region

2. **Short read alignment (query-in-reference)**: When aligning a short sequence (e.g., a gene, primer, or sequencing read) to a much longer sequence (e.g., a chromosome or reference genome). The short sequence should be **globally aligned** (fully matched) while only a local portion of the long sequence participates.

(Wikipedia: Sequence alignment)

### Algorithm Variants

Multiple configurations exist depending on which sequence ends are "free" (no gap penalty):

| Variant | Free End Gaps | Typical Use Case |
|---------|---------------|------------------|
| Query-in-reference | Start and end of reference | Short read mapping, primer alignment |
| Overlap | End of seq1, start of seq2 | Sequence assembly |
| Full semi-global | All ends free | General substring finding |

The query-in-reference variant is most common in bioinformatics for mapping short reads to references.

### Relationship to Needleman–Wunsch

Semi-global alignment is a modification of the Needleman–Wunsch (global) algorithm with these differences:

| Aspect | Global (NW) | Semi-Global |
|--------|-------------|-------------|
| First row initialization | Gap penalties accumulate | Zero (free start gaps in reference) |
| First column initialization | Gap penalties accumulate | Gap penalties (query fully aligned) |
| Traceback start | Bottom-right cell | Maximum in last row |
| End-gap handling | Penalized | Free in designated sequence |

(Derived from Wikipedia: Needleman–Wunsch algorithm)

### Scoring Model

The scoring uses the same recurrence as Needleman–Wunsch:

$$F_{i,j} = \max\left(F_{i-1,j-1} + s(a_i, b_j),\; F_{i-1,j} + d,\; F_{i,j-1} + d\right)$$

Where:
- $s(a,b)$ is the substitution score (match/mismatch)
- $d$ is the gap penalty (typically negative)

The key difference is in **initialization** and **traceback**:
- First row: $F_{0,j} = 0$ for all $j$ (free leading gaps in reference)
- First column: $F_{i,0} = i \cdot d$ (query must be fully aligned)
- Traceback: Start from $\max_j(F_{m,j})$ where $m$ = query length

### Complexity

For sequences of lengths $n$ (query) and $m$ (reference):
- **Time**: $O(n \cdot m)$
- **Space**: $O(n \cdot m)$ for scoring matrix

(Wikipedia: Needleman–Wunsch algorithm)

---

## Implementation Notes (Seqeron.Genomics)

**Implementation Location:** [Seqeron.Genomics/SequenceAligner.cs](../../../src/Seqeron/Seqeron.Genomics/SequenceAligner.cs)

### Methods
- `SequenceAligner.SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` – canonical method
- `SemiGlobalAlignCore(string, string, ScoringMatrix)` – core implementation

### Implementation Details (Query-in-Reference Variant)

1. **Matrix Initialization**:
   - First row: All zeros (free start gaps in seq2/reference)
   - First column: Accumulating gap penalties (seq1/query fully aligned)

2. **Traceback Start**:
   - Finds maximum score in last row (not bottom-right corner)
   - This allows query to align anywhere within reference

3. **Trailing Gap Handling**:
   - Appends remaining reference bases as gaps in aligned seq1
   - No penalty for unaligned trailing reference

4. **AlignmentResult Properties**:
   - `AlignmentType`: Set to `SemiGlobal`
   - `AlignedSequence1`, `AlignedSequence2`: Aligned sequences with gaps
   - `Score`: Alignment score (can be negative unlike local alignment)

### Key Invariants

| ID | Invariant | Description |
|----|-----------|-------------|
| INV-1 | AlignmentType = SemiGlobal | Type marker set correctly |
| INV-2 | Query fully represented | Removing gaps from aligned seq1 = original seq1 |
| INV-3 | Equal alignment length | len(aligned1) = len(aligned2) |
| INV-4 | Reference substring | Removing gaps from aligned seq2 = substring of original seq2 |

---

## Deviations / Assumptions

1. **ASSUMPTION**: The implementation uses the query-in-reference variant (free end gaps in seq2 only). Other variants (overlap, full semi-global) are not currently implemented.

2. **ASSUMPTION**: Null sequence inputs throw `ArgumentNullException`. Empty sequences may not be explicitly handled.

3. **ASSUMPTION**: The score can be negative (unlike Smith–Waterman) because there is no zero floor.

---

## References

1. Wikipedia. "Sequence alignment." https://en.wikipedia.org/wiki/Sequence_alignment
2. Wikipedia. "Needleman–Wunsch algorithm." https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
3. Brudno M et al. (2003). "Glocal alignment: finding rearrangements during alignment." Bioinformatics 19 Suppl 1:i54-62.
