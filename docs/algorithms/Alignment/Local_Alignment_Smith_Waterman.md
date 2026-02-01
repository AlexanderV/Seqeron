# Local Alignment (Smith–Waterman)

## Documented Theory (Authoritative Sources)

### Purpose
Smith–Waterman is a dynamic programming algorithm that computes an **optimal local alignment** between two sequences, finding the highest-scoring pair of subsequences. Unlike global alignment (Needleman–Wunsch), it identifies regions of similarity within longer sequences that may be dissimilar overall. (Smith–Waterman algorithm; Sequence alignment)

### Algorithm Definition
The algorithm was proposed by Temple F. Smith and Michael S. Waterman in 1981. It is a variation of Needleman–Wunsch that performs local alignment by:
1. Setting negative scores to zero (floor at 0)
2. Initializing first row and column to zero (no end-gap penalty)
3. Starting traceback from the maximum score cell
4. Ending traceback when a zero is encountered

(Smith–Waterman algorithm)

### Scoring Model
An alignment score is computed using a substitution matrix $s(a,b)$ and gap penalty $W_k$. The recurrence relation is:

$$H_{i,j} = \max\left(0,\; H_{i-1,j-1} + s(a_i,b_j),\; \max_{k\geq 1}\{H_{i-k,j} - W_k\},\; \max_{l\geq 1}\{H_{i,j-l} - W_l\}\right)$$

The key difference from Needleman–Wunsch is the **zero floor**: negative scores are set to 0, indicating no similarity up to that point and allowing alignment to restart. (Smith–Waterman algorithm)

### Initialization
Both the first row and first column are initialized to zero:

$$H_{k,0} = H_{0,l} = 0 \quad \text{for } 0 \leq k \leq n \text{ and } 0 \leq l \leq m$$

This enables alignment of any segment without end-gap penalty. (Smith–Waterman algorithm)

### Linear Gap Penalty Simplification
When linear gap penalty is used ($W_k = k \cdot W_1$), the recurrence simplifies to:

$$H_{i,j} = \max\left(0,\; H_{i-1,j-1} + s(a_i,b_j),\; H_{i-1,j} - W_1,\; H_{i,j-1} - W_1\right)$$

(Smith–Waterman algorithm)

### Traceback
Traceback begins at the cell with the **highest score** in the entire matrix $H$ and proceeds backwards following the path of maximum scores until a cell with score zero is encountered. (Smith–Waterman algorithm)

### Key Properties
1. **Score ≥ 0**: All matrix cells have non-negative values due to the zero floor
2. **Optimal local subsequence**: The algorithm finds the highest-scoring pair of subsequences
3. **No end-gap penalty**: The algorithm can align any internal region without penalizing unaligned ends
4. **Guaranteed optimality**: Like Needleman–Wunsch, Smith–Waterman is guaranteed to find the optimal local alignment with respect to the scoring system used

(Smith–Waterman algorithm)

### Comparison with Needleman–Wunsch

| Aspect | Smith–Waterman | Needleman–Wunsch |
|--------|----------------|------------------|
| Initialization | First row/column = 0 | Subject to gap penalty |
| Scoring | Negative scores → 0 | Scores can be negative |
| Traceback | Start at max score, end at 0 | Start at bottom-right, end at top-left |
| Alignment type | Local (best subsequences) | Global (full sequences) |

(Smith–Waterman algorithm)

### Complexity
For sequences of lengths $n$ and $m$:
- Time: $O(nm)$ with linear gap penalty (Gotoh optimization)
- Space: $O(nm)$ standard, reducible to $O(n)$ (Myers and Miller, 1988)

(Smith–Waterman algorithm)

### Example from Wikipedia
Sequences: `TGTTACGG` and `GGTTGACTA`
Scoring: Match +3, Mismatch -3, Gap penalty $W_k = 2k$ (linear, $W_1 = 2$)

The algorithm finds the optimal local alignment:
```
G T T - A C
G T T G A C
```

(Smith–Waterman algorithm)

---

## Implementation Notes (Seqeron.Genomics)

**Implementation Location:** [Seqeron.Genomics/SequenceAligner.cs](../../../Seqeron.Genomics/SequenceAligner.cs)

### Methods
- `SequenceAligner.LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` – canonical method for local alignment
- `SequenceAligner.LocalAlign(string, string, ScoringMatrix?)` – string wrapper with null/empty handling
- `LocalAlignCore(string, string, ScoringMatrix)` – core Smith–Waterman implementation

### Implementation Details
1. **Zero floor**: Implemented via `Math.Max(0, score)` when filling matrix cells
2. **Linear gap penalty**: Uses `scoring.GapExtend` for gaps
3. **Traceback from maximum**: Tracks `maxScore` and `maxI`, `maxJ` during matrix filling
4. **Returns positions**: `AlignmentResult` includes `StartPosition1`, `EndPosition1`, `StartPosition2`, `EndPosition2`
5. **Empty input handling**: Returns `AlignmentResult.Empty` for null or empty inputs

### AlignmentResult Properties
- `AlignedSequence1`, `AlignedSequence2`: The aligned subsequences with gaps inserted
- `Score`: The alignment score (always ≥ 0 for local alignment)
- `AlignmentType`: Set to `Local`
- `StartPosition1/2`, `EndPosition1/2`: Positions of aligned regions in original sequences

---

## Deviations / Assumptions

1. **Linear gap model only**: The implementation uses linear gap penalty (`GapExtend`), not the original arbitrary gap penalty or affine gaps (Gotoh). This is an implementation choice affecting performance with biological data containing long gaps.

2. **Single optimal alignment**: The implementation returns a single traceback path. When multiple optimal paths exist, the chosen alignment is implementation-dependent.

3. **Empty-input behavior**: Returning `AlignmentResult.Empty` for null/empty inputs is an implementation behavior not specified in the cited algorithm sources.

4. **Scoring matrix defaults**: Default scoring values are implementation-specific and may differ from standard biological matrices (e.g., BLOSUM, PAM).

---

## Sources

- https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
- https://en.wikipedia.org/wiki/Sequence_alignment
- Smith, T.F. & Waterman, M.S. (1981). "Identification of Common Molecular Subsequences". Journal of Molecular Biology. 147 (1): 195–197.
