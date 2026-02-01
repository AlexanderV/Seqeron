# Global Alignment (Needleman–Wunsch)

## Documented Theory (Authoritative Sources)

### Purpose
Needleman–Wunsch is a dynamic programming algorithm that computes an **optimal global alignment** between two sequences, aligning them end-to-end by maximizing an alignment score. (Needleman–Wunsch algorithm; Sequence alignment)

### Scoring Model
An alignment score is computed as the sum of per-position scores for:
- **Match** (same character)
- **Mismatch** (different character)
- **Indel/Gap** (character aligned to a gap)

The standard recurrence uses a linear gap penalty $d$ and a similarity function $S(a,b)$:

$$F_{i,j} = \max\left(F_{i-1,j-1} + S(A_i,B_j),\; F_{i-1,j} + d,\; F_{i,j-1} + d\right)$$

with boundary initialization:

$$F_{0,j} = d \cdot j, \quad F_{i,0} = d \cdot i$$

Traceback follows the max choices to build the aligned sequences, inserting gaps on vertical/horizontal moves. (Needleman–Wunsch algorithm)

### Global Alignment Properties
Global alignment attempts to align the **entire lengths** of both sequences, inserting gaps so the aligned sequences have equal length. (Sequence alignment)

### Complexity
For sequences of lengths $n$ and $m$, Needleman–Wunsch uses $O(nm)$ time and $O(nm)$ space. (Needleman–Wunsch algorithm)

---

## Implementation Notes (Seqeron.Genomics)

**Implementation Location:** [Seqeron.Genomics/SequenceAligner.cs](Seqeron.Genomics/SequenceAligner.cs)

- `SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` implements global alignment via dynamic programming.
- The recurrence uses **linear gap penalties** via `GapExtend` for insertions/deletions.
- Initialization includes `GapOpen` added once to the first row/column for non-zero lengths.
- `GlobalAlign(string, string, ScoringMatrix?)` normalizes inputs to uppercase and returns `AlignmentResult.Empty` when either input is null or empty.
- `DnaSequence` inputs are validated to contain only A/C/G/T and normalized to uppercase.

---

## Deviations / Assumptions

1. **Affine gaps:** The implementation uses `GapOpen` only in initialization and `GapExtend` in the recurrence, which is not a full affine gap model (Gotoh). This is an implementation-specific choice.
2. **Empty-input behavior:** Returning `AlignmentResult.Empty` for empty string inputs is an implementation behavior and not specified in the cited sources.
3. **Multiple optimal alignments:** The implementation returns a single traceback path; when multiple optimal paths exist, the chosen alignment is implementation-dependent.

---

## Sources

- https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
- https://en.wikipedia.org/wiki/Sequence_alignment
