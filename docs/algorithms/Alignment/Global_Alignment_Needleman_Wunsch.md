# Global Alignment (Needleman-Wunsch)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-GLOBAL-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Reference |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Needleman-Wunsch computes an optimal global alignment of two sequences under a dynamic-programming scoring model. It aligns both sequences end-to-end, inserting gaps as needed so the aligned strings have equal length. In this repository, the public pairwise global-alignment entry points are `SequenceAligner.GlobalAlign(...)`, with overloads for `DnaSequence`, raw `string`, and cancellation-aware variants for both typed and raw-string inputs. The repository implementation uses a linear gap penalty derived from `ScoringMatrix.GapExtend`.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Global alignment is used when the full length of both input sequences is intended to participate in one alignment rather than searching for a highest-scoring internal region. The resulting alignment represents both sequences in a common coordinate space by inserting gap characters where necessary. This is the standard "end-to-end" formulation described in the Needleman-Wunsch and sequence-alignment references.

### 2.2 Core Model

With linear gap penalty $d$ and substitution score $S(a,b)$, the standard Needleman-Wunsch recurrence is:

$$
F_{i,j} = \max\left(F_{i-1,j-1} + S(A_i, B_j),\; F_{i-1,j} + d,\; F_{i,j-1} + d\right)
$$

with boundary conditions:

$$
F_{0,j} = d \cdot j, \qquad F_{i,0} = d \cdot i
$$

Traceback starts at $F_{m,n}$ and follows the maximizing predecessor choices to recover one optimal alignment. The total alignment score is the sum of the per-column match, mismatch, and gap contributions implied by that traceback. (Needleman-Wunsch algorithm; Sequence alignment)

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The two aligned output strings have equal length. | Global alignment inserts gaps until both sequences occupy the same alignment coordinate system. |
| INV-02 | Removing gap characters from the aligned outputs recovers the original input sequences. | Traceback emits original characters or gap symbols; it does not rewrite sequence characters. |
| INV-03 | Under a linear scoring model, the reported score equals the sum of the per-column match, mismatch, and gap scores for the returned alignment. | The dynamic program maximizes exactly that additive objective. |
| INV-04 | Traceback covers the full lengths of both inputs, from $F_{m,n}$ back to the origin. | Needleman-Wunsch is defined for full-length global alignment rather than local or semi-global matching. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence1` | `DnaSequence` or `string` | required | First sequence to align. | `DnaSequence` overload: non-null. `string` overload: returns `AlignmentResult.Empty` if null or empty. |
| `sequence2` | `DnaSequence` or `string` | required | Second sequence to align. | Same handling as `sequence1`. |
| `scoring` | `ScoringMatrix` | `SequenceAligner.SimpleDna` | Match, mismatch, and gap values used by the dynamic program. | Optional in all public overloads. |
| `cancellationToken` | `CancellationToken` | required on cancellation-aware overloads | Allows cancellation of the long-running matrix fill and traceback. | Used by both `GlobalAlign(DnaSequence, DnaSequence, ..., CancellationToken, ...)` and `GlobalAlign(string, string, ..., CancellationToken, ...)`. |
| `progress` | `IProgress<double>` | `null` | Optional progress reporter for the cancellation-aware overloads. | Values are reported during matrix fill and before/after traceback. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `AlignedSequence1` | `string` | First aligned sequence, including gap characters (`-`). |
| `AlignedSequence2` | `string` | Second aligned sequence, including gap characters (`-`). |
| `Score` | `int` | Alignment score under the chosen linear scoring matrix. |
| `AlignmentType` | `AlignmentType` | `Global` for successful global alignments; `AlignmentResult.Empty` carries `Global` as its type tag as well. |
| `StartPosition1` | `int` | Current implementation returns `0` for non-empty global alignments. |
| `StartPosition2` | `int` | Current implementation returns `0` for non-empty global alignments. |
| `EndPosition1` | `int` | Current implementation returns `sequence1.Length - 1` for non-empty global alignments. |
| `EndPosition2` | `int` | Current implementation returns `sequence2.Length - 1` for non-empty global alignments. |

### 3.3 Preconditions and Validation

`SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ...)` throws `ArgumentNullException` when either `DnaSequence` argument is null. `DnaSequence` construction normalizes input to uppercase and validates that each character is one of `A`, `C`, `G`, or `T`; invalid characters raise `ArgumentException`. `SequenceAligner.GlobalAlign(string, string, ...)` uppercases the supplied strings but does not apply `DnaSequence` validation; if either string is null or empty, it returns `AlignmentResult.Empty` instead of throwing. Empty `DnaSequence` values remain legal inputs. In the non-cancellation typed overload they are aligned through the normal traceback path, while the cancellation-aware typed overload delegates to the cancellation-aware string path and therefore returns `AlignmentResult.Empty` on empty typed inputs. Coordinates in `AlignmentResult` are 0-based.

## 4. Algorithm

### 4.1 High-Level Steps

1. Initialize the first row and first column of the score matrix with cumulative linear gap penalties.
2. Fill the remaining matrix cells with the maximum of diagonal, up, and left predecessor scores.
3. Start traceback from the bottom-right cell.
4. On each traceback step, emit a character pair or a character-gap pair according to the predecessor that realizes the current score.
5. Reverse the accumulated traceback characters and return the aligned strings with the terminal score.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The recurrence uses one substitution score per aligned character pair and one linear gap penalty per gap position. In the repository implementation, diagonal moves are preferred over up moves, and up moves are preferred over left moves when multiple predecessors yield the same score. That tie order affects which optimal alignment is returned when multiple optima exist, but it does not change the optimal score.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Global alignment of sequences of lengths `m` and `n` | `O(mn)` | `O(mn)` | Standard Needleman-Wunsch complexity under the full dynamic-programming matrix. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)`: validates non-null `DnaSequence` inputs and delegates to `GlobalAlignCore`.
- `SequenceAligner.GlobalAlign(string, string, ScoringMatrix?)`: uppercases raw strings and returns `AlignmentResult.Empty` on null or empty input.
- `SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?, CancellationToken, IProgress<double>?)`: validates non-null typed inputs and then delegates through the cancellation-aware string path.
- `SequenceAligner.GlobalAlign(string, string, ScoringMatrix?, CancellationToken, IProgress<double>?)`: cancellation-aware variant that fills a 2D score matrix and reports progress.
- `SequenceAligner.GlobalAlignCore(string, string, ScoringMatrix)`: pooled-array implementation of the core dynamic program.
- `SequenceAligner.Traceback(...)`: reconstructs the final aligned strings from the score matrix.

**Supporting types:** [AlignmentTypes.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

### 5.2 Current Behavior

The repository exposes four global-alignment entry points, but all compute the same linear-gap algorithm. `GlobalAlignCore` rents a flat integer buffer from `ArrayPool<int>`, fills it as a flattened matrix, copies it into an `int[,]`, and then calls `Traceback`. The cancellation-aware overload uses an `int[,]` directly and checks cancellation periodically during matrix fill and traceback. The traceback code is deterministic because it tests diagonal first, then up, then left. `ScoringMatrix.GapOpen` is not used by the global dynamic program; the gap cost comes from `ScoringMatrix.GapExtend`. Only the non-cancellation typed overload sends empty `DnaSequence` inputs through `Traceback`, producing end coordinates derived as `sequence.Length - 1` and therefore `-1` for empty typed inputs. The cancellation-aware typed overload delegates to the raw-string cancellation path and returns `AlignmentResult.Empty` on empty typed inputs.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Linear-gap Needleman-Wunsch boundary initialization `F(i,0)=d*i` and `F(0,j)=d*j`.
- The standard three-way recurrence over diagonal, up, and left predecessors.
- Traceback from the bottom-right cell to recover one optimal global alignment.

**Intentionally simplified:**

- The public `ScoringMatrix` includes `GapOpen`, but `GlobalAlign` uses only `GapExtend`; **consequence:** the repository implements a linear gap model rather than an affine gap model.

**Not implemented:**

- (none)

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty raw string input | Returns `AlignmentResult.Empty`. | The string overload exits early when either input is null or empty. |
| Empty `DnaSequence` input on the non-cancellation typed overload | Returns a typed `Global` result whose end coordinates are `-1`. | That overload delegates to the core routine and `Traceback` sets end coordinates to `sequence.Length - 1`. |
| Empty `DnaSequence` input on the cancellation-aware typed overload | Returns `AlignmentResult.Empty`. | That overload delegates to the cancellation-aware raw-string path, which short-circuits on empty input. |
| Null `DnaSequence` input | Throws `ArgumentNullException`. | Public `DnaSequence` overloads guard with `ThrowIfNull`. |
| Identical sequences | Returns a full-length gap-free alignment with score equal to the sum of match scores. | Verified by `SequenceAligner_GlobalAlign_Tests` and the Needleman-Wunsch recurrence. |
| Unequal-length sequences | Uses initialized border cells to absorb leading or trailing gaps in the optimal alignment. | This is the purpose of the boundary conditions `F(i,0)` and `F(0,j)`. |
| Completely different equal-length sequences under negative mismatch penalties | Can return a negative score. | Global alignment has no zero floor; every position still participates in the alignment. |

### 6.2 Limitations

The implementation returns one optimal alignment, not all optimal alignments, because traceback uses a fixed tie order. The public string overload does not validate the alphabet beyond uppercasing the input. The full dynamic-programming matrix keeps the asymptotic memory requirement at `O(mn)`. Affine gap penalties are not implemented in the pairwise global-alignment API.

## 7. Examples and Related Material

### 7.1 Worked Example

With `Match = +1`, `Mismatch = -1`, and `GapExtend = -1`, the repository tests validate the standard `GCATGCG` versus `GATTACA` example with optimal score `0`. The exact aligned strings are allowed to vary among optimal solutions, but the tests require the global-alignment invariants and the final score described by the cited Needleman-Wunsch example.

## 8. References

1. [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
2. [Sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment)
