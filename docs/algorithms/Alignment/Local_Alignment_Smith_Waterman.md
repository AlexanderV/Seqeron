# Local Alignment (Smith-Waterman)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-LOCAL-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Smith-Waterman computes an optimal local alignment, meaning it returns the highest-scoring pair of aligned subsequences rather than forcing both full sequences into one end-to-end alignment. The distinguishing feature is a zero floor: negative-running scores are reset to zero so alignment can restart when similarity is lost. In this repository, `SequenceAligner.LocalAlign(...)` implements the dynamic program with a linear gap penalty and returns 0-based start and end coordinates for the aligned subsequences. Compared with the global Needleman-Wunsch method, this is the repository's exact local-alignment entry point for pairwise DNA-style inputs.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Local alignment is used when only a region of similarity is expected to be shared between the two inputs. This is appropriate when larger prefixes or suffixes may be unrelated, while one internal region still needs an optimal alignment under the chosen scoring model. Smith-Waterman is the classical dynamic-programming formulation for that task.

### 2.2 Core Model

The general Smith-Waterman formulation uses substitution score $s(a,b)$ together with a gap-cost function $W_k$:

$$
H_{i,j} = \max\left(0,\; H_{i-1,j-1} + s(a_i,b_j),\; \max_{k \ge 1}\{H_{i-k,j} - W_k\},\; \max_{l \ge 1}\{H_{i,j-l} - W_l\}\right)
$$

The defining property is the zero floor: negative values are replaced by `0`, which allows the alignment to restart after a low-scoring region. The first row and first column are initialized to zero, traceback begins from the highest-valued cell in the matrix, and traceback stops when a zero-valued cell is reached. (Smith-Waterman algorithm)

When a linear gap penalty is used, the recurrence simplifies to:

$$
H_{i,j} = \max\left(0,\; H_{i-1,j-1} + s(a_i,b_j),\; H_{i-1,j} - W_1,\; H_{i,j-1} - W_1\right)
$$

This simplified linear-gap form is the variant implemented in the repository. (Smith-Waterman algorithm)

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every matrix cell, and therefore every reported local-alignment score, is non-negative. | Smith-Waterman clamps each cell with the zero floor. |
| INV-02 | The aligned output strings have equal length. | Traceback emits one alignment column per step, inserting gaps where required. |
| INV-03 | Removing gaps from the aligned outputs yields substrings of the original inputs. | Local alignment returns subsequences, not arbitrary reordered characters. |
| INV-04 | Traceback starts at the matrix maximum and ends at a zero-valued cell. | That is the defining Smith-Waterman traceback rule. |

### 2.5 Comparison with Related Methods

| Aspect | Smith-Waterman | Needleman-Wunsch |
|--------|----------------|------------------|
| Boundary initialization | First row and column are `0` | First row and column accumulate gap penalties |
| Score floor | Negative scores are reset to `0` | Scores can remain negative |
| Traceback start | Highest-scoring cell in the matrix | Bottom-right cell |
| Alignment scope | Best local subsequences | Full sequence lengths |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence1` | `DnaSequence` or `string` | required | First sequence to align locally. | `DnaSequence` overload: non-null. `string` overload: returns `AlignmentResult.Empty` if null or empty. |
| `sequence2` | `DnaSequence` or `string` | required | Second sequence to align locally. | Same handling as `sequence1`. |
| `scoring` | `ScoringMatrix` | `SequenceAligner.SimpleDna` | Scoring values used for matches, mismatches, and gaps. | Optional in both public overloads. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `AlignedSequence1` | `string` | Aligned subsequence from input 1, with gap characters where needed. |
| `AlignedSequence2` | `string` | Aligned subsequence from input 2, with gap characters where needed. |
| `Score` | `int` | Non-negative local-alignment score under the chosen scoring matrix. |
| `AlignmentType` | `AlignmentType` | `Local` for successful local alignments; `AlignmentResult.Empty` carries `Global` because it is a shared sentinel value. |
| `StartPosition1` | `int` | 0-based inclusive start index of the aligned subsequence in input 1. |
| `StartPosition2` | `int` | 0-based inclusive start index of the aligned subsequence in input 2. |
| `EndPosition1` | `int` | 0-based inclusive end index of the aligned subsequence in input 1. |
| `EndPosition2` | `int` | 0-based inclusive end index of the aligned subsequence in input 2. |

### 3.3 Preconditions and Validation

`SequenceAligner.LocalAlign(DnaSequence, DnaSequence, ...)` throws `ArgumentNullException` when either `DnaSequence` argument is null. `DnaSequence` itself uppercases and validates `A`, `C`, `G`, and `T`, throwing `ArgumentException` on invalid characters. `SequenceAligner.LocalAlign(string, string, ...)` uppercases raw strings and returns `AlignmentResult.Empty` when either input is null or empty; it does not apply `DnaSequence` validation. Empty `DnaSequence` values are still legal inputs and go through the normal dynamic-programming path rather than the raw-string sentinel path. Start and end coordinates in successful results are 0-based.

## 4. Algorithm

### 4.1 High-Level Steps

1. Initialize the score matrix with a zero first row and zero first column.
2. Fill each interior cell with the maximum of zero, diagonal, up, and left candidate scores.
3. Track the highest-scoring cell while filling the matrix.
4. Start traceback from that maximum-scoring cell.
5. Stop traceback when a zero-valued cell is reached and return the aligned subsequences with their source coordinates.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository implementation uses `scoring.GapExtend` as a linear gap penalty inside the recurrence. During traceback, ties are resolved deterministically by preferring diagonal moves, then up moves, then left moves. The implementation records the matrix coordinates of the maximal cell and converts them to 0-based inclusive source positions when building `AlignmentResult`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Local alignment of sequences of lengths `m` and `n` | `O(mn)` | `O(mn)` | The repository stores the full score matrix and performs one traceback from the maximal cell. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)`: validated public API for `DnaSequence` inputs.
- `SequenceAligner.LocalAlign(string, string, ScoringMatrix?)`: raw-string overload with null and empty short-circuit behavior.
- `SequenceAligner.LocalAlignCore(string, string, ScoringMatrix)`: fills the Smith-Waterman score matrix and tracks the maximum-scoring cell.
- `SequenceAligner.TracebackLocal(...)`: reconstructs the local alignment and source positions.

**Supporting types:** [AlignmentTypes.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

### 5.2 Current Behavior

The repository implements the linear-gap form of Smith-Waterman with `Math.Max(0, ...)` in the cell update. `LocalAlignCore` keeps the first strictly larger maximum it encounters while scanning the matrix in row-major order; `TracebackLocal` then reconstructs one alignment by testing diagonal, then up, then left predecessors. The string overload uppercases its inputs and returns `AlignmentResult.Empty` when either input is null or empty. Empty `DnaSequence` inputs are not short-circuited, so the typed overload can return a `Local` result with empty aligned strings and coordinate fields derived as `-1` from the zero-size traceback endpoints. The coordinate fields in `AlignmentResult` are otherwise populated as 0-based inclusive indices.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Zero-initialized first row and first column.
- Smith-Waterman zero floor for every matrix cell.
- Traceback from the matrix maximum until a zero-valued cell is reached.

**Intentionally simplified:**

- The repository uses `GapExtend` as a linear gap penalty and does not evaluate the general $W_k$ gap-cost form; **consequence:** separate gap-open and gap-extension behavior is not modeled in `LocalAlign`.

**Not implemented:**

- Affine-gap or arbitrary-gap-cost Smith-Waterman variants; **users should rely on:** no current alternative in the pairwise local-alignment API.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty raw string input | Returns `AlignmentResult.Empty`. | The string overload exits before matrix construction on null or empty input. |
| Empty `DnaSequence` input | Returns a typed `Local` result with empty aligned strings and `-1` coordinates. | The typed overload delegates to the core routine instead of using the raw-string sentinel path. |
| Null `DnaSequence` input | Throws `ArgumentNullException`. | Public `DnaSequence` overloads guard with `ThrowIfNull`. |
| Identical sequences | Returns the full sequence with no gaps and score equal to the total match reward. | The full diagonal dominates under the repository's positive-match / negative-gap scoring examples. |
| Completely dissimilar sequences under negative mismatch and gap penalties | Returns score `0` and an empty alignment. | The zero floor prevents negative local scores from propagating. |
| Known Wikipedia example (`TGTTACGG` vs `GGTTGACTA`, score model `+3/-3/-2`) | Returns score `13` with alignment `GTT-AC` / `GTTGAC`. | Validated directly in `SequenceAligner_LocalAlign_Tests`. |

### 6.2 Limitations

The implementation returns one optimal local alignment, not the full set of equally optimal tracebacks. The algorithm keeps the full `O(mn)` score matrix and does not implement the linear-space refinements mentioned in the literature. The public string overload does not validate the alphabet beyond uppercasing its input. The pairwise local-alignment API does not expose affine-gap scoring.

## 7. Examples and Related Material

### 7.1 Worked Example

The repository tests validate the standard `TGTTACGG` versus `GGTTGACTA` example with `Match = +3`, `Mismatch = -3`, and `GapExtend = -2`. The reported optimal local alignment is:

```text
GTT-AC
GTTGAC
```

with score `13`, matching the expected Smith-Waterman example used in the test suite.

## 8. References

1. [Smith-Waterman algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
2. [Sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment)
3. Smith, T. F.; Waterman, M. S. (1981). "Identification of Common Molecular Subsequences." Journal of Molecular Biology 147(1): 195-197.
