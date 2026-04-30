# Semi-Global Alignment (Fitting / Query-in-Reference)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-SEMI-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Semi-global alignment is the family of dynamic-programming alignments that leave one or more sequence ends unpenalized. The repository implementation is the fitting, or query-in-reference, variant: it aligns the entire query sequence against the best-scoring placement inside the reference while leaving leading and trailing reference context free. `SequenceAligner.SemiGlobalAlign(...)` therefore behaves like a modified Needleman-Wunsch algorithm rather than like local alignment with a zero floor. The implementation is exact for its linear-gap recurrence, but it does not expose the broader overlap or all-ends-free semi-global variants.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Semi-global alignment is appropriate when one sequence should be aligned in full but the other may contribute only a best-matching region. The current document's references describe this as a hybrid between global and local alignment: the query is aligned globally, while selected ends of the reference are left unpenalized. This is the classical use case for matching a shorter sequence against a longer reference.

### 2.2 Core Model

The repository's fitting variant keeps the standard Needleman-Wunsch recurrence

$$
F_{i,j} = \max\left(F_{i-1,j-1} + s(a_i, b_j),\; F_{i-1,j} + d,\; F_{i,j-1} + d\right)
$$

but changes the initialization and traceback start:

$$
F_{0,j} = 0, \qquad F_{i,0} = d \cdot i, \qquad \text{score} = \max_j F_{m,j}
$$

The zero first row leaves leading reference gaps unpenalized, the first column still penalizes gaps needed to align the full query, and the final score is taken from the maximum cell in the last row rather than from the bottom-right corner. Unlike Smith-Waterman, this recurrence has no zero floor, so the optimal fitting score may be negative. (Sequence alignment; Needleman-Wunsch algorithm; Rosalind SIMS; Rosalind SMGB)

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The query is fully represented after removing gaps from `AlignedSequence1`. | Fitting alignment aligns the query globally rather than allowing it to drop unmatched ends. |
| INV-02 | The two aligned output strings have equal length. | Traceback emits one alignment column per step, adding gaps where required. |
| INV-03 | The reported score is the maximum value in the last row of the score matrix. | Fitting alignment starts traceback from `max_j F(m,j)`. |
| INV-04 | The score may be negative. | The recurrence does not clamp cell values to zero. |

### 2.5 Comparison with Related Methods

| Variant | Free end gaps | Typical use |
|---------|---------------|-------------|
| Query-in-reference (implemented here) | Start and end of the reference | Short query against a longer reference |
| Overlap alignment | One suffix and one prefix | Sequence overlap / assembly-style problems |
| Full semi-global | All four ends | General ends-free alignment |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence1` | `DnaSequence` | required | Query sequence, typically the shorter sequence. | Must be non-null. Validation and uppercasing come from `DnaSequence`. |
| `sequence2` | `DnaSequence` | required | Reference sequence, typically the longer sequence. | Must be non-null. Validation and uppercasing come from `DnaSequence`. |
| `scoring` | `ScoringMatrix` | `SequenceAligner.SimpleDna` | Match, mismatch, and gap values used in the fitting dynamic program. | Optional. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `AlignedSequence1` | `string` | Query alignment string, including gap characters. |
| `AlignedSequence2` | `string` | Reference alignment string, including the matched region and any preserved unmatched reference context. |
| `Score` | `int` | Fitting-alignment score taken from the maximum of the last matrix row. |
| `AlignmentType` | `AlignmentType` | `SemiGlobal`. |
| `StartPosition1` | `int` | Current implementation returns `0`. |
| `StartPosition2` | `int` | Current implementation returns `0`. |
| `EndPosition1` | `int` | Current implementation returns `sequence1.Length - 1`. |
| `EndPosition2` | `int` | Current implementation returns `sequence2.Length - 1`. |

### 3.3 Preconditions and Validation

`SequenceAligner.SemiGlobalAlign(...)` accepts only `DnaSequence` inputs and throws `ArgumentNullException` when either argument is null. `DnaSequence` uppercases input and validates `A`, `C`, `G`, and `T`, throwing `ArgumentException` on invalid characters. The returned coordinate fields are 0-based, but for semi-global alignment the current implementation leaves them at full-input bounds rather than reporting the fitted interval explicitly.

## 4. Algorithm

### 4.1 High-Level Steps

1. Initialize the first row to zero and the first column to cumulative linear gap penalties.
2. Fill the score matrix with the standard three-way Needleman-Wunsch recurrence.
3. Find the maximum score in the last row.
4. Start traceback from that last-row maximum.
5. Preserve any unmatched trailing reference suffix as free end gaps in the query alignment and continue traceback until the matrix origin is reached.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository implementation uses `scoring.GapExtend` as the linear gap penalty and ignores `GapOpen` inside the semi-global dynamic program. `Traceback` is shared with the global-alignment path, but for semi-global results it first appends the unmatched trailing reference suffix as characters in `AlignedSequence2` against gaps in `AlignedSequence1`. During predecessor selection, ties are resolved by the same deterministic order used for global alignment: diagonal, then up, then left.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Fitting alignment of query length `m` against reference length `n` | `O(mn)` | `O(mn)` | The repository builds and stores the full DP matrix before traceback. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)`: public semi-global fitting entry point.
- `SequenceAligner.SemiGlobalAlignCore(string, string, ScoringMatrix)`: constructs the fitting-alignment score matrix and finds the last-row maximum.
- `SequenceAligner.Traceback(...)`: shared traceback routine used for both global and semi-global alignment.

**Supporting types:** [AlignmentTypes.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

### 5.2 Current Behavior

The implementation is specifically query-in-reference fitting alignment. The first row is left at zero, the first column is filled with cumulative `GapExtend` penalties, and the score returned in `AlignmentResult` is the maximum value in the last row. `Traceback` preserves unmatched reference suffix characters by placing them in `AlignedSequence2` against gaps in `AlignedSequence1`, and it also preserves unmatched reference prefix characters when the traceback reaches the top row. As a result, `AlignedSequence2` may retain the full reference string with unmatched context shown explicitly, while `StartPosition*` and `EndPosition*` remain full-input bounds rather than fitted-interval coordinates.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Fitting-alignment initialization with `F(0,j) = 0` and `F(i,0) = d * i`.
- Standard Needleman-Wunsch recurrence without a Smith-Waterman-style zero floor.
- Selection of the optimal fitting score from the maximum of the last row.

**Intentionally simplified:**

- The repository implements the query-in-reference fitting subset of semi-global alignment rather than the broader overlap or all-ends-free family; **consequence:** `SemiGlobalAlign` is targeted at fully aligning the query against a best-scoring placement in the reference.
- The returned coordinate fields do not identify the fitted interval; **consequence:** users must inspect the aligned strings and the score, not `StartPosition*` / `EndPosition*`, to locate the fitted region.

**Not implemented:**

- Overlap alignment and fully ends-free semi-global alignment; **users should rely on:** no current alternative in the repository's public pairwise alignment API.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null query or reference | Throws `ArgumentNullException`. | Public API guard. |
| Query embedded exactly inside the reference | Returns the embedded match score from the last-row maximum. | Validated by `SequenceAligner_SemiGlobalAlign_Tests`. |
| Query at the start or end of the reference | Leaves the unmatched reference prefix or suffix free. | This is the point of the fitting initialization and last-row maximum. |
| All mismatches under negative mismatch penalties | Can return a negative score. | The recurrence has no zero floor. |
| Query longer than the reference | Still aligns the full query, possibly introducing gaps into the reference alignment. | The query is always globally represented in the fitting variant. |

### 6.2 Limitations

The repository implements only one member of the broader semi-global family. The coordinate fields in `AlignmentResult` do not expose the fitted interval directly. Like the global path, the method uses a full `O(mn)` matrix and a linear gap model based on `GapExtend`. The public API does not offer a raw-string overload for semi-global alignment.

## 8. References

1. [Sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment)
2. [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
3. [Rosalind: Finding a Motif with Modifications (SIMS)](https://rosalind.info/problems/sims/)
4. [Rosalind: Semiglobal Alignment (SMGB)](https://rosalind.info/problems/smgb/)
5. Brudno, M. et al. (2003). "Glocal alignment: finding rearrangements during alignment." Bioinformatics 19 Suppl 1: i54-i62.
