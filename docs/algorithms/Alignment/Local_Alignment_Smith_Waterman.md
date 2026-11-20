# Local Alignment (Smith-Waterman)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-LOCAL-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Complete |
| Last Reviewed | 2026-11-19 |

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

This linear-gap form is implemented by `LocalAlign`. The affine form $W_k = o + k\,e$ (Gotoh 1982) is implemented by `LocalAlignAffine` with three states: $M_{i,j} = \max(0, M_{i-1,j-1}, X_{i-1,j-1}, Y_{i-1,j-1}) + s(a_i,b_j)$, $X_{i,j} = \max(M_{i-1,j}+o+e,\; X_{i-1,j}+e,\; Y_{i-1,j}+o+e)$, $Y_{i,j}$ symmetric, all borders $-\infty$, score $= \max M_{i,j}$ (o = `GapOpen`, e = `GapExtend`; Biopython `open_gap_score = o+e`, parasail `open = -(o+e)`). (Smith-Waterman algorithm; Gotoh 1982)

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
| `sequence1` | `DnaSequence` or `string` | required | First sequence to align locally. | Non-null (null throws `ArgumentNullException` on every overload). Empty is legal and gives the empty local result (§6.1). |
| `sequence2` | `DnaSequence` or `string` | required | Second sequence to align locally. | Same handling as `sequence1`. |
| `scoring` | `ScoringMatrix` | `SequenceAligner.SimpleDna` | Scoring values used for matches, mismatches, and gaps. | Optional in both public overloads. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `AlignedSequence1` | `string` | Aligned subsequence from input 1, with gap characters where needed. |
| `AlignedSequence2` | `string` | Aligned subsequence from input 2, with gap characters where needed. |
| `Score` | `int` | Non-negative local-alignment score under the chosen scoring matrix. |
| `AlignmentType` | `AlignmentType` | `Local` — also for the empty result (no positive-scoring pair), which is `("", "", 0, Local, -1, -1, -1, -1)`, not the `Global`-tagged `AlignmentResult.Empty` sentinel. |
| `StartPosition1` | `int` | 0-based inclusive start index of the aligned subsequence in input 1. |
| `StartPosition2` | `int` | 0-based inclusive start index of the aligned subsequence in input 2. |
| `EndPosition1` | `int` | 0-based inclusive end index of the aligned subsequence in input 1. |
| `EndPosition2` | `int` | 0-based inclusive end index of the aligned subsequence in input 2. |

### 3.3 Preconditions and Validation

`SequenceAligner.LocalAlign(DnaSequence, DnaSequence, ...)` throws `ArgumentNullException` when either `DnaSequence` argument is null. `DnaSequence` itself uppercases and validates `A`, `C`, `G`, and `T`, throwing `ArgumentException` on invalid characters. `SequenceAligner.LocalAlign(string, string, ...)` uppercases raw strings and throws `ArgumentNullException` when either input is null (empty input flows through the DP and yields the tagged-`Local` empty result); it does not apply `DnaSequence` validation. Empty `DnaSequence` values are still legal inputs and go through the normal dynamic-programming path rather than the raw-string sentinel path. Start and end coordinates in successful results are 0-based.

## 4. Algorithm

### 4.1 High-Level Steps

1. Initialize the score matrix with a zero first row and zero first column.
2. Fill each interior cell with the maximum of zero, diagonal, up, and left candidate scores.
3. Track the highest-scoring cell while filling the matrix.
4. Start traceback from that maximum-scoring cell.
5. Stop traceback when a zero-valued cell is reached and return the aligned subsequences with their source coordinates.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository implementation uses `scoring.GapExtend` as a linear gap penalty inside the recurrence. The end cell is the first maximal cell in row-major order (smallest end in input 1, then input 2; parasail's `sw` instead prefers the smallest end in input 2). During traceback, ties are resolved deterministically by preferring diagonal moves, then up moves, then left moves, and the traceback stops at the first zero cell, so the reported alignment is always one of the co-optimal alignments Biopython's local `PairwiseAligner` enumerates (verified on 3192 random/edge cases, 2026-09 review). The implementation records the matrix coordinates of the maximal cell and converts them to 0-based inclusive source positions when building `AlignmentResult`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Local alignment of sequences of lengths `m` and `n` | `O(mn)` | `O(mn)` up to 25 000 000 cells, `O(m+n)` above | Full score matrix and one traceback from the maximal cell up to `LinearSpaceCellThreshold`; above it the exact linear-space program (locate the optimal interval with an `O(n)`-memory forward pass, then solve it globally: Huang-Hardison-Miller 1990 / SIM, EMBOSS `matcher`) runs automatically with identical rows ([Linear_Space_Alignment.md](Linear_Space_Alignment.md)). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)`: validated public API for `DnaSequence` inputs.
- `SequenceAligner.LocalAlign(string, string, ScoringMatrix?)`: raw-string overload (null throws `ArgumentNullException`; empty is aligned normally).
- `SequenceAligner.LocalAlignCore(string, string, ScoringMatrix)`: fills the Smith-Waterman score matrix and tracks the maximum-scoring cell.
- `SequenceAligner.TracebackLocal(...)`: reconstructs the local alignment and source positions.
- `SequenceAligner.LocalAlignAffine(DnaSequence|string, DnaSequence|string, ScoringMatrix?)`: Smith-Waterman-Gotoh affine-gap local alignment; shares the three-state core `AffineAlignCore` with `GlobalAlignAffine`.
- `SequenceAligner.LocalAlignTopN` / `LocalAlignAffineTopN(seq1, seq2, scoring, maxAlignments, minScore)` ([SequenceAligner.TopN.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TopN.cs)): Waterman-Eggert non-intersecting local alignments; `LocalAlignCoOptimal` / `LocalAlignAffineCoOptimal(seq1, seq2, scoring, maxAlignments)`: all co-optimal local alignments (capped); `AlignedPairs(alignment)`: the forbidden pairs. See §5.5.

- **Further entry points (partial files of `SequenceAligner`):** `LocalAlignAffine` (Smith-Waterman-Gotoh, shared affine core with the global aligner), the substitution-matrix overloads with `MatrixScoring` (Global page §5.4), `LocalAlignTopN` / the Waterman-Eggert declumped alignments including the exact linear-space variant above 25 M cells (`SequenceAligner.TopN.cs`, `SequenceAligner.TopNLinearSpace.cs`; audit A1-7 / F77), co-optimal local alignments (`LocalAlignCoOptimal`, `LocalAlignAffineCoOptimal`; `SequenceAligner.TopN.cs`), `SequenceAligner.Cancellation.cs` (cancellable overloads, audit A2-40).

**Supporting types:** [AlignmentTypes.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

### 5.2 Current Behavior

The repository implements the linear-gap form of Smith-Waterman with `Math.Max(0, ...)` in the cell update. `LocalAlignCore` keeps the first strictly larger maximum it encounters while scanning the matrix in row-major order; `TracebackLocal` then reconstructs one alignment by testing diagonal, then up, then left predecessors. The string overload uppercases its inputs and returns `AlignmentResult.Empty` when either input is null or empty. Empty `DnaSequence` inputs are not short-circuited, so the typed overload can return a `Local` result with empty aligned strings and coordinate fields derived as `-1` from the zero-size traceback endpoints. The coordinate fields in `AlignmentResult` are otherwise populated as 0-based inclusive indices.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Zero-initialized first row and first column.
- Smith-Waterman zero floor for every matrix cell.
- Traceback from the matrix maximum until a zero-valued cell is reached.

**Model split:**

- `LocalAlign` uses `GapExtend` as a linear gap penalty (`GapOpen` ignored); affine gap costs (`GapOpen + k·GapExtend`) are provided by `LocalAlignAffine` (Smith-Waterman-Gotoh), cross-checked against Biopython 1.88 and parasail 1.3.4.

**Arbitrary gap-cost functions (A1-5, 2026-11, F45):** `SequenceAligner.LocalAlignGeneralGap(…, Func<int,double> gapScore, …)` implements the general $W_k$ form (Waterman-Smith-Beyer 1976) with the single-run model of Biopython's callable `gap_score` (every maximal gap run is charged $w(k)$ once, so a convex $w$ is honoured and gaps are never split into cheaper pieces); see [Global_Alignment_Needleman_Wunsch.md §5.5](Global_Alignment_Needleman_Wunsch.md). Linear and affine functions are recognised and run on the O(mn) `LocalAlign` / `LocalAlignAffine` programs.

**Not implemented:** none.

### Substitution-matrix scoring (A1-2, B09 F41)

`SequenceAligner.LocalAlign` and `LocalAlignAffine` also accept a `MatrixScoring` (BLOSUM62 / PAM250 / EDNAFULL residue-pair scores, fractional gap scores; `string`, `ProteinSequence`, `DnaSequence`) and run the same dynamic-program core as the integer overloads; see [Global_Alignment_Needleman_Wunsch.md §5.4](Global_Alignment_Needleman_Wunsch.md) for the conversion, the EMBOSS gap convention (`MatrixScoring.Emboss`) and the oracle (EMBOSS water 6.6.0, Biopython 1.88, parasail 1.3.4).

### 5.5 Suboptimal (Waterman-Eggert) and co-optimal local alignments (A1-7, B09 F47)

**Waterman-Eggert declumping** (Waterman & Eggert 1987, J Mol Biol 197:723-728; the algorithm of FASTA `lalign` and of EMBOSS `matcher -alternatives`): `LocalAlignTopN` / `LocalAlignAffineTopN` report the best local alignment, then the best local alignment that shares no *aligned pair* with it, and so on, up to `maxAlignments` alignments of score >= `minScore` (always > 0; scores are non-increasing).

- **Exact definition.** An *aligned pair* is a diagonal column, `a_i` opposite `b_j` (match or mismatch); a column that aligns a residue with a gap is not a pair. Alignment k is a maximum-score local alignment among those containing none of the pairs of alignments 1..k-1 (`SequenceAligner.AlignedPairs` returns them, 1-based). Two reported alignments can therefore cross the same matrix cell through a gap column but never align the same two residues. This is the definition implemented by `matcher`'s `MATCHERDIAG` macro (a pair used by an earlier alignment is removed from the diagonal move only, gap moves are untouched), checked below.
- **Algorithm (the 1987 recomputation scheme).** The full DP matrices (linear: H; affine: M, X, Y) are kept. After an alignment is reported its pairs are forbidden - the diagonal / M-state move into that cell is removed - and only the cells that can change are recomputed: a cell is recomputed when it is a new forbidden pair or one of its diagonal / upper / left predecessors changed, row by row, until a row has no change. The recurrences are the very functions the plain aligners fill the matrix with (`LocalLinearCell`, `AffineCell`) and the traceback is that of `LocalAlign` / `LocalAlignAffine` (`TracebackLocal`, `AffineTraceback`), so alignment 1 is exactly `LocalAlign` / `LocalAlignAffine`. Best cell per row is cached, so selecting the next best cell is `O(m)`.
- **Ties.** The alignment chosen among co-optimal ones is the fixed `LocalAlign` choice (first maximal cell in row-major order, diagonal > up > left / M > X > Y, stop at the first zero cell), applied to the matrix with the earlier pairs removed. A later round depends on which co-optimal alignment was taken earlier, so another tool may report a different (equally valid) sequence on tie-rich input; every reported score is the true maximum of the remaining problem. **`TopNTieOrder.Sim`** (A2-28d, F145) reproduces EMBOSS `matcher`'s own order exactly: a line-by-line port of `matcher.c` (K-vertex LIST keyed by start point, forward `BigPass`, `Locate` / `SmallPass` recomputation, linear-space `Diff`), 1 138 of 1 138 oracle lists equal alignment by alignment (the 59 further rows are those where the program aborts because fewer than K alignments exist; the port returns the alignments found); the list depends on K as in the program; integer scores only (a matrix with fractional gaps is converted to its exact integer model).
- **Scoring/overloads.** `ScoringMatrix` (string, `DnaSequence`) and `MatrixScoring` (string, `ProteinSequence`, `DnaSequence`; fractional gaps, `ExactScore`, real-valued `minScore`), linear and affine.
- **Memory.** `O(m*n)` (4 bytes per cell linear / 12 affine) as in the original algorithm up to `TopNMaxCells` (25 000 000 cells). Above it (A1-7c, B09 F77) `LocalAlignTopN` / `LocalAlignAffineTopN` (all overloads) run the **linear-space Waterman-Eggert** program (`LocalAlignTopNLinearSpace` / `LocalAlignAffineTopNLinearSpace` force it for any size): round k is one run of the exact linear-space local program (`LinearSpaceAlign`, Huang-Miller locate-then-solve, see `Linear_Space_Alignment.md` §2.4) on the matrix in which the aligned pairs of rounds 1..k-1 have no diagonal move - the scheme of SIM / EMBOSS `matcher`. The forbidden pairs are threaded (absolute coordinates, per-row sorted columns) through `LocateLinearLocal` / `LocateAffine`, `LinearMidpoint` / `AffineMidpoint`, `LinearSolve` / `AffineSolve` and the full-matrix base cases, so the result is **identical** (rows, score, coordinates, order, ties) to the in-memory declumper; memory `O(m + n)` per round plus the forbidden pairs (the total length of the reported alignments), time one local alignment per reported alignment (+1 for the empty final round). It needs non-positive gap scores (`ArgumentException` above the limit with a positive gap score). 1.9 kb x 1.3 kb, 50 alignments in memory: 1.2-1.4 s.

**Co-optimal enumeration** (`LocalAlignCoOptimal` / `LocalAlignAffineCoOptimal`): every distinct local alignment of the optimum score - the set Biopython's `PairwiseAligner(mode='local').align(a, b)` iterates - capped at `maxAlignments` (the count grows exponentially with repeats). An alignment is a depth-first traceback path from a maximal cell reached by a diagonal move, following every predecessor that reproduces the value (diagonal > up > left / M > X > Y), until the predecessor value is 0; as in Biopython ("zero-score extensions to alignments are not allowed", `_pairwisealigner.c`) a path may not pass through another end cell of an optimal alignment. Order: end cell in row-major order, then the branch order above, so the first element is `LocalAlign`'s result.

**Oracle** (`TestData/Alignment/topn_oracle.tsv`, `scripts/alignment/generate_topn_oracle.py`, 1 350 pairs: tandem-repeat-rich DNA, repeated protein domains, random pairs; linear and affine; scalar scoring and BLOSUM62 with fractional gaps):

| Reference | Result |
|---|---|
| Independent Waterman-Eggert written from the definition (full recompute of the matrix after every alignment, Python) | identical alignments (score, coordinates, both rows) on all 1 350 rows, and on 300 further long pairs (up to 174 residues, up to 40 alignments, 6 405 alignments) |
| Exhaustive enumeration of every local alignment of tiny inputs (<= 5 residues, 198 + 600 pairs) | every round of every row is the exhaustive maximum among alignments avoiding the earlier pairs; the list stops exactly when none is left |
| EMBOSS `matcher` 6.6.0 `-alternatives k` (Huang-Miller SIM; `matcher_BigPass`/`matcher_Diff` source opened, `kimrutherford/EMBOSS` mirror) | 397 rows / 1 695 alignments (plus 300 more random pairs / 1 207 alignments in a first check): every alignment matcher prints re-scores to its reported score, uses no pair of an earlier one and is the optimum of the matrix with the earlier pairs removed (so matcher = the definition above); Seqeron's score list equals matcher's on 390 of 397 rows, the 7 others (lines 46, 487, 502, 787, 793, 865, 1348 of `topn_oracle.tsv`) diverge only through tie choice - `scripts/alignment/analyze_topn_matcher_ties.py` re-runs matcher and shows in each of them an alignment of the same score but different pairs in an earlier or the same round (audit A3-14c, B09 F143); with `TopNTieOrder.Sim` all 397 rows and 800 further pairs are exact (A2-28d, F145); round 1 equal on all 397 |
| Biopython 1.88 `align()` co-optimal set | identical set on all 1 345 rows with <= 60 co-optimal alignments |

The tool `matcher` reports a different alignment than `water` only after the first one; `water` reports only the best.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input (either side or both; string or `DnaSequence`; linear or affine) | Returns `("", "", 0, Local, -1, -1, -1, -1)`. | No positive-scoring pair exists, so the zero-floored matrix has maximum 0 and the traceback is empty; all overloads share this result (the former `string` overload returned the `Global`-tagged `AlignmentResult.Empty`). |
| Null input (string or `DnaSequence`) | Throws `ArgumentNullException`. | Every public overload guards with `ThrowIfNull`. |
| Identical sequences | Returns the full sequence with no gaps and score equal to the total match reward. | The full diagonal dominates under the repository's positive-match / negative-gap scoring examples. |
| Completely dissimilar sequences under negative mismatch and gap penalties | Returns score `0` and an empty alignment. | The zero floor prevents negative local scores from propagating. |
| Known Wikipedia example (`TGTTACGG` vs `GGTTGACTA`, score model `+3/-3/-2`) | Returns score `13` with alignment `GTT-AC` / `GTTGAC`. | Validated directly in `SequenceAligner_LocalAlign_Tests`. |

### 6.2 Limitations

`LocalAlign` / `LocalAlignAffine` return one optimal local alignment (fixed tie order); the suboptimal alignments of Waterman-Eggert and the full co-optimal set are `LocalAlignTopN` and `LocalAlignCoOptimal` (§5.5; `LocalAlignTopN` runs in linear space above 25 000 000 cells (§5.5), `LocalAlignCoOptimal` keeps the `O(mn)` matrices and rejects more than 25 000 000 cells - an exponential-size output has no linear-space form). Up to 25 000 000 cells the full `O(mn)` score matrix is kept; above that `LocalAlign` / `LocalAlignAffine` run the exact linear-space refinement (Huang-Hardison-Miller 1990 / Huang-Miller 1991, the algorithm of SIM and EMBOSS `matcher`: locate the optimal local interval in `O(m+n)` memory, then align it globally with Hirschberg / Myers-Miller) with the same score, rows, coordinates and ties (`LocalAlignLinearSpace` / `LocalAlignAffineLinearSpace` force it; they need non-positive gap scores). The public string overload does not validate the alphabet beyond uppercasing its input. Affine-gap scoring is available through `LocalAlignAffine`; like the linear version it expects non-positive gap scores.

## 7. Examples and Related Material

### 7.1 Worked Example

The repository tests validate the standard `TGTTACGG` versus `GGTTGACTA` example with `Match = +3`, `Mismatch = -3`, and `GapExtend = -2`. The reported optimal local alignment is:

```text
GTT-AC
GTTGAC
```

with score `13`, matching the expected Smith-Waterman example used in the test suite.

### 5.6 EMBOSS `dotmatcher` (A2-28d, B09 F152)

`DotMatcher.Compute(sequence1, sequence2, windowSize = 10, threshold = 23, matrix)` / MCP `dotmatcher`: the windowed dot plot of EMBOSS `dotmatcher` - the summed matrix score of a window of `windowSize` residue pairs slides along each diagonal and every maximal run of windows at or above `threshold` is one segment; 600 / 600 plots equal the program (the program draws the second sequence on the x axis).

## 8. References

1. [Smith-Waterman algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
2. [Sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment)
3. Smith, T. F.; Waterman, M. S. (1981). "Identification of Common Molecular Subsequences." Journal of Molecular Biology 147(1): 195-197.
4. Gotoh, O. (1982). "An improved algorithm for matching biological sequences." Journal of Molecular Biology 162: 705-708.
5. Reference implementations: Biopython `Bio/Align/_pairwisealigner.c` (GOTOH_LOCAL_SCORE), parasail `src/sw.c`.
6. Waterman, M. S.; Eggert, M. (1987). "A new algorithm for best subsequence alignments with application to tRNA-rRNA comparisons." Journal of Molecular Biology 197(4): 723-728.
7. Huang, X.; Hardison, R. C.; Miller, W. (1990). "A space-efficient algorithm for local similarities." CABIOS 6: 373-381; Huang, X.; Miller, W. (1991). Adv. Appl. Math. 12: 337-357. EMBOSS 6.6.0 `emboss/matcher.c` (SIM).
