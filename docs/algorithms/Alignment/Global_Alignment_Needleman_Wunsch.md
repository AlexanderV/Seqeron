# Global Alignment (Needleman-Wunsch)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-GLOBAL-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Production (linear + affine, matrix scoring, ends-free policies via `SemiGlobalAlign`, linear-space, co-optimal enumeration) |
| Last Reviewed | 2026-11-19 |

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
| `sequence1` | `DnaSequence` or `string` | required | First sequence to align. | Non-null (null throws `ArgumentNullException` on every overload). An empty sequence is legal: the result is the all-gap alignment (see §6.1). |
| `sequence2` | `DnaSequence` or `string` | required | Second sequence to align. | Same handling as `sequence1`. |
| `scoring` | `ScoringMatrix` | `SequenceAligner.SimpleDna` | Match, mismatch, and gap values used by the dynamic program. | Optional in all public overloads; `null` selects the default (the `ScoringMatrix?` overloads carry `[OverloadResolutionPriority(1)]`, so `GlobalAlign(a, b, null)` and `LocalAlignTopN(a, b, null, 3, 1)` compile - they used to be CS0121 against the `MatrixScoring` overloads). `GapOpen` and `GapExtend` must be <= 0 and greater than `int.MinValue` (a positive gap score is an `ArgumentException`, like `MatrixScoring`); the score of any path, `max(\|Match\|, \|Mismatch\|, \|GapOpen\| + \|GapExtend\|) * (m + n + 2)`, must stay below 2^28, otherwise `ArgumentException` "Score range too large for the integer dynamic program" (audit A2-37). |
| `cancellationToken` | `CancellationToken` | required on cancellation-aware overloads | Allows cancellation of the long-running matrix fill and traceback. | Used by both `GlobalAlign(DnaSequence, DnaSequence, ..., CancellationToken, ...)` and `GlobalAlign(string, string, ..., CancellationToken, ...)`. |
| `progress` | `IProgress<double>` | `null` | Optional progress reporter for the cancellation-aware overloads. | Values are reported during matrix fill and before/after traceback. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `AlignedSequence1` | `string` | First aligned sequence, including gap characters (`-`). |
| `AlignedSequence2` | `string` | Second aligned sequence, including gap characters (`-`). |
| `Score` | `int` | Alignment score under the chosen linear scoring matrix. |
| `AlignmentType` | `AlignmentType` | `Global` (also for empty inputs). |
| `StartPosition1` | `int` | `0`. |
| `StartPosition2` | `int` | `0`. |
| `EndPosition1` | `int` | `sequence1.Length - 1` (`-1` for an empty sequence). |
| `EndPosition2` | `int` | `sequence2.Length - 1` (`-1` for an empty sequence). |

### 3.3 Preconditions and Validation

`SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ...)` throws `ArgumentNullException` when either `DnaSequence` argument is null. `DnaSequence` construction normalizes input to uppercase and validates that each character is one of `A`, `C`, `G`, or `T`; invalid characters raise `ArgumentException`. `SequenceAligner.GlobalAlign(string, string, ...)` uppercases the supplied strings but does not apply `DnaSequence` validation; a null string throws `ArgumentNullException` (as a null `DnaSequence` does). **One contract for every `GlobalAlign`/`GlobalAlignAffine` overload (string, `DnaSequence`, cancellable): null throws `ArgumentNullException`; an empty sequence is not special-cased and is aligned by the normal DP.** Coordinates in `AlignmentResult` are 0-based.

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
| Global alignment of sequences of lengths `m` and `n` | `O(mn)` | `O(mn)` up to `LinearSpaceCellThreshold` = 25 000 000 cells, `O(m+n)` above | Standard Needleman-Wunsch complexity under the full dynamic-programming matrix; above the threshold the call switches automatically to the exact linear-space program (Hirschberg 1975 / Myers-Miller 1988, `GlobalAlignLinearSpace` / `GlobalAlignAffineLinearSpace`), identical rows, about 1.3-2x the time ([Linear_Space_Alignment.md](Linear_Space_Alignment.md)). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)`: validates non-null `DnaSequence` inputs and delegates to `GlobalAlignCore`.
- `SequenceAligner.GlobalAlign(string, string, ScoringMatrix?)`: uppercases raw strings and uppercases and delegates to `GlobalAlignCore` (null throws `ArgumentNullException`; empty is aligned normally).
- `SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?, CancellationToken, IProgress<double>?)`: validates non-null typed inputs and delegates to the same cancellation-aware `GlobalAlignCore`.
- `SequenceAligner.GlobalAlign(string, string, ScoringMatrix?, CancellationToken, IProgress<double>?)`: cancellation-aware variant that fills a 2D score matrix and reports progress.
- `SequenceAligner.GlobalAlignCore(string, string, ScoringMatrix)`: pooled-array implementation of the core dynamic program.
- `SequenceAligner.Traceback(...)`: reconstructs the final aligned strings from the score matrix.

**Supporting types:** [AlignmentTypes.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

### 5.2 Current Behavior

The repository exposes four linear-gap global-alignment entry points that share one core (`GlobalAlignCore`): it fills an `int[,]` matrix (checking cancellation every 100 rows and reporting progress for the cancellation-aware overloads) and calls `Traceback`. The traceback code is deterministic because it tests diagonal first, then up, then left. `ScoringMatrix.GapOpen` is not used by the linear dynamic program; the gap cost comes from `ScoringMatrix.GapExtend`. The raw-string overloads return `AlignmentResult.Empty` when either input is empty; the typed overloads (with or without cancellation) run the DP, yielding the all-gap border alignment `F(0,n) = n·d`. Affine gap costs are provided by `GlobalAlignAffine` (Gotoh 1982 three-state DP with the Flouri et al. 2015 border initialization; gap of length k = `GapOpen + k·GapExtend`, the NCBI BLAST existence + extension convention).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Linear-gap Needleman-Wunsch boundary initialization `F(i,0)=d*i` and `F(0,j)=d*j`.
- The standard three-way recurrence over diagonal, up, and left predecessors.
- Traceback from the bottom-right cell to recover one optimal global alignment.

**Intentionally simplified:**

- The public `ScoringMatrix` includes `GapOpen`, but `GlobalAlign` uses only `GapExtend` (linear model); the affine model is the separate `GlobalAlignAffine`.

**Ends-free / overlap alignment (A1-3, F44):** implemented as `SemiGlobalAlign(…, EndGaps)` / `SemiGlobalAlignAffine(…, EndGaps)` (all 16 end-gap policies, `EndGaps.All` = EMBOSS needle's default `-endweight N`; `EndGaps.None` = this global alignment); see [Semi_Global_Alignment.md §5.4](Semi_Global_Alignment.md).

**Arbitrary gap-cost functions (A1-5, F45):** implemented, §5.5.

**Not implemented:** none.

### 5.4 Substitution-matrix scoring and fractional gap scores (A1-2, B09 F41)

**Entry points:** `SequenceAligner.GlobalAlign` / `GlobalAlignAffine` / `LocalAlign` / `LocalAlignAffine` / `SemiGlobalAlign` / `SemiGlobalAlignAffine` with a `MatrixScoring(SubstitutionMatrix Matrix, double GapOpen, double GapExtend)` argument, for `string`, `ProteinSequence` and `DnaSequence` inputs (additive overloads; every `ScoringMatrix` overload, its behaviour and the tie order are unchanged). `MatrixScoring.Emboss(matrix, gapOpenPenalty = 10, gapExtendPenalty = 0.5)` builds the scoring from the EMBOSS needle/water numbers. `AlignmentResult.ExactScore` is the real-valued score (`Score` is that value rounded half away from zero; for every integer-scored aligner `ExactScore == Score`).

**One DP, no second implementation.** The matrix overloads convert the model to the exact integer model "every score = k/scale" (`PairScoring`: matrix cells x scale, gap scores resolved to 1e-6 and divided by their gcd with 10^d, so scale = 2 for the EMBOSS default 0.5 and 10 for 0.1) and run the *same* `LinearAlignCore` / `LocalAlignCore` / `AffineAlignCore` / `Traceback` as the integer aligners; the integer score is divided by the scale afterwards. Ties and tracebacks are therefore exact (no floating-point drift) and identical to the integer aligners (test `IntegerMatrix_IntegerGaps_IsIdenticalToScalarAligners_RowForRow`: 6 aligners x 400 random pairs, rows/score/coordinates equal to the scalar overloads when the matrix is "match 2 / mismatch -3"). Guards: a residue outside the matrix alphabet, a gap score > 0, NaN/infinite, finer than 1e-6, or a scaled score range that could overflow the `int` dynamic program (|score| x (m + n + 2) >= 2^28) throws `ArgumentException`.

**Gap convention.** `GapOpen + k*GapExtend` for a gap of length k (affine) / `GapExtend` per position (linear), both as negative scores. EMBOSS needle/water `-gapopen G -gapextend E` (positive penalties) charge G for a gap of length 1 and G + (k-1)E for length k, i.e. `GapOpen = -(G - E)`, `GapExtend = -E` (`MatrixScoring.Emboss`); the EMBOSS linear model is `gapopen == gapextend`. Biopython `open_gap_score = GapOpen + GapExtend`, `extend_gap_score = GapExtend`; parasail `open = -(GapOpen + GapExtend)`, `extend = -GapExtend`.

**needle's end-gap default.** `needle` leaves end gaps free on *both* sequences by default (`-endweight N`); the unit's global alignment charges them, which is `needle -endweight Y -endopen G -endextend E` (used for the oracle). The needle default (ends-free alignment of both sequences, `-endweight N`) is a different variant, implemented since audit A1-3 (B09 F44) as `SemiGlobalAlign` / `SemiGlobalAlignAffine` with `EndGaps.All` ([Semi_Global_Alignment.md](Semi_Global_Alignment.md)).

**Oracle** (`scripts/alignment/generate_matrix_align_oracle.py`, `TestData/Alignment/matrix_align_oracle.tsv`): 331 random protein (BLOSUM62/PAM250, incl. B/Z/X), DNA and IUPAC-DNA (EDNAFULL) pairs x {global, local, semi-global} x {linear, affine} = 1 986 rows with integer and fractional gap scores (0.5, 0.1, 0.25, 0.3, 0.75, 1.5, 2.5, ...), plus the HBA_HUMAN/HBB_HUMAN example. Biopython 1.88 `PairwiseAligner(substitution_matrix=...)` equals Seqeron on **all** 1 986 rows (1e-6), parasail 1.3.4 on all 616 integer-gap rows, EMBOSS needle (`-endweight Y`) / water 6.6.0 on 1 103 of 1 121 rows. The 18 EMBOSS disagreements are EMBOSS-internal, not Seqeron defects (Biopython and parasail, independent implementations, agree with Seqeron on every one): 14 needle rows are *lower* than the optimum (never higher: `embAlignPathCalcWithEndGapPenalties`, `nucleus/embaln.c`, opened, opens a last-row/last-column end gap only from the match matrix and never from the other gap state, so its end-gap-penalised DP is not exactly optimal) and 4 water rows differ (3 lower, 1 higher on EDNAFULL; water uses a different EMBOSS routine that was not analysed further). HBA/HBB with EMBOSS defaults (gapopen 10, gapextend 0.5, EBLOSUM62): global (affine) 292.5, local (water) 293.5, linear 366.5.

### 5.5 General gap-cost functions (Waterman-Smith-Beyer, A1-5, B09 F45)

**Entry points:** `SequenceAligner.GlobalAlignGeneralGap` / `LocalAlignGeneralGap` / `SemiGlobalAlignGeneralGap(…, EndGaps)` with a `Func<int,double> gapScore` (score of a gap of length k >= 1, negative = penalty) and either a `ScoringMatrix` (match/mismatch) or a `SubstitutionMatrix`; `string`, `DnaSequence` (and `ProteinSequence` with a matrix) overloads. File `SequenceAligner.GeneralGap.cs`.

**Model.** Waterman, Smith & Beyer (1976): a gap of length k scores an arbitrary $w(k)$ — concave (logarithmic, square root), convex (quadratic), capped, step, two-slope, even non-monotone. Three states $M, X, Y$ with $X_{i,j} = \max_{1 \le k \le i} \max(M, Y)_{i-k,j} + w(k)$, $Y_{i,j} = \max_{1 \le k \le j} \max(M, X)_{i,j-k} + w(k)$, $M_{i,j} = \max(M,X,Y)_{i-1,j-1} + s(a_i,b_j)$ (local: $\max(0,\cdot)$ on $M$; ends-free: the border/end rules of [Semi_Global_Alignment.md §2.2](Semi_Global_Alignment.md)). Every maximal gap run is charged $w(k)$ **once** — a run is never split, so a convex $w$ is honoured (the textbook single-matrix recurrence $H_{i,j} = \max(\ldots, \max_k H_{i-k,j} + w(k), \ldots)$ would cut a gap into cheaper pieces; measured against Biopython 1.88 on 200 random global pairs of length <= 9 per family it disagrees on 12 % (logarithmic), 57 % (quadratic) and 27 % (step function) of the rows, the three-state single-run model on none). Two runs of the two kinds may be adjacent and are charged separately. Complexity $O(mn(m+n))$ time (each cell scans every gap length), $O(mn)$ memory (back-pointers); `w` is evaluated once per length $1 \ldots \max(m,n)$ and must be pure. Ties: $M > X > Y$, then the shortest run (local: a non-positive predecessor is replaced by a fresh start).

**Fast path.** When $w(k) = w(1) + (k-1)e$ on $k = 1 \ldots \max(m,n)$ with $e \le 0$ and $w(1) - e \le 0$ (linear or affine, representable at 1e-6; integers for the scalar scoring) the call is delegated to the O(mn) programs of `GlobalAlign`/`GlobalAlignAffine` (open = 0 → linear), `LocalAlign(Affine)`, `SemiGlobalAlign(Affine)` — same score; the alignment is the one of that program (its tie order, identical to calling it directly). The general program gives the same score (tests).

**Oracle** (`scripts/alignment/generate_general_gap_oracle.py`, `TestData/Alignment/general_gap_oracle.tsv`, 1 440 rows): **Biopython 1.88** `PairwiseAligner(gap_score=callable)` (its "Waterman-Smith-Beyer global/local alignment algorithm"; `align().score` and `score()` agree on all rows) on 960 random global/local pairs (60 per function family x 8 families x 2 modes; scalar match/mismatch, BLOSUM62 and EDNAFULL scoring; lengths <= 26, including planted long deletions/insertions) = an independent three-state reference written for the oracle = Seqeron (general program and public API, 1e-8) on 960/960; 480 ends-free rows (policies 1..15) against the independent reference (Biopython 1.88 ignores the end-gap properties when `gap_score` is a callable); an exhaustive enumeration of every alignment of short sequences (60 pairs, 6 gap functions incl. a non-monotone one, all 16 policies, local on 15 pairs) equals the program.

**Limitations.** $O(mn(m+n))$ time: practical to a few hundred residues per sequence (affine/linear functions use the O(mn) fast path); doubles are compared exactly in the traceback (integer-valued scores are exact; a transcendental $w$ may in principle break an exact tie by a rounding error, giving an alignment within 1e-15 of the optimum).

### 5.6 Co-optimal global and ends-free alignments (A1-7b, B09 F59)

`SequenceAligner.GlobalAlignCoOptimal` / `GlobalAlignAffineCoOptimal` (and `SemiGlobalAlignCoOptimal` / `SemiGlobalAlignAffineCoOptimal` for any `EndGaps` policy; `ScoringMatrix` on `string` / `DnaSequence`, `MatrixScoring` on `string` / `ProteinSequence` / `DnaSequence`; [SequenceAligner.GlobalCoOptimal.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.GlobalCoOptimal.cs)) enumerate **every** optimal alignment, capped at `maxAlignments` (the count grows exponentially with tandem repeats and long gap-placement ties). MCP: `global_align_co_optimal`.

- **Definition.** An alignment is a gapped column sequence; its score is the sum of the substitution columns plus the gap cost (linear: `GapExtend` per gap position; affine: `GapOpen + k·GapExtend` per maximal run of gap columns of the same row, a gap in row 1 directly followed by a gap in row 2 being two runs); the leading / trailing gap columns of a free end (`EndGaps`) cost 0. The set is every column sequence of maximum score, each distinct pair of gapped rows once - the set Biopython's `PairwiseAligner(mode='global').align(a, b)` iterates.
- **No second DP.** The enumerator walks the very matrices of the single aligners: `LinearAlignCore` and `AffineAlignCore` were split into the shared helpers `LinearFillMatrix`, `AffineNewMatrices` / `AffineInitBorders` / `AffineFill` and `EndCellCandidates` (the Flouri border initialisation and the `AffineCell` recurrence are the same code; the local top-N enumerator `TopN.cs` now calls `AffineInitBorders` too), and only the traceback differs: a depth-first walk from every end cell of optimal value over **all** predecessors that reproduce a cell's value (diagonal > up > left; affine M > X > Y with the own gap state before the other), to the origin or to the first border cell of a free start. The first alignment is therefore the single aligner's (same end cell, same tie order).
- **Distinct rows.** With zero gap costs a free flank is also a run of zero-cost gap columns of the path (the same rows end at several cells); a set of already-listed rows removes the duplicates. Coordinates (`StartPosition/EndPosition`) are the consumed intervals of each alignment's rows (first / last column decide which side's flank is free); with zero gap costs and a free flank the single aligner reports the cell its traceback ends in, so the coordinates of the first alignment can differ from `SemiGlobalAlign`'s there (rows, score and type are identical).
- **Cost.** `O(m·n)` memory (the full matrices; `TopNMaxCells` = 25 000 000, larger input throws `ArgumentException`; one optimal alignment of larger input is `GlobalAlign`, which switches to linear space); the walk is output-sensitive.
- **Oracle** (`TestData/Alignment/global_cooptimal_oracle.tsv`, `scripts/alignment/generate_global_cooptimal_oracle.py`; Biopython 1.88): 1 920 pairs (14 scalar schemes + BLOSUM62 with fractional gaps; linear and affine; **all 16 end-gap policies**; long gaps, tandem repeats, zero gap costs) - the C# set equals Biopython's on every row with <= 150 alignments (1 859 rows, 935 with several alignments, up to 150 compared; the 61 larger rows check the capped prefix: distinct, all optimal), and on every row of **exhaustive enumeration** of every column sequence of pairs with <= 5 residues (written from the definition) for all 16 policies, linear and affine; every returned alignment re-scores to Biopython's optimum.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| One side empty (any overload) | The all-gap alignment: `""` vs `ACGT` -> `----` / `ACGT`, score `d·len` (`GapExtend` per position; `GlobalAlignAffine`: `GapOpen + len·GapExtend`). End coordinate of the empty side is `-1`. | Needleman-Wunsch border `F(0,j) = d·j`, `F(i,0) = d·i` (Needleman & Wunsch 1970; Durbin et al. 1998 §2.3). Biopython 1.88 `PairwiseAligner`, parasail 1.3.4, EMBOSS needle 6.6 and scikit-bio 0.7 reject zero-length input, so they are no oracle here; the recurrence, the Hirschberg base case and `ApproximateMatcher` indel cost agree. |
| Both sequences empty | Score `0`, empty rows, `Global`, coordinates `0, 0, -1, -1`. | Empty alignment of empty sequences. |
| Null input (string or `DnaSequence`, any overload) | Throws `ArgumentNullException`. | Every public overload guards with `ThrowIfNull`; no silent `AlignmentResult.Empty`. |
| Identical sequences | Returns a full-length gap-free alignment with score equal to the sum of match scores. | Verified by `SequenceAligner_GlobalAlign_Tests` and the Needleman-Wunsch recurrence. |
| Unequal-length sequences | Uses initialized border cells to absorb leading or trailing gaps in the optimal alignment. | This is the purpose of the boundary conditions `F(i,0)` and `F(0,j)`. |
| Completely different equal-length sequences under negative mismatch penalties | Can return a negative score. | Global alignment has no zero floor; every position still participates in the alignment. |

### 6.2 Limitations

`GlobalAlign` / `GlobalAlignAffine` return one optimal alignment (fixed tie order); the full co-optimal set is `GlobalAlignCoOptimal` / `GlobalAlignAffineCoOptimal` (§5.6, capped), the Waterman-Eggert suboptimal alignments exist for the local problem (`LocalAlignTopN`, [Local_Alignment_Smith_Waterman.md §5.5](Local_Alignment_Smith_Waterman.md)). The public string overload does not validate the alphabet beyond uppercasing the input. Up to 25 000 000 cells (`(m+1)(n+1)`, two 5 kb sequences) the full dynamic-programming matrix is kept (`O(mn)` memory, 4 bytes per cell linear / 12 affine); above that the aligners run the exact linear-space program on their own — Hirschberg (1975) for the linear model, Myers-Miller (1988) for affine gaps, same score, rows, coordinates and ties — with `O(m+n)` memory (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB; see [Linear_Space_Alignment.md](Linear_Space_Alignment.md)); `GlobalAlignLinearSpace` / `GlobalAlignAffineLinearSpace` force it for any size. Affine gap penalties are available through `GlobalAlignAffine`.

## 7. Examples and Related Material

### 7.1 Worked Example

With `Match = +1`, `Mismatch = -1`, and `GapExtend = -1`, the repository tests validate the standard `GCATGCG` versus `GATTACA` example with optimal score `0`. The exact aligned strings are allowed to vary among optimal solutions, but the tests require the global-alignment invariants and the final score described by the cited Needleman-Wunsch example.

## 8. References

1. [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
2. [Sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment)
