# Semi-Global Alignment (Fitting, Overlap and Ends-Free Variants)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-SEMI-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Production |
| Last Reviewed | 2026-11-18 |

## 1. Overview

Semi-global alignment is the family of dynamic-programming alignments that leave one or more sequence ends unpenalized. `SequenceAligner.SemiGlobalAlign(...)` / `SemiGlobalAlignAffine(...)` implement the whole family: each of the four end gaps (leading/trailing residues of sequence 1 and of sequence 2) is free or penalised according to an `EndGaps` flag set, with linear or affine (Gotoh) gap costs. The default (two/three-argument overloads) is the fitting, or query-in-reference, variant: it aligns the entire query sequence against the best-scoring placement inside the reference while leaving leading and trailing reference context free (`EndGaps.Fitting`). The overlap (suffix-prefix) variants in both orientations, sequence 2 fitted inside sequence 1, all four ends free (EMBOSS needle's default `-endweight N`, parasail `sg`), any other combination of the 16 policies, and `EndGaps.None` (= global alignment) use the same two dynamic programs (`LinearAlignCore`, `AffineAlignCore`) — only the border initialisation and the end cell change — so there is no second implementation. The algorithms behave like a modified Needleman-Wunsch algorithm rather than like local alignment with a zero floor.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Semi-global alignment is appropriate when one sequence should be aligned in full but the other may contribute only a best-matching region. The current document's references describe this as a hybrid between global and local alignment: the query is aligned globally, while selected ends of the reference are left unpenalized. This is the classical use case for matching a shorter sequence against a longer reference. Freeing the end gaps of both sequences in opposite orientation gives the *overlap* (dovetail, suffix-prefix) alignment of sequence assembly; freeing all four gives the general ends-free alignment (Rosalind SMGB; EMBOSS needle's default).

### 2.2 Core Model

The repository's fitting variant keeps the standard Needleman-Wunsch recurrence

$$
F_{i,j} = \max\left(F_{i-1,j-1} + s(a_i, b_j),\; F_{i-1,j} + d,\; F_{i,j-1} + d\right)
$$

but changes the initialization and traceback start:

$$
F_{0,j} = 0, \qquad F_{i,0} = d \cdot i, \qquad \text{score} = \max_j F_{m,j}
$$

Only end gaps *in the query* (overhanging reference residues) are free; a query overhang is charged like any other gap. Reference-implementation equivalents: Biopython `PairwiseAligner(mode="global", end_insertion_score=0)` with target = query, and parasail `sg_dx` (s1 = query) — both confirmed on 2428 random/edge cases in the 2026-09 review.

The zero first row leaves leading reference gaps unpenalized, the first column still penalizes gaps needed to align the full query, and the final score is taken from the maximum cell in the last row rather than from the bottom-right corner. Unlike Smith-Waterman, this recurrence has no zero floor, so the optimal fitting score may be negative. (Sequence alignment; Needleman-Wunsch algorithm; Rosalind SIMS; Rosalind SMGB)

**Ends-free generalisation (`EndGaps`).** Four independent flags `FreeStart1`, `FreeEnd1`, `FreeStart2`, `FreeEnd2` (leading/trailing residues of sequence 1 / sequence 2 that overhang the alignment cost nothing). The recurrence is unchanged; the initialisation is $F_{i,0} = 0$ if `FreeStart1` else $d\,i$, $F_{0,j} = 0$ if `FreeStart2` else $d\,j$, and the score is $F_{m,n}$ when no end is free, otherwise the maximum over the last row (`FreeEnd2`) and/or the last column (`FreeEnd1`). The affine program gives a free start the border $M_{i,0} = 0$ / $M_{0,j} = 0$ (penalised starts keep $X_{i,0} = o + i\,e$, $Y_{0,j} = o + j\,e$, Flouri et al. 2015) and takes the end maximum of $\max(M,X,Y)$. The fitting variant of this section is `FreeStart2 | FreeEnd2`.

| `EndGaps` | Free ends | Biopython 1.88 | parasail 1.3.4 | EMBOSS / textbook |
|-----------|-----------|----------------|----------------|-------------------|
| `None` | none | `mode='global'` | `nw` | needle `-endweight Y` |
| `Fitting` (default) | `FreeStart2`, `FreeEnd2` | `end_insertion_score = 0` | `sg_dx` | Rosalind SIMS |
| `FittingSequence2` | `FreeStart1`, `FreeEnd1` | `end_deletion_score = 0` | `sg_qx` | |
| `OverlapSuffix1Prefix2` | `FreeStart1`, `FreeEnd2` | `left_deletion_score = right_insertion_score = 0` | `sg_qb_de` | overlap / dovetail |
| `OverlapSuffix2Prefix1` | `FreeEnd1`, `FreeStart2` | `right_deletion_score = left_insertion_score = 0` | `sg_qe_db` | overlap, other orientation |
| `All` | all four | `end_gap_score = 0` | `sg` | needle default, Rosalind SMGB |
| single flags | one | `left/right_insertion/deletion_score = 0` | `sg_qb`, `sg_qe`, `sg_db`, `sg_de` | |

(sequence 1 = Biopython target = parasail s1/query; a free overhang of sequence 1 is a gap in sequence 2, a Biopython *deletion*.)

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The query is fully represented after removing gaps from `AlignedSequence1`. | Fitting alignment aligns the query globally rather than allowing it to drop unmatched ends. |
| INV-02 | The two aligned output strings have equal length. | Traceback emits one alignment column per step, adding gaps where required. |
| INV-03 | The reported score is the maximum value in the last row of the score matrix. | Fitting alignment starts traceback from `max_j F(m,j)`. |
| INV-04 | The score may be negative. | The recurrence does not clamp cell values to zero. |
| INV-05 | `EndGaps.Fitting` equals the classic fitting alignment and `EndGaps.None` equals the global alignment, row for row. | Same dynamic programs, same tie order. |
| INV-06 | Freeing an end gap never lowers the score; `All` is the maximum, `None` the minimum. | A free end only relaxes a constraint. |
| INV-07 | The aligned strings always contain both complete inputs; the free flanks are shown against `-`. | Traceback appends the unconsumed prefix/suffix. |

### 2.5 Comparison with Related Methods

| Variant | Free end gaps | Typical use |
|---------|---------------|-------------|
| Query-in-reference (`Fitting`, the default) | Start and end of the reference | Short query against a longer reference |
| Overlap alignment (`OverlapSuffix1Prefix2` / `OverlapSuffix2Prefix1`) | One suffix and one prefix | Sequence overlap / assembly-style problems |
| Full semi-global (`All`) | All four ends | General ends-free alignment; EMBOSS needle default |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence1` | `DnaSequence` | required | Query sequence, typically the shorter sequence. | Must be non-null. Validation and uppercasing come from `DnaSequence`. |
| `sequence2` | `DnaSequence` | required | Reference sequence, typically the longer sequence. | Must be non-null. Validation and uppercasing come from `DnaSequence`. |
| `scoring` | `ScoringMatrix` | `SequenceAligner.SimpleDna` | Match, mismatch, and gap values used in the fitting dynamic program. | Optional. |
| `endGaps` | `EndGaps` | `EndGaps.Fitting` (overloads without the parameter) | Which of the four end gaps are free: `FreeStart1`, `FreeEnd1`, `FreeStart2`, `FreeEnd2`, or the presets `None`, `Fitting`, `FittingSequence2`, `OverlapSuffix1Prefix2`, `OverlapSuffix2Prefix1`, `All`. | Overloads `SemiGlobalAlign(string\|DnaSequence, …, ScoringMatrix?, EndGaps)`, `SemiGlobalAlignAffine(…)` and `(string\|ProteinSequence\|DnaSequence, …, MatrixScoring, EndGaps)`. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `AlignedSequence1` | `string` | Query alignment string, including gap characters. |
| `AlignedSequence2` | `string` | Reference alignment string, including the matched region and any preserved unmatched reference context. |
| `Score` | `int` | Fitting-alignment score taken from the maximum of the last matrix row. |
| `AlignmentType` | `AlignmentType` | `SemiGlobal`; `Global` when `endGaps` is `None`. |
| `StartPosition1` | `int` | `0` (the query is aligned in full); with a free sequence-1 flank (`FreeStart1`) the number of free leading sequence-1 residues. |
| `StartPosition2` | `int` | First index of the reference consumed by the alignment = number of free leading reference residues (the fitted interval start). |
| `EndPosition1` | `int` | `sequence1.Length - 1` (`-1` for an empty query); with `FreeEnd1` the last sequence-1 index consumed by the path. |
| `EndPosition2` | `int` | Last index of the reference consumed by the alignment (the free trailing residues lie after it); an empty fitted interval (empty query, or query aligned entirely against gaps) has `EndPosition2 = StartPosition2 - 1`. |

For every `EndGaps` value the coordinates are the **consumed intervals** of the two sequences (0-based inclusive: from the first to the last residue between the free flanks; an empty interval has `End = Start − 1`); everything outside them is a free overhang shown against `-`. At most one of the two leading flanks and at most one of the two trailing flanks is non-empty. Overlap example (Biopython 1.88, match 1 / mismatch −1 / gap −1): `GGGGACGTAC` vs `ACGTACTTTT` with `OverlapSuffix1Prefix2` → `GGGGACGTAC----` / `----ACGTACTTTT`, score 6, coordinates (4, 9, 0, 5).

### 3.3 Preconditions and Validation

`SequenceAligner.SemiGlobalAlign(...)` (and `SemiGlobalAlignAffine`) accept `DnaSequence` or raw `string` inputs and throw `ArgumentNullException` when either argument is null; an empty side is not special-cased (empty query: score 0, the whole reference is free; empty reference: the query is one gap run). `DnaSequence` uppercases input and validates `A`, `C`, `G`, and `T`, throwing `ArgumentException` on invalid characters; the string overloads uppercase only. The returned coordinate fields are 0-based inclusive; `StartPosition2`/`EndPosition2` are the **fitted interval** of the reference (Biopython `alignment.coordinates`/`aligned` span between the free end gaps; parasail `end_ref`), `StartPosition1`/`EndPosition1` span the whole query.

## 4. Algorithm

### 4.1 High-Level Steps

1. Initialize the first row to zero (`FreeStart2`) or to cumulative gap penalties, the first column likewise (`FreeStart1`).
2. Fill the score matrix with the standard three-way Needleman-Wunsch recurrence.
3. Find the end cell: the bottom-right cell, or the first maximum of the last row (`FreeEnd2`, smallest j) and/or last column (`FreeEnd1`, smallest i; row candidates first).
4. Start traceback from that cell, emitting the unconsumed trailing flank against gaps.
5. Trace back (diagonal > up > left; affine M > X > Y) until the origin, or until row 0 / column 0 when the matching start is free, and emit the unconsumed leading flank against gaps.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository implementation uses `scoring.GapExtend` as the linear gap penalty and ignores `GapOpen` inside the semi-global dynamic program. `Traceback` is shared with the global-alignment path, but for semi-global results it first appends the unmatched trailing reference suffix as characters in `AlignedSequence2` against gaps in `AlignedSequence1`. During predecessor selection, ties are resolved by the same deterministic order used for global alignment: diagonal, then up, then left. Among tied last-row maxima the smallest column j wins, including j = 0 (the query aligned entirely against gaps, e.g. AAAA / CCCC with 1/−1/−1 → `AAAA----` / `----CCCC`, score −4).

**Affine variant (`SemiGlobalAlignAffine`).** Gotoh three-state recurrences shared with `GlobalAlignAffine`/`LocalAlignAffine` (`AffineAlignCore`); a gap of length k scores `GapOpen + k·GapExtend` (Biopython open_gap_score = GapOpen + GapExtend; parasail open = −(GapOpen + GapExtend)). Borders: M(0,j) = 0 for all j, X(i,0) = GapOpen + i·GapExtend, other border cells −∞; score = max_j max(M, X, Y)(m, j), first maximal column wins, state ties M > X > Y.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Fitting alignment of query length `m` against reference length `n` | `O(mn)` | `O(mn)` up to 25 000 000 cells, `O(m+n)` above | The full DP matrix is built before traceback up to `LinearSpaceCellThreshold`; above it (and for every `EndGaps` policy) the exact linear-space program runs automatically — locate the start vertex and the end cell with an `O(n)`-memory forward pass, then Hirschberg / Myers-Miller on the rectangle — with identical rows ([Linear_Space_Alignment.md](Linear_Space_Alignment.md)). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)`: public semi-global fitting entry point.
- `SequenceAligner.SemiGlobalAlignAffine(DnaSequence|string, …, ScoringMatrix?)`: affine-gap fitting alignment (via the shared Gotoh `AffineAlignCore`).
- `SequenceAligner.SemiGlobalAlign(string|DnaSequence, …, ScoringMatrix?, EndGaps)` / `SemiGlobalAlignAffine(…, EndGaps)` and the `MatrixScoring` overloads (`string`, `ProteinSequence`, `DnaSequence`): the ends-free family (A1-3, F44).
- `SequenceAligner.LinearAlignCore(string, string, PairScoring, EndGaps, …)`: linear-gap DP shared with `GlobalAlign`; the end-gap flags zero the first row/column and choose the end cell (last row and/or column maximum).
- `SequenceAligner.Traceback(...)`: shared traceback routine used for both global and ends-free alignment; `AffineAlignCore(…, EndGaps)` is its Gotoh counterpart.

- **Further entry points (partial files of `SequenceAligner`):** `SemiGlobalAlign(…, EndGaps)` / `SemiGlobalAlignAffine(…, EndGaps)` (all 16 end-gap policies, §5.4), `SemiGlobalAlignLinearSpace` / `SemiGlobalAlignAffineLinearSpace` (`SequenceAligner.LinearSpace.cs`), `SemiGlobalAlignCoOptimal` / `SemiGlobalAlignAffineCoOptimal` (`SequenceAligner.GlobalCoOptimal.cs`), `SemiGlobalAlignGeneralGap` (`SequenceAligner.GeneralGap.cs`), the substitution-matrix overloads (Global page §5.4).

**Supporting types:** [AlignmentTypes.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

### 5.2 Current Behavior

For the default `EndGaps.Fitting` the behaviour is the query-in-reference fitting alignment described here (other policies: §2.2 and §5.4). The first row is left at zero, the first column is filled with cumulative `GapExtend` penalties, and the score returned in `AlignmentResult` is the maximum value in the last row. `Traceback` preserves unmatched reference suffix characters by placing them in `AlignedSequence2` against gaps in `AlignedSequence1`, and it also preserves unmatched reference prefix characters when the traceback reaches the top row. As a result, `AlignedSequence2` may retain the full reference string with unmatched context shown explicitly, while `StartPosition2`/`EndPosition2` report the fitted interval (first/last reference index consumed between the free end gaps). Co-optimal fittings can occupy different intervals; the reported one is that of the returned alignment (see §6.1).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Fitting-alignment initialization with `F(0,j) = 0` and `F(i,0) = d * i`.
- Standard Needleman-Wunsch recurrence without a Smith-Waterman-style zero floor.
- Selection of the optimal fitting score from the maximum of the last row.

**Implemented since the first review (A1-3, 2026-11, F44):** the overlap (suffix-prefix) variants in both orientations and the general ends-free family — all 16 end-gap policies — via `EndGaps`, linear and affine, scalar and substitution-matrix scoring (§5.4). **Intentionally simplified:** none. **Not implemented:** none. Per-end and per-direction gap *scores* (Biopython's `left_gap_score` ... `{open,extend}_{internal,left,right}_{insertion,deletion}_score`) are implemented since audit round 3 (A3-10, B09 F137): see §5.5.

**Affine gaps:** available through `SemiGlobalAlignAffine` (Gotoh 1982; cross-checked against Biopython 1.88 and parasail 1.3.4 `sg_dx`).

### 5.5 Per-direction and per-end gap scores (A3-10, B09 F137)

`GapScoreProfile` (Infrastructure/AlignmentTypes.cs) carries, for each gap direction, an interior, a start (left) and an end (right) `GapScore(Open, Extend)`: `Interior1` / `Start1` / `End1` for a gap in sequence 2 opposite residues of sequence 1 (Biopython *deletion* = `query_*`), `Interior2` / `Start2` / `End2` for a gap in sequence 1 (*insertion* = `target_*`); a `null` member is the scoring's own gap score, `GapScore.Free` a free overhang, so `GapScoreProfile.FromEndGaps(flags)` is every `EndGaps` policy. `GapScoreProfile.FromBiopython(dictionary)` accepts both name families of Biopython 1.88 (`open_left_insertion_score`, ... and the deprecated `target_left_open_gap_score`, ...) and every aggregate of its property tree (`open_gap_score`, `end_extend_gap_score`, `left_gap_score`, `insertion_score`, ...; the more specific name wins). Entry points: `GlobalAlign` / `GlobalAlignAffine` / `SemiGlobalAlign` / `SemiGlobalAlignAffine` / `LocalAlign` / `LocalAlignAffine` (string, `ScoringMatrix?` | `MatrixScoring`, `GapScoreProfile`); MCP `align_gap_profile`.

*No second dynamic program.* A vertical run (gap in sequence 2) stays in one matrix column and a horizontal run in one row, so the existing linear (Needleman-Wunsch), Gotoh and Smith-Waterman cores take the gap scores as a function of the position (`PairScoring.XExt/XOpen/YExt/YOpen`): a run in column / row 0 is a start run, in the last column / row an end run, otherwise interior. Local alignments use the interior scores only (they never start or end in a gap; Biopython's local mode ignores the left / right scores). Linear gaps use the per-position `Extend`; affine gaps score `Open + k·Extend` (a Biopython pair converts as `Open = open - extend`). Scores may be fractional with a substitution matrix (exact integer scaling as for `MatrixScoring`); `Open` itself may be a bonus as long as `Open + Extend <= 0`. Coordinates: a free start / end run is excluded from `StartPosition*` / `EndPosition*` (the rows still show it), a priced one is inside. Inputs above `LinearSpaceCellThreshold` cells are aligned by the exact linear-space program with the same gap profile (A4-8, B09 F163; `GlobalAlignLinearSpace` / `GlobalAlignAffineLinearSpace` / `SemiGlobalAlign*LinearSpace` / `LocalAlign*LinearSpace` with a `GapScoreProfile`, scalar and `MatrixScoring`): a sub-problem of the Hirschberg / Myers-Miller recursion sits at (row offset, column offset) of the whole matrix and resolves its start / interior / end gap score by that position (`PairScoring.Window`), a local alignment uses the interior scores only (`PairScoring.InteriorOnly`); rows, score and coordinates are identical to the full-matrix program (600 forced-recursion cases with the base-case size 1 ... 300 cells and the 1 000 Biopython oracle rows with the base-case size 1; an 8 700 x 7 999 pair, 69.6 M cells, with Biopython scores).

**Oracle** (`scripts/alignment/generate_gap_profile_oracle.py`, `TestData/Alignment/gap_profile_oracle.tsv`): 1000 random pairs x global | local x linear | affine x scalar | EDNAFULL | BLOSUM62 (fractional gap scores) with random values for all twelve properties (free, equal-to-interior and independently priced ends); Biopython 1.88 `PairwiseAligner.score` = Seqeron on 998 / 1000 rows, an independent Python three-state reference on 1000 / 1000. **Biopython 1.88 corner (the 2 differing rows, source `Bio/Align/_pairwisealigner.c` read on raw.githubusercontent.com):** the LINEAR Needleman-Wunsch score initialises the first column of the last row as `row[0] = nA * right_gap_extend_B` (a leading run that consumes ALL of sequence 1 is charged like an end gap), whereas its Gotoh routine of the same file uses the LEFT score (`Ix_row[0] = left_gap_open_B + left_gap_extend_B * (i - 1)`); Seqeron follows the position rule. The generator reproduces the two rows with that rule.

### 5.4 Ends-free and overlap variants (A1-3, B09 F44)

**Oracle** (`scripts/alignment/generate_endgaps_oracle.py`, `TestData/Alignment/endgaps_oracle.tsv`, 3 590 rows): 70 random/structured DNA pairs (dovetail overlaps, containments, random) x 16 policies x {linear, affine} x random schemes = 2 240 scalar rows (alphabets A/AC/ACGT, GapOpen 0 and GapExtend 0 included) plus 1 350 substitution-matrix rows (BLOSUM62 / EDNAFULL, fractional gaps). **Biopython 1.88** `PairwiseAligner(left/right_insertion/deletion_score = 0)`, `align().score` = Seqeron on 3 590 / 3 590 rows (also = an independent exhaustive reference and a brute force over every alignment in the unit tests); the returned alignment is one of Biopython's co-optimal alignments on 2 050 / 2 050 rows with <= 40 co-optimal alignments (adjacent gap columns of opposite type are order-normalised). **parasail 1.3.4** `nw`/`sg_*` (integer rows, 1 680): equal on all 12 variants for the exact ones (`nw`, `sg_qb`, `sg_db`, `sg_qx`, `sg_dx`, `sg_qb_db`); the variants with a free *end* flag (`sg`, `sg_qe`, `sg_de`, `sg_qb_de`, `sg_qe_db`, `sg_qe_de`) are lower than the optimum on 23 rows and never higher — parasail does not consider the border end cells (m,0) / (0,n) (verified: an exhaustive reference that skips those two cells reproduces parasail on 300/300 random rows). **EMBOSS needle 6.6.0 default** (`-endweight N`, free end gaps on both sequences = `EndGaps.All`, 287 protein/DNA rows with fractional gaps via `MatrixScoring.Emboss`): equal on 277, lower (never higher) on 10 — the known embaln.c quirk that opens a last-row/column end gap only from the match matrix (F41). **Biopython quirk:** `PairwiseAligner.score()` disagrees with its own `align().score` on 25 rows (degenerate paths with one sequence entirely against gaps, e.g. `score('AC','T')` = −2 with `left_deletion_score = right_insertion_score = 0` while the alignment `AC-`/`--T` it returns scores 0); the oracle is `align().score`, which equals the exhaustive optimum on every row.

**Degenerate end paths.** A path may consume nothing of one sequence (zero-length overlap): `SemiGlobalAlign("AC", "T", d = −1, OverlapSuffix1Prefix2)` → `AC-`/`--T`, score 0, coordinates (2, 1, 0, −1).

### Substitution-matrix scoring (A1-2, B09 F41)

`SequenceAligner.SemiGlobalAlign` and `SemiGlobalAlignAffine` also accept a `MatrixScoring` (BLOSUM62 / PAM250 / EDNAFULL residue-pair scores, fractional gap scores; `string`, `ProteinSequence`, `DnaSequence`) and run the same dynamic-program core as the integer overloads; see [Global_Alignment_Needleman_Wunsch.md §5.4](Global_Alignment_Needleman_Wunsch.md) for the conversion, the EMBOSS gap convention (`MatrixScoring.Emboss`) and the oracle (Biopython end_insertion_score = 0 / parasail sg_dx; EMBOSS has no one-sided free-end mode, Biopython 1.88, parasail 1.3.4); the `EndGaps` overloads take the same `MatrixScoring` (needle's default = `EndGaps.All`).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null query or reference | Throws `ArgumentNullException`. | Public API guard. |
| Query embedded exactly inside the reference | Returns the embedded match score from the last-row maximum. | Validated by `SequenceAligner_SemiGlobalAlign_Tests`. |
| Query at the start or end of the reference | Leaves the unmatched reference prefix or suffix free. | This is the point of the fitting initialization and last-row maximum. |
| Query aligned entirely against gaps (e.g. all-mismatch under a harsh mismatch penalty) or empty query | The fitted interval is empty: `EndPosition2 = StartPosition2 - 1` (the path ends in the first maximal last-row column, `j = 0` when it ties). | Biopython enumerates the empty placement among several co-optimal ones; Seqeron reports the leftmost. |
| Tied fittings (query matches several places equally well) | One co-optimal interval is reported (first maximal last-row column, tie order diagonal > up > left). | Verified to be in the set of intervals of all Biopython 1.88 co-optimal alignments (`PairwiseAligner.align`). |
| All mismatches under negative mismatch penalties | Can return a negative score. | The recurrence has no zero floor. |
| Query longer than the reference | Still aligns the full query, possibly introducing gaps into the reference alignment. | The query is always globally represented in the fitting variant. |

### 6.2 Limitations

Like the global path, the method keeps a full `O(mn)` matrix up to 25 000 000 cells and switches to the exact linear-space program above it (`SemiGlobalAlignLinearSpace` / `SemiGlobalAlignAffineLinearSpace`, all 16 `EndGaps` policies, identical rows); `SemiGlobalAlignGeneralGap` (non-affine `w(k)`) keeps `O(mn)` memory and is limited to 4 000 000 cells. (`SemiGlobalAlign`) a linear gap model based on `GapExtend` (use `SemiGlobalAlignAffine` for affine costs, `SemiGlobalAlignGeneralGap` for an arbitrary gap-cost function). Every end is either free or an ordinary gap (no separate end-gap penalty value). Raw-string overloads exist for the linear and affine variants (uppercasing only, no alphabet validation).

## 8. References

1. [Sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment)
2. [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
3. [Rosalind: Finding a Motif with Modifications (SIMS)](https://rosalind.info/problems/sims/)
4. [Rosalind: Semiglobal Alignment (SMGB)](https://rosalind.info/problems/smgb/)
5. Gotoh, O. (1982). "An improved algorithm for matching biological sequences." J Mol Biol 162: 705-708.
6. Brudno, M. et al. (2003). "Glocal alignment: finding rearrangements during alignment." Bioinformatics 19 Suppl 1: i54-i62.
