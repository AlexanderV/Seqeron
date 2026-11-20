# Linear-Space Pairwise Alignment (Hirschberg / Myers-Miller / Huang-Miller)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-GLOBAL-001, ALIGN-LOCAL-001, ALIGN-SEMI-001 (audit A1-6, F46) |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-11-18 |

## 1. Overview

The full-matrix pairwise aligners (`GlobalAlign`, `GlobalAlignAffine`, `LocalAlign`, `LocalAlignAffine`, `SemiGlobalAlign`, `SemiGlobalAlignAffine`, every ends-free policy, scalar and `MatrixScoring` overloads) need `O(m n)` memory: 4 bytes per cell for the linear gap model and 12 bytes for the affine (Gotoh) one — two 20 kb sequences need 4.8 GB, and `new int[m+1, n+1]` fails with `OutOfMemoryException` at about 46 k x 46 k. This unit adds the exact linear-space counterparts of all of them:

* **Automatic switch.** Above `SequenceAligner.LinearSpaceCellThreshold` = 25 000 000 cells (`(m+1)(n+1)`, e.g. two 5 kb sequences; 100 MB linear / 300 MB affine; the limit `AnchorBasedAligner` already used) every aligner above runs the linear-space program on its own. **The result is identical** — score, aligned rows, coordinates, alignment type, tie behaviour — so no caller or MCP output changes.
* **MCP (audit A3-17, decision recorded).** The MCP pairwise tools (`global_align`, `local_align`, `semi_global_align` and their `_affine` variants, `local_align_top_n`) run these programs through the automatic switch above; because the result is identical, a separate `*_linear_space` tool would only duplicate each tool's output under another name and is not exposed. `GlobalAlignAffineCoOptimal` (all co-optimal affine alignments) is exposed as `global_align_co_optimal` with `affine = true`; the co-optimal enumerators have no linear-space variant (their output is the set of all alignments, which is not O(m + n)).
* **Explicit API** (any size): `GlobalAlignLinearSpace`, `GlobalAlignAffineLinearSpace`, `LocalAlignLinearSpace`, `LocalAlignAffineLinearSpace`, `SemiGlobalAlignLinearSpace(…, EndGaps)`, `SemiGlobalAlignAffineLinearSpace(…, EndGaps)`, each for `string`, `DnaSequence` (scalar `ScoringMatrix`) and `string`, `ProteinSequence`, `DnaSequence` (`MatrixScoring`, fractional gaps), with an optional `maxFullMatrixCells` (size of the sub-problems solved with the full-matrix core; default `LinearSpaceBaseCells` = 2^22).

## 2. Scientific / Formal Basis

### 2.1 Sources

| Method | Reference | Role here |
|--------|-----------|-----------|
| Hirschberg | D. S. Hirschberg, "A linear space algorithm for computing maximal common subsequences", *Commun. ACM* 18:341-343 (1975) | divide and conquer on the middle row; one DP row alive |
| Myers-Miller | E. W. Myers, W. Miller, "Optimal alignments in linear space", *CABIOS* 4:11-17 (1988) — Gotoh's three-state affine recurrence under Hirschberg's scheme (EMBOSS `stretcher`) | affine gaps: the split vertex carries its state (M/X/Y), a gap running through the middle row is opened once |
| Huang-Miller | X. Huang, W. Miller, "A time-efficient, linear-space local similarity algorithm", *Adv. Appl. Math.* 12:337-357 (1991); Huang, Hardison, Miller, *CABIOS* 6:373 (1990) (SIM; EMBOSS `matcher`) | local: locate the optimal interval in `O(n)` memory, then solve it globally |

Sources actually used: WebSearch records of the Myers-Miller (PMID 3382986), Huang-Miller (DOI 10.1016/0196-8858(91)90017-D) and SIM descriptions; the EMBOSS 6.6.0 binaries `stretcher` (Myers-Miller global) and `matcher` (SIM local) run locally as oracles; Biopython 1.88 `PairwiseAligner` and parasail 1.3.4 for scores and co-optimal alignment sets. On 2026-11-19 Hirschberg 1975 was opened in full (author's page), the Huang-Hardison-Miller 1990 / Huang-Miller 1991 abstracts (Michigan Tech repository: "k best non-intersecting local alignments ... O(M + N + K) space") and Waterman-Eggert 1987 (Waterman's page); the Myers-Miller 1988 paper and the Huang-Miller full texts are unreachable (OUP / Elsevier HTTP 403, no mirror) - see `docs/Evidence/ALIGN-GLOBAL-001-Evidence.md` "Primary sources opened (audit A1-22)".

### 2.2 The program (two phases)

**Phase 1 — locate.** Global alignment (no free end): nothing to locate, the rectangle is the whole matrix `(0,0)–(m,n)`. Ends-free alignment and local alignment: one forward pass with two rolling rows finds the end cell (rules of the full-matrix core: first maximum of the free last row, then the first strictly larger maximum of the free last column; local: first row-major maximum, state M) and — carried through the same pass — the **start vertex** of the traceback of the full-matrix program (payload per cell and state = payload of the predecessor the traceback would choose; stop vertices: row 0 / column 0 when the start is free, the origin otherwise; local: the cell of the first aligned pair, i.e. the M cell whose predecessor value is 0).

**Phase 2 — solve the rectangle `[start, end]` globally** (Hirschberg / Myers-Miller). A sub-problem is "global alignment of `a[i0..i1)` against `b[j0..j1)` from arrival state `s0` at the origin to end state `s1`" (`s0 = M` and `s1` = best end state at the top level). For `m ≥ 2` rows and more than `maxFullMatrixCells` cells:

1. run the forward DP over rows `0..mid` (`mid = m/2`) with `O(n)` memory;
2. continue over rows `mid+1..m`; every cell (and affine state) carries the column — and state — at which the traceback of the full-matrix program leaves through row `mid` (`payload[cell] = payload[chosen predecessor]`); at row `mid` the payload is the cell itself;
3. the payload of the end cell is the split vertex `(mid, j*, s*)`; recurse on the top rectangle (`… → (mid, j*)`, end state `s*`) and the bottom rectangle (`(mid, j*) → …`, start state `s*`); the aligned rows are the concatenation;
4. sub-problems of at most `maxFullMatrixCells` cells (or one row) are solved by the full-matrix core with the forced start/end state.

Memory: `O(n)` per recursion level (freed before recursing) + the `O(m+n)` output; time: about `2 m n` cell updates (each level is half the area of the previous) — 1.3-2x the full-matrix fill.

### 2.3 Why the rows are identical to the full-matrix program (ties included)

The traceback of every full-matrix aligner is a deterministic walk from the end vertex that takes, at each vertex, the **first tight predecessor** in a fixed order (linear: diagonal > up > left; affine: M > X > Y). Let `opt` be the optimum. A predecessor edge is tight iff it lies on some optimal path (`F[p] + w = F[v]` and `F[v] + G[v] = opt` imply `F[p] + G[p] = opt`, where `G` is the best score from a vertex to the end). Hence the returned path `P` is the greedy walk through the DAG `D` of optimal edges.

*Lemma.* Let `v` be a vertex of `P` and `D_v` the DAG of optimal edges of the sub-problem "start at `v`, end at the original end". Then (i) every edge of `P` after `v` is in `D_v` (the tail of `P` is optimal for the sub-problem: a better tail would beat the optimum), and (ii) `D_v ⊆ D` (an optimal sub-path prefixed by an optimal path to `v` is an optimal path). Greedy-first among a subset that contains the greedy choice of the superset is the same choice, so the greedy walk of the sub-problem is the tail of `P`. The same holds for the head (same DP values). Therefore splitting at any vertex of `P` and solving the halves with the same tie rules reproduces `P`; the carried payload is exactly that vertex (the first vertex of `P` on row `mid` seen from the end).

For the local program the sub-problem starts after the first aligned pair (a local alignment never starts with a gap when gaps score ≤ 0; every scoring model - scalar and matrix - rejects a positive gap score with an `ArgumentException` since audit A2-37, so the premise always holds). For the ends-free programs the start vertex is the traceback's stop vertex, and the free leading/trailing overhangs are emitted against gaps exactly like the full-matrix traceback.

The statement is also what the tests check empirically: 6 000+ random, tie-rich (alphabets `A`, `AC`, `ACGT`, zero gap scores) and long-indel pairs with forced recursion down to 1-cell base cases, **all** 16 end-gap policies x linear/affine x local x matrix scoring: 0 differences in score, rows, coordinates or type.

### 2.4 Waterman-Eggert top-N in linear space (A1-7c, B09 F77)

`LocalAlignTopN` / `LocalAlignAffineTopN` (Waterman-Eggert 1987 non-intersecting local alignments, `Local_Alignment_Smith_Waterman.md` §5.5) keep the matrices up to `TopNMaxCells` cells. Above it, and for `LocalAlignTopNLinearSpace` / `LocalAlignAffineTopNLinearSpace` at any size, each round is one run of the local linear-space program of §2.2 on the matrix in which the aligned residue pairs of the earlier rounds are *forbidden* (the diagonal / M-state move into that cell is removed, gap moves are unchanged) - the recomputation scheme of Huang-Miller SIM (EMBOSS `matcher`). The predicate `ForbiddenPairs` (1-based absolute `(row, column)`, sorted columns per row, cursor `ForbiddenRow` so a DP row without forbidden pairs costs one dictionary lookup) is threaded through

* the locate passes `LocateLinearLocal` / `LocateAffine` (a forbidden cell has no diagonal candidate, so the first row-major maximum, its traceback tie order and its first aligned pair are those of the matrix the in-memory declumper holds after the same rounds);
* the interval solves `LinearSolve` / `AffineSolve`, their midpoint passes `LinearMidpoint` / `AffineMidpoint` and the full-matrix base cases (`LinearAlignCore` / `AffineAlignCore` with `LinearFillMatrix` / `AffineFill` / `Traceback`); sub-problems know their row / column offset so the predicate is always addressed by absolute coordinates.

*Why the rows are identical.* Removing diagonal edges only deletes edges of the DAG of §2.3; the greedy walk of the traceback, the Lemma (tail of the path is optimal for the sub-problem, `D_v ⊆ D`) and the "first tight predecessor" argument do not depend on which edges exist. The sub-problem after the first aligned pair never contains that (non-forbidden) pair, and a forbidden pair has value `-inf` in M, so no tight predecessor edge enters it. Memory is `O(m + n)` per round plus the forbidden pairs (the total length of the reported alignments); time is `(N + 1)` local-program runs for `N` reported alignments (the final run finds no positive alignment, or none >= `minScore`).

Tests (`SequenceAligner_TopNLinearSpace_Tests`): 480 random tie-rich pairs (320 scalar + 160 matrix) (alphabets `AC`/`ACG`/`ACGT`, planted repeats, 5 scalar schemes + BLOSUM62 / EDNAFULL with fractional gaps, k up to 11, minScore; base-case sizes 1, 9, 40, 200 so the divide and conquer recurses on tiny rectangles) forced through the linear path = the in-memory declumper row for row; all 1 350 rows of `TestData/Alignment/topn_oracle.tsv` (independent full-recompute reference, EMBOSS matcher 6.6.0 score lists on 397 rows) through the linear path; large inputs above `TopNMaxCells` (5 300 x 5 300 linear, 5 100 x 5 100 affine, planted similar segments) = the in-memory declumper (forced with `fullMatrixCells = long.MaxValue`), MCP `local_align_top_n` on 5 100 x 5 100.

## 3. Contract

* Same inputs, validation, case folding, empty-input and `null` contract as the full-matrix aligner of the family (see the Global / Local / Semi-Global documents).
* `*LinearSpace` with `ScoringMatrix`: `GapOpen`/`GapExtend` as the full-matrix family (linear: only `GapExtend`); `MatrixScoring`: residues validated against the matrix alphabet, fractional gap scores (exact integer model).
* Every aligner requires non-positive gap scores (`ArgumentException` otherwise, scalar `ScoringMatrix` and `MatrixScoring` alike; audit A2-37), which the local linear-space program relies on.
* Cancellation: `GlobalAlign(…, CancellationToken, IProgress)` honours the token in the linear-space program too (checked every 64 rows); progress is reported 1.0 at the end.
* Integer range: one check for both models (`PairScoring.EnsureRange`: `max|step| · (m + n + 2) < 2^28`, `ArgumentException` "Score range too large for the integer dynamic program"; the scalar model used to wrap around silently, audit A2-37).

## 4. Limits and memory guard

| Aligner | Memory | Limit |
|---------|--------|-------|
| every linear / affine / local / ends-free aligner | full matrix up to 25 M cells (100 MB / 300 MB), then `O(m + n)` | none beyond the integer score range and the `int` string length (~2^31 residues); time `O(m n)` |
| `LocalAlignTopN` / `LocalAlignAffineTopN` (all overloads) | in-memory declumper up to `TopNMaxCells` = 25 M cells, then `O(m + n)` per round (§2.4) | none; time `O(m n)` per reported alignment (+1) |
| `LocalAlignCoOptimal` / `LocalAlignAffineCoOptimal` | `O(m n)` matrices | 25 M cells: `ArgumentException` above (the enumeration is output-exponential, no linear-space form) |
| `*GeneralGap` (arbitrary `w(k)`, linear/affine `w` are delegated) | `O(m n)`, 33 bytes per cell, time `O(m n (m + n))` | `SequenceAligner.GeneralGapMaxCells` = 4 000 000 cells: larger inputs throw `ArgumentException` (no linear-space variant exists for a non-affine `w`) |
| progressive / iterative / consistency MSA profile DP | `O(m n)` doubles | `SequenceAligner.ProfileMatrixMaxCells` = 50 000 000 cells: `ArgumentException`; the exact center-star `MultipleAlign` uses the linear-space pairwise stage and has no such matrix |

Measured (Release build, one core, random DNA with 10 % edits, BLAST DNA scoring 2/-3/5/2, 2026-11-18; peak working set of the process):

| Pair | Aligner | Time | Peak working set | Full-matrix memory |
|------|---------|------|------------------|--------------------|
| 20 000 x 19 994 | `GlobalAlignAffine` (auto) | 16.0 s | 111 MB | 4.8 GB |
| 20 000 x 19 994 | `LocalAlignAffine` (auto) | 27.8 s | 118 MB | 4.8 GB |
| 20 000 x 19 994 | `SemiGlobalAlignAffine` (auto) | 28.3 s | 113 MB | 4.8 GB |
| 20 000 x 19 994 | `GlobalAlign` (linear, auto) | 6.0 s | 57 MB | 1.6 GB |
| 100 000 x 100 000 | `GlobalAlignAffine` (auto) | 391.7 s | 185 MB | 120 GB |
| 100 000 x 100 000 | `GlobalAlign` (linear, auto) | 153.9 s | 73 MB | 40 GB |

Speed against the full matrix (Release, 10 000 x 10 000 = 100 M cells, 1.2 GB full matrix): `GlobalAlignAffine` 2.6 s full / 3.3 s linear-space; `LocalAlignAffine` 3.1 s / 6.4 s (locate pass + solve); `GlobalAlign` (linear) 1.0 s / 1.5 s. EMBOSS 6.6.0 `stretcher` (hand-tuned C) needs 1.6 s for 20 000 x 20 000 and 10.8 s for 50 000 x 50 000 (this implementation 15.6 s and 89.5 s: about 8x slower, same score, exact tie behaviour of the full-matrix aligner instead of Myers-Miller's).

## 5. Verification

`SequenceAligner_LinearSpace_Tests` (20 tests):

* differential against the full-matrix core, row for row: 1 530 + 2 800 + 420 random/tie-rich/matrix cases with base cells 1…100 (forced deep recursion), all 16 `EndGaps` policies, linear + affine, scalar + BLOSUM62/EDNAFULL with fractional gaps; 2.6-3 kb related pairs with planted indels at the default base and with 50 000-cell bases; the automatic switch on a 5 300 x 5 200 pair (27.6 M cells) for all six aligners and protein BLOSUM62 above the threshold;
* memory: allocation of `GlobalAlignAffineLinearSpace` on a 5 600 x 5 600 pair is below 1/8 of the 12-byte-per-cell matrices; cancellation of an 8 000 x 8 000 linear alignment;
* oracle (`TestData/Alignment/linear_space_oracle.tsv`, `scripts/alignment/generate_linear_space_oracle.py`): 1 280 Biopython 1.88 / parasail 1.3.4 rows (scores of all six families incl. ends-free policies, long indels, ties; the returned alignment re-scores independently to the score and is one of Biopython's co-optimal alignments when ≤ 40 exist); 160 EMBOSS `stretcher` 6.6.0 rows (score equal, EMBOSS rows optimal, rows equal wherever Biopython finds a unique optimum) and 160 EMBOSS `matcher` 6.6.0 rows (score equal, rows and 1-based coordinates equal wherever the optimum is unique).

See `docs/Evidence/ALIGN-GLOBAL-001-Evidence.md` for the numbers (including the large `stretcher`/`matcher` score comparisons).

## 6. Tie behaviour

Identical to the full-matrix aligner (linear: diagonal > up > left; affine: M > X > Y; first maximal end cell). EMBOSS `stretcher`/`matcher` break ties differently: where several optimal alignments exist their rows differ (55 / 160 and 99 / 160 identical rows in the oracle sets, all of them where the optimum is unique, always the same score).

## 7. References

1. Hirschberg DS (1975) A linear space algorithm for computing maximal common subsequences. Commun ACM 18:341-343.
2. Myers EW, Miller W (1988) Optimal alignments in linear space. CABIOS 4:11-17 (PMID 3382986).
3. Huang X, Miller W (1991) A time-efficient, linear-space local similarity algorithm. Adv Appl Math 12:337-357.
4. Huang X, Hardison RC, Miller W (1990) A space-economical algorithm for computing local similarity. CABIOS 6:373-381.
5. Gotoh O (1982) An improved algorithm for matching biological sequences. J Mol Biol 162:705-708.
6. EMBOSS 6.6.0 `stretcher` (Myers-Miller), `matcher` (SIM) — local binaries; Biopython 1.88 `PairwiseAligner`; parasail 1.3.4.
