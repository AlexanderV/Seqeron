# Evidence Artifact: ALIGN-GLOBAL-001

**Test Unit ID:** ALIGN-GLOBAL-001
**Algorithm:** Global Alignment (Needleman–Wunsch)
**Date Collected:** 2026-03-06

---

## Online Sources

### Wikipedia: Needleman–Wunsch Algorithm
**URL:** https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
**Accessed:** 2026-03-06 · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)

**Key Extracted Points:**

1. **Purpose:**
   - Algorithm for computing an optimal global alignment between two sequences.
   - One of the first applications of dynamic programming to biological sequences.
   - Aligns the entire sequences end-to-end.

2. **Scoring System (from "Choosing a scoring system"):**
   - Match: +1 (or other positive value depending on substitution matrix)
   - Mismatch: −1 (or negative value from substitution matrix)
   - Indel (gap): −1 (or specified gap penalty)
   - Quote: "The score of the whole alignment candidate is the sum of the scores of all the pairings."

3. **Linear Gap Penalty Model (from "Advanced presentation of algorithm"):**
   - Uses a similarity function $S(a,b)$ and a single linear gap penalty $d$.
   - Basis:
     - $F_{0,j} = d \cdot j$
     - $F_{i,0} = d \cdot i$
   - Recurrence:
     - $F_{i,j} = \max(F_{i-1,j-1} + S(A_i,B_j),\ F_{i,j-1} + d,\ F_{i-1,j} + d)$
   - The entry $F_{n,m}$ gives the maximum score among all possible alignments.

4. **Traceback Rules (from "Tracing arrows back to origin"):**
   - Diagonal path (top-left cell) → match or mismatch.
   - Horizontal path (left cell) → gap in the "side" sequence.
   - Vertical path (top cell) → gap in the "top" sequence.
   - "If there are multiple arrows to choose from, they represent a branching of the alignments."

5. **Complexity:**
   - Time: $O(nm)$ for sequences of length $n$ and $m$.
   - Space: $O(nm)$ for the full matrix.

6. **Example:**
   - Sequences: GCATGCG and GATTACA
   - Scoring: match = +1, mismatch/indel = −1
   - One possible optimal alignment: `GCATG-CG` / `G-ATTACA`
   - Score breakdown: `+−++−−+−` → $4 \times (+1) + 4 \times (-1) = 0$
   - Border initialization shown: first row is 0, −1, −2, −3, −4, −5, −6, −7.

7. **Gap Penalty Extensions (from "Gap penalty"):**
   - Affine gap penalties (gap-open + gap-extend) are mentioned as an *extension* to the basic model.
   - Quote: "The simple and common way to do this is via a large gap-start score for a new indel and a smaller gap-extension score for every letter which extends the indel."
   - The standard NW pseudocode uses **only the linear penalty** $d$.

### Wikipedia: Sequence Alignment (Global vs. Local)
**URL:** https://en.wikipedia.org/wiki/Sequence_alignment
**Accessed:** 2026-02-01 · re-fetched 2026-10-14 (SRC citation audit: see corrections below)

**Key Extracted Points:**

1. **Global Alignment:**
   - Attempts to align the entire length of all sequences (every residue).
   - Most useful when sequences are similar and of roughly equal size.
   - Can start and/or end in gaps.
   - Needleman–Wunsch is the canonical global alignment method.

2. **Multiple Optimal Alignments** (not on the Sequence alignment page — re-attributed 2026-10-14, SRC citation audit):
   - Multiple traceback paths can exist for the same optimal score. Stated on the Needleman–Wunsch page (above):
     "If there are multiple arrows to choose from, they represent a branching of the alignments. If two or more
     branches all belong to paths from the bottom right to the top left cell, they are equally viable alignments."
     (likewise the Smith–Waterman page: "An element can receive score from more than one element, each will form a
     different path if this element is traced back.")
   - Common in sequences with repeated motifs or low complexity regions — an editorial observation, not stated on
     any of the cited pages.

---

## Test Dataset: Wikipedia Example

| Parameter | Value |
|-----------|-------|
| Sequence 1 | GCATGCG |
| Sequence 2 | GATTACA |
| Match Score | +1 |
| Mismatch Score | −1 |
| Gap Penalty (d) | −1 |
| Optimal Score | 0 |
| One Optimal Alignment (seq1) | `GCATG-CG` |
| One Optimal Alignment (seq2) | `G-ATTACA` |
| Score Breakdown | 4 matches (+4) + 2 mismatches (−2) + 2 gaps (−2) = 0 |

### Wikipedia DP Matrix (Border Initialization)

Standard NW border for the example (d = −1):

| | - | G | A | T | T | A | C | A |
|---|---|---|---|---|---|---|---|---|
| **-** | 0 | −1 | −2 | −3 | −4 | −5 | −6 | −7 |
| **G** | −1 | | | | | | | |
| **C** | −2 | | | | | | | |
| **A** | −3 | | | | | | | |
| **T** | −4 | | | | | | | |
| **G** | −5 | | | | | | | |
| **C** | −6 | | | | | | | |
| **G** | −7 | | | | | | | |

Border values: $F(i,0) = -i$, $F(0,j) = -j$ — directly from Wikipedia.

---

## Affine gaps (added 2026-09 review)

- **Gotoh O (1982)** J Mol Biol 162:705-708 — three-matrix DP for w(k) = v + u·k. **Flouri T, Kobert K, Rognes T, Stamatakis A (2015)** "Are all global alignment algorithms and implementations correct?" bioRxiv 10.1101/031500 — Gotoh's P/Q border initialization is wrong (opened via search snippets; bioRxiv PDF blocked). Implemented as `GlobalAlignAffine` with M/X/Y states, X(i,0) = o + i·e, Y(0,j) = o + j·e, all other border cells −∞, state-switching traceback.
- **NCBI BLAST+ User Manual** (NBK279684, blastn options; via search snippet): blastn task default reward 2 / penalty −3 / gapopen 5 / gapextend 2 → `SequenceAligner.BlastDna`.
- Cross-check: Biopython 1.88 `PairwiseAligner` (`open_gap_score = GapOpen + GapExtend`) and parasail 1.3.4 `nw_trace` (`open = -(GapOpen + GapExtend)`) agree with `GlobalAlignAffine` on 2008 cases (random 1–40 nt, 8 scoring schemes, edge cases). A vs CAC (1, −10, −5, −1): −11 (Gotoh's original border: −9).

## Deviations and Assumptions

**None.**

- The implementation follows the standard Needleman–Wunsch linear gap penalty model exactly as given in the Wikipedia pseudocode.
- `ScoringMatrix.GapExtend` acts as the linear gap penalty $d$. `ScoringMatrix.GapOpen` is not used by `GlobalAlign`; it is used by `GlobalAlignAffine` (see above).
- When multiple optimal alignments exist, one is returned deterministically. This is explicitly allowed by the source.
- Empty input (2026-11, A1-1): an empty sequence is not special-cased by any overload; the NW border `F(0,j) = d·j`, `F(i,0) = d·i` (Needleman & Wunsch 1970; Durbin, Eddy, Krogh, Mitchison 1998 §2.3) gives the all-gap alignment (`""` vs `ACGT` → `----`/`ACGT`, score `4·d`; affine `o + 4·e`), and two empty sequences score 0 with empty rows. **No reference tool is an oracle for this case:** Biopython 1.88 `PairwiseAligner.score` raises `ValueError: sequence has zero length`, parasail 1.3.4 `nw_trace` raises `_s1Len must be > 0`, EMBOSS needle 6.6 rejects a zero-length sequence, scikit-bio 0.7.4 `pair_align` raises `Sequence 1 has a length of zero`, and `Bio.pairwise2` returns no alignment (all run 2026-11-18). The earlier audit remark "Biopython score = −8 for CCAC/\"\"" is therefore not reproducible in 1.88; the value follows from the recurrence (and equals the indel cost of the Hirschberg base case and of `ApproximateMatcher` weighted edit distance). Null input throws `ArgumentNullException` on every overload.

---

## References

1. Needleman, Saul B. & Wunsch, Christian D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology. 48(3): 443–53.
2. Wikipedia contributors. "Needleman–Wunsch algorithm". In: Wikipedia, The Free Encyclopedia. Accessed: 2026-03-06.
3. Wikipedia contributors. "Sequence alignment". In: Wikipedia, The Free Encyclopedia. Accessed: 2026-02-01.

---

## Substitution-matrix alignment with fractional gap scores (audit A1-2, B09 F41)

- **Opened/used:** EMBOSS 6.6 source `nucleus/embaln.c` (`embAlignPathCalcWithEndGapPenalties`, raw.githubusercontent.com) for the end-gap behaviour of needle; EMBOSS 6.6.0 `needle`/`water` binaries; Biopython 1.88 `Bio.Align.PairwiseAligner` + `substitution_matrices` (BLOSUM62, PAM250, NUC.4.4); parasail 1.3.4 (`nw`, `sw`, `sg_dx`, `blosum62`, `pam250`, `nuc44`).
- **Gap convention:** EMBOSS `-gapopen G -gapextend E` = cost G for a gap of length 1, G + (k-1)E for length k (verified: HBA/HBB score identical to Biopython with `open_gap_score = -G`, `extend_gap_score = -E`) -> `GapOpen = -(G - E)`, `GapExtend = -E`. `needle` default `-endweight N` leaves end gaps free on both sequences; the oracle uses `-endweight Y -endopen G -endextend E` (end gaps charged like the unit's global alignment).
- **Numbers (`scripts/alignment/generate_matrix_align_oracle.py`, 1 986 rows):** Biopython = Seqeron on 1 986/1 986 (1e-6); parasail = Seqeron on 616/616 integer-gap rows; EMBOSS needle/water = Seqeron on 1 103/1 121 (18 EMBOSS-internal disagreements, see the algorithm doc §5.4). HBA_HUMAN vs HBB_HUMAN, EBLOSUM62, gapopen 10 / gapextend 0.5: needle -endweight Y 292.5, water 293.5, linear gap 0.5: 366.5.
- **Deviations:** gap scores must be <= 0 and are resolved to 1e-6 (EMBOSS keeps float; Biopython compares with an epsilon); sequences must consist of symbols of the matrix alphabet (EMBOSS maps unknown characters to a zero-score code); the ends-free needle default is not implemented (A1-3).

## General gap-cost functions (audit A1-5, B09 F45)

**Sources actually used (2026-11-18):** Waterman, Smith & Beyer (1976) (WebSearch snippets: "Some biological sequence metrics", arbitrary gap weights, O(m n (m+n))), Gusfield 1997 §11.2, Biopython 1.88 `PairwiseAligner` with a callable `gap_score(position, length)` (selects the "Waterman-Smith-Beyer global/local alignment algorithm"; for a callable the end-gap properties have no effect — verified: identical scores with and without `end_insertion_score = 0`).

**Oracle** (`scripts/alignment/generate_general_gap_oracle.py`, 1 440 rows): 960 global/local pairs x 8 gap families (linear, affine, logarithmic, square-root, quadratic, capped, step, two-slope; scalar / BLOSUM62 / EDNAFULL): Biopython (`align().score` and `score()`) = independent three-state implementation on 960/960, = Seqeron general program and public API on 960/960 (1e-8); 480 ends-free rows vs the independent implementation. **Model fixed by the oracle:** every maximal gap run is charged w(k) once (single-run); the textbook single-matrix recurrence max_k H(i-k,j) + w(k) lets a gap be split into cheaper adjacent pieces and disagrees with Biopython on 12 % (log) / 57 % (quadratic) / 27 % (step) of 200 random global pairs per family.


## Exact linear-space alignment (audit A1-6, B09 F46)

- **Sources actually used (2026-11-18):** Hirschberg 1975, Myers & Miller 1988 (CABIOS 4:11, PMID 3382986), Huang & Miller 1991 (Adv Appl Math 12:337, DOI 10.1016/0196-8858(91)90017-D), Huang-Hardison-Miller 1990 / SIM (WebSearch records in 2026-11-18; on 2026-11-19 Hirschberg 1975 was opened in full and the Huang et al. abstracts - Myers-Miller and the Huang-Miller full texts stay unreachable, see "Primary sources opened (audit A1-22)" at the end of this file). Reference tools run locally: EMBOSS 6.6.0 `stretcher` (Myers-Miller linear-space global, end gaps penalised; gap of length k = gapopen + (k-1) gapextend, i.e. `GapOpen = -(G - E)`, `GapExtend = -E`) and `matcher` (SIM local, `-alternatives 1`), Biopython 1.88 `PairwiseAligner`, parasail 1.3.4.
- **Contract verified:** the linear-space result equals the full-matrix aligner row for row (score, rows, coordinates, type, ties): 4 750 random / tie-rich / long-indel / matrix cases with base cells 1-100 over all 16 end-gap policies, local, linear + affine; 2.6-3 kb pairs; the automatic switch on a 5 300 x 5 200 pair for all six aligners.
- **Biopython / parasail** (`scripts/alignment/generate_linear_space_oracle.py`, 1 280 pairs): scores equal on all rows; parasail equal on every row where a variant exists; returned alignments re-score independently and lie in Biopython's co-optimal set whenever <= 40 alignments exist.
- **EMBOSS stretcher / matcher oracle rows (160 + 160):** scores equal on all rows; EMBOSS rows re-score to their score under the Seqeron model; where Biopython finds exactly one optimal alignment (49 stretcher rows, 93 matcher rows) the EMBOSS and Seqeron rows (and matcher's 1-based coordinates) are equal; elsewhere ties are broken differently (55 / 160 and 99 / 160 identical rows).
- **Large EMBOSS cross-check (DNA, 10 % edits + 3 planted indels of 20-300 nt, EDNAFULL, gapopen 16 / gapextend 4):** `stretcher` score = `GlobalAlignAffine` and `matcher` score = `LocalAlignAffine` for 100, 500, 1 000, 2 000, 5 000, 20 000, 50 000 nt (14 / 14): 50 000 x 50 360: 191 903 / 192 525; 20 000 x 19 965: 75 454 / 75 477. Timing: stretcher 1.6 s / 10.8 s for 20 kb / 50 kb, this implementation 15.6 s / 89.5 s.
- **Memory (Release, process peak working set):** 20 000 x 19 994 affine global 111 MB (4.8 GB full matrices), local 118 MB, fitting 113 MB, linear global 57 MB; 100 000 x 100 029 (10^10 cells) affine global 185 MB / 391.7 s, linear global 73 MB / 153.9 s; 10 000 x 10 000 affine: 1 179 MB full matrix vs 103 MB.
- **Limits:** `LinearSpaceCellThreshold` 25 000 000 cells; `GeneralGapMaxCells` 4 000 000; `ProfileMatrixMaxCells` 50 000 000 (progressive / iterative / consistency MSA profile DP); see `docs/algorithms/Alignment/Linear_Space_Alignment.md` §4.


## Waterman-Eggert non-intersecting local alignments and co-optimal enumeration (audit A1-7, B09 F47)

- **Sources actually used (2026-11-18):** Waterman & Eggert 1987, J Mol Biol 197:723-728 (the paper was opened 2026-11-19, see "Primary sources opened (audit A1-22)" at the end of this file; the definition is also confirmed from the implementation of the tool that carries it, below); EMBOSS 6.6.0 `emboss/matcher.c` (raw.githubusercontent.com `kimrutherford/EMBOSS` mirror, opened: `matcher_Sim` "Smith Waterman Eggert local alignment", `matcher_BigPass`, `MATCHERDIAG`, `matcher_Diff`, `matcher_Locate`/`matcher_SmallPass`); Biopython 1.88 `Bio/Align/_pairwisealigner.c` (raw.githubusercontent.com `biopython-188`, opened: `SELECT_TRACE_SMITH_WATERMAN_*`, ENDPOINT / STARTPOINT handling, "we don't allow zero-score extensions to alignments"); installed `matcher`/`water` 6.6.0.
- **Definition confirmed from `matcher`:** `MATCHERDIAG(i, j, c, value)` assigns the diagonal value only when `(i, j)` was never used before (`row[i]` list of used columns); gap moves are unaffected; `matcher_Diff` marks `(A[I], B[J])` used for every REP (diagonal) column of the reported alignment. So the forbidden set = the aligned residue pairs of earlier alignments, gap columns do not forbid anything. The `water` program reports only the best alignment.
- **Tie-independent verification of matcher (`scripts/alignment/generate_topn_oracle.py` `matcher_scores`, run on 397 oracle rows = 1 695 alignments, and on a separate 300-case check = 1 207 alignments):** each printed alignment is parsed (rows + coordinates), re-scored with the EMBOSS gap convention (gapopen G, gapextend E = GapOpen -(G - E), GapExtend -E), checked to use no pair of an earlier alignment, and its score is compared with the maximum of the independently recomputed matrix in which the earlier pairs are forbidden: 0 deviations (assertion in the generator). Hence matcher implements exactly the Waterman-Eggert definition used here.
- **Seqeron vs matcher:** round-1 score equal on 397 / 397 rows; the whole score list equal on 390 / 397 (the 7 others are tie rows: a later round depends on which co-optimal alignment was taken earlier; matcher's alignment is still optimal for its own forbidden set).
- **Independent Waterman-Eggert reference** (Python, full recomputation of the matrix with the forbidden pairs after every alignment - the C# code recomputes only the affected cells): equal on all 1 350 oracle rows (score, 0-based coordinates, both gapped rows) and on 300 extra long pairs (up to 174 residues, up to 40 alignments each, 6 405 alignments; tandem-rich ACGT / AC, linear + affine, six scoring schemes).
- **Exhaustive definition check:** for 198 oracle pairs and 600 random C# pairs of <= 5 residues every reported alignment is the exhaustive maximum (all gapped column sequences enumerated) among alignments avoiding the earlier pairs, and the list ends exactly when no positive alignment is left.
- **Biopython 1.88 co-optimal set** (`PairwiseAligner(mode='local').align(a, b)`, scalar scoring + BLOSUM62 with gaps in quarters, linear + affine): `LocalAlignCoOptimal` / `LocalAlignAffineCoOptimal` return the identical set (score, coordinates, rows) on all 1 345 rows with <= 60 co-optimal alignments (5 rows with more are `NA`). Two Biopython rules had to be matched and are documented: an alignment ends with a diagonal column, and zero-score extensions of a shorter optimal alignment are not enumerated (`ENDPOINT` cells are not passed through).
- **Speed:** 1.9 kb x 1.3 kb tandem-rich DNA, 50 alignments: 1.2 s (affine) / 1.4 s (linear), including process start.


## Co-optimal global / ends-free alignments (audit A1-7b, B09 F59)

- **Reference:** Biopython 1.88 `PairwiseAligner(mode='global').align(a, b)` iterates every optimal alignment (rows rebuilt from the coordinates; `left/right_insertion/deletion_score = 0` for the free ends; affine `open_gap_score = GapOpen + GapExtend`, `extend_gap_score = GapExtend`). `scripts/alignment/generate_global_cooptimal_oracle.py` writes `global_cooptimal_oracle.tsv`: 1 920 pairs (14 scalar schemes + BLOSUM62 with fractional gaps; linear + affine; all 16 end-gap policies; tandem repeats, long gaps, zero gap costs), 61 rows above the 150-alignment cap.
- **Independent check inside the script:** exhaustive enumeration of every column sequence of pairs with <= 5 residues (score written from the definition: substitution columns + gap cost, free leading/trailing gap columns 0, affine run = same row) is asserted equal to Biopython's set on every tiny row.
- **Seqeron vs Biopython:** identical sets on all 1 859 comparable rows (935 with several alignments); every returned alignment re-scores to Biopython's optimum; exhaustive C# enumeration equal for all 16 policies, linear and affine; capped prefix of the 61 large rows distinct and optimal.
- **Finding during the resume of the work package:** with zero gap costs a free flank is also a zero-cost gap run of the path, so the same rows ended at several cells; the first implementation skipped the "dominated" end cells and then listed a different first alignment than `SemiGlobalAlign` on 1 oracle row. Replaced by a set of already-listed rows (the first alignment is always the single aligner's traceback) and coordinates computed from the rows (documented in the algorithm doc §5.6).

## Primary sources opened (audit A1-22, 2026-11-19)

Earlier entries above cite these papers from WebSearch snippets only; on 2026-11-19 (user-agent `Seqeron-review/1.0`) the following full texts were actually opened (PDF -> `pdftotext`; the Duke / USC copies are OCR scans) and the definitions used by the code were checked against them:

| Source (opened copy) | What it states | Checked against the code |
|---|---|---|
| Needleman & Wunsch 1970, J Mol Biol 48:443 (`courses.cs.duke.edu/spring21/compsci260/resources/AlignmentPapers/1970.needleman.wunsch.pdf`) | the original global comparison matrix | `GlobalAlign` (additive gap model; NW recurrence in the form of Wikipedia/Gusfield) |
| Waterman, Smith & Beyer 1976, Adv Math 20:367 (`.../1976.waterman.smith.beyer.pdf`) | arbitrary gap weights, matrix recurrence with a maximum over all gap lengths | `*GeneralGap` (F45) |
| Gotoh 1982, J Mol Biol 162:705 (`.../1982.gotoh.pdf`, 4 pages) | the WSB algorithm of M²N steps reduced to essentially MN steps for w_k = u·k + v (three recurrences D, P, Q) | `GlobalAlignAffine` / `LocalAlignAffine`: the same three states (here M, X, Y with a maximum instead of a minimum, `GapOpen + k·GapExtend`); the border initialisation is the Flouri et al. 2015 correction, unchanged |
| Hirschberg 1975, CACM 18:341 (`ics.uci.edu/~dhirschb/pubs/p341-hirschberg.pdf`, the author's page) | linear-space LCS by divide and conquer on the middle row | `LinearSolve` / `LinearMidpoint` (F46) |
| Waterman & Eggert 1987, J Mol Biol 197:723 (`dornsife.usc.edu/msw/wp-content/uploads/sites/236/2023/09/msw-079.pdf`, Waterman's page) | "Two alignments are said to intersect if they have a match or mismatch in common" and "the simplest idea is to recompute the matrix, simply not allowing matches or mismatches that were involved in alignments already output" (the paper then restricts the recomputation to the cells that can change) | `LocalAlignTopN` / `LocalAlignAffineTopN` (F47) and `LocalAlignTopNLinearSpace` (F77): forbidden aligned pair = diagonal column; EMBOSS `matcher` source (already opened) implements the same |
| Huang, Hardison & Miller 1990, CABIOS 6:373 and Huang & Miller 1991, Adv Appl Math 12:337 (abstracts only: `digitalcommons.mtu.edu/michigantech-p/9773`, `/5622`) | local similarity in space proportional to M + N; "k best non-intersecting local alignments ... O(M + N + K) space" (K = total length of the alignments) | F46 (locate-then-solve) and F77 (rounds of the linear program with the earlier pairs forbidden) |

**Not reachable (proof):** Myers & Miller 1988, CABIOS 4:11 (`academic.oup.com/bioinformatics/article-pdf/4/1/11/495701/4-1-11.pdf` HTTP 403; `cs.arizona.edu/~gene/PAPERS/linear.pdf` 404; a mirror was not found by WebSearch; Gusfield's UC Davis notes `cs.ucdavis.edu/~gusfield/cs224f11/spacenotes.pdf` were opened but only cite it); full texts of Huang & Miller 1991 / Huang-Hardison-Miller 1990 (`sciencedirect.com/science/article/pii/019688589190017D` 403; OUP 403; the Michigan Tech repository holds the abstract only). The Myers-Miller state-carrying midpoint is therefore verified against the EMBOSS `stretcher` binary and the full-matrix aligners (F46), not against the paper's text.

## 2026-11 audit round 2 (A2-37 / A2-43, B09 F80 / F85): scoring contract, `null` scoring overloads
- Defect values (exact integer arithmetic, no external oracle needed): `GlobalAlign("A"x50, "A"x50, new ScoringMatrix(1_000_000_000, -1, -2, -1)).Score` returned 1 999 999 904 (true score 50 x 1e9 = 5e10), `LocalAlign` 2 000 000 000, `GlobalAlignAffine` 1 999 999 900; `GlobalAlignAffine(.., new ScoringMatrix(1, -1, int.MinValue, int.MinValue))` raised `IndexOutOfRangeException` (negation overflow); `GlobalAlign("ACGT","AC", new ScoringMatrix(1,-1,5,5))` = 30 (a positive gap "score" rewards gaps; `MatrixScoring` always rejected it). The reported `AlignmentResult.Score` is a 32-bit integer, so a score beyond the range cannot be returned: the one honest contract is the existing matrix-path guard (`max|step| x (m + n + 2) < 2^28`, the affine -inf sentinel is `int.MinValue / 4`), now applied to every scoring model at the point where the dynamic program is about to run (`PairScoring.EnsureRange`).
- Gap scores are penalties (the convention already documented for both models in `ScoringMatrix` / `MatrixScoring`: EMBOSS `-gapopen` / `-gapextend` and BLAST+ `-gapopen` / `-gapextend` are positive penalties, i.e. scores <= 0 here; see the F1 evidence of this unit): the matrix path always rejected a positive gap, the scalar path did not. A local alignment (Smith-Waterman 1981) never starts with a gap only when gaps do not score positively, which the linear-space local program relies on (Linear_Space_Alignment.md). Hence `GapOpen`, `GapExtend` <= 0 for both models (no external source defines a meaningful alignment with rewarded gaps).
- `null` for the scoring argument: C# 13 `OverloadResolutionPriorityAttribute` (Microsoft Learn, "Overload resolution priority", .NET 9+; the repository targets net10.0): among applicable candidates only those of the highest priority remain, so `null` binds to the `ScoringMatrix?` overload (documented default `SimpleDna`) instead of failing with CS0121; calls with a `MatrixScoring` argument are unaffected (the `ScoringMatrix?` overload is not applicable).
