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

- **Sources actually used (2026-11-18):** Hirschberg 1975, Myers & Miller 1988 (CABIOS 4:11, PMID 3382986), Huang & Miller 1991 (Adv Appl Math 12:337, DOI 10.1016/0196-8858(91)90017-D), Huang-Hardison-Miller 1990 / SIM (WebSearch records only; papers not opened). Reference tools run locally: EMBOSS 6.6.0 `stretcher` (Myers-Miller linear-space global, end gaps penalised; gap of length k = gapopen + (k-1) gapextend, i.e. `GapOpen = -(G - E)`, `GapExtend = -E`) and `matcher` (SIM local, `-alternatives 1`), Biopython 1.88 `PairwiseAligner`, parasail 1.3.4.
- **Contract verified:** the linear-space result equals the full-matrix aligner row for row (score, rows, coordinates, type, ties): 4 750 random / tie-rich / long-indel / matrix cases with base cells 1-100 over all 16 end-gap policies, local, linear + affine; 2.6-3 kb pairs; the automatic switch on a 5 300 x 5 200 pair for all six aligners.
- **Biopython / parasail** (`scripts/alignment/generate_linear_space_oracle.py`, 1 280 pairs): scores equal on all rows; parasail equal on every row where a variant exists; returned alignments re-score independently and lie in Biopython's co-optimal set whenever <= 40 alignments exist.
- **EMBOSS stretcher / matcher oracle rows (160 + 160):** scores equal on all rows; EMBOSS rows re-score to their score under the Seqeron model; where Biopython finds exactly one optimal alignment (49 stretcher rows, 93 matcher rows) the EMBOSS and Seqeron rows (and matcher's 1-based coordinates) are equal; elsewhere ties are broken differently (55 / 160 and 99 / 160 identical rows).
- **Large EMBOSS cross-check (DNA, 10 % edits + 3 planted indels of 20-300 nt, EDNAFULL, gapopen 16 / gapextend 4):** `stretcher` score = `GlobalAlignAffine` and `matcher` score = `LocalAlignAffine` for 100, 500, 1 000, 2 000, 5 000, 20 000, 50 000 nt (14 / 14): 50 000 x 50 360: 191 903 / 192 525; 20 000 x 19 965: 75 454 / 75 477. Timing: stretcher 1.6 s / 10.8 s for 20 kb / 50 kb, this implementation 15.6 s / 89.5 s.
- **Memory (Release, process peak working set):** 20 000 x 19 994 affine global 111 MB (4.8 GB full matrices), local 118 MB, fitting 113 MB, linear global 57 MB; 100 000 x 100 029 (10^10 cells) affine global 185 MB / 391.7 s, linear global 73 MB / 153.9 s; 10 000 x 10 000 affine: 1 179 MB full matrix vs 103 MB.
- **Limits:** `LinearSpaceCellThreshold` 25 000 000 cells; `GeneralGapMaxCells` 4 000 000; `ProfileMatrixMaxCells` 50 000 000 (progressive / iterative / consistency MSA profile DP); see `docs/algorithms/Alignment/Linear_Space_Alignment.md` §4.


## Waterman-Eggert non-intersecting local alignments and co-optimal enumeration (audit A1-7, B09 F47)

- **Sources actually used (2026-11-18):** Waterman & Eggert 1987, J Mol Biol 197:723-728 (WebSearch record; the paper itself was not opened - the definition is taken from the implementation of the tool that carries it, below); EMBOSS 6.6.0 `emboss/matcher.c` (raw.githubusercontent.com `kimrutherford/EMBOSS` mirror, opened: `matcher_Sim` "Smith Waterman Eggert local alignment", `matcher_BigPass`, `MATCHERDIAG`, `matcher_Diff`, `matcher_Locate`/`matcher_SmallPass`); Biopython 1.88 `Bio/Align/_pairwisealigner.c` (raw.githubusercontent.com `biopython-188`, opened: `SELECT_TRACE_SMITH_WATERMAN_*`, ENDPOINT / STARTPOINT handling, "we don't allow zero-score extensions to alignments"); installed `matcher`/`water` 6.6.0.
- **Definition confirmed from `matcher`:** `MATCHERDIAG(i, j, c, value)` assigns the diagonal value only when `(i, j)` was never used before (`row[i]` list of used columns); gap moves are unaffected; `matcher_Diff` marks `(A[I], B[J])` used for every REP (diagonal) column of the reported alignment. So the forbidden set = the aligned residue pairs of earlier alignments, gap columns do not forbid anything. The `water` program reports only the best alignment.
- **Tie-independent verification of matcher (`scripts/alignment/generate_topn_oracle.py` `matcher_scores`, run on 397 oracle rows = 1 695 alignments, and on a separate 300-case check = 1 207 alignments):** each printed alignment is parsed (rows + coordinates), re-scored with the EMBOSS gap convention (gapopen G, gapextend E = GapOpen -(G - E), GapExtend -E), checked to use no pair of an earlier alignment, and its score is compared with the maximum of the independently recomputed matrix in which the earlier pairs are forbidden: 0 deviations (assertion in the generator). Hence matcher implements exactly the Waterman-Eggert definition used here.
- **Seqeron vs matcher:** round-1 score equal on 397 / 397 rows; the whole score list equal on 390 / 397 (the 7 others are tie rows: a later round depends on which co-optimal alignment was taken earlier; matcher's alignment is still optimal for its own forbidden set).
- **Independent Waterman-Eggert reference** (Python, full recomputation of the matrix with the forbidden pairs after every alignment - the C# code recomputes only the affected cells): equal on all 1 350 oracle rows (score, 0-based coordinates, both gapped rows) and on 300 extra long pairs (up to 174 residues, up to 40 alignments each, 6 405 alignments; tandem-rich ACGT / AC, linear + affine, six scoring schemes).
- **Exhaustive definition check:** for 198 oracle pairs and 600 random C# pairs of <= 5 residues every reported alignment is the exhaustive maximum (all gapped column sequences enumerated) among alignments avoiding the earlier pairs, and the list ends exactly when no positive alignment is left.
- **Biopython 1.88 co-optimal set** (`PairwiseAligner(mode='local').align(a, b)`, scalar scoring + BLOSUM62 with gaps in quarters, linear + affine): `LocalAlignCoOptimal` / `LocalAlignAffineCoOptimal` return the identical set (score, coordinates, rows) on all 1 345 rows with <= 60 co-optimal alignments (5 rows with more are `NA`). Two Biopython rules had to be matched and are documented: an alignment ends with a diagonal column, and zero-score extensions of a shorter optimal alignment are not enumerated (`ENDPOINT` cells are not passed through).
- **Speed:** 1.9 kb x 1.3 kb tandem-rich DNA, 50 alignments: 1.2 s (affine) / 1.4 s (linear), including process start.
