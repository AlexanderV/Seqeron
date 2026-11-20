# Evidence: ALIGN-LOCAL-001

**Test Unit ID:** ALIGN-LOCAL-001
**Algorithm:** Local Alignment (Smith-Waterman, Smith-Waterman-Gotoh, Waterman-Eggert top-N)
**Date Collected:** 2026-03-06 (Wikipedia pass); 2026-11-18 / 2026-11-19 (B09 audit round 1: primary sources, oracles)
**Last Audited:** 2026-11-19

This file was missing until audit item A1-8 (the other ten units of batch B09 had one); it collects the sources and the oracle numbers that, before, lived only in `docs/Validation/review-2026-09/B09.md` (F2, F39, F41, F45-F47, F77) and in the TestSpec.

---

## 1. Authoritative Sources

| Source | Location | Opened |
|--------|----------|--------|
| Smith TF, Waterman MS (1981) "Identification of common molecular subsequences", J Mol Biol 147:195-197 | Duke course copy `courses.cs.duke.edu/spring21/compsci260/resources/AlignmentPapers/` (`www.cs.duke.edu/courses/fall20/compsci260/resources/AlignmentPapers/1981.smith.waterman.pdf`, 3-page OCR scan; also Waterman's page `dornsife.usc.edu/msw/.../msw-042.pdf`) | 2026-11-19, full text |
| Gotoh O (1982) "An improved algorithm for matching biological sequences", J Mol Biol 162:705-708 | `courses.cs.duke.edu/spring21/compsci260/resources/AlignmentPapers/1982.gotoh.pdf` (4 pages) | 2026-11-19, full text |
| Waterman MS, Eggert M (1987) "A new algorithm for best subsequence alignments with application to tRNA-rRNA comparisons", J Mol Biol 197:723-728 | `dornsife.usc.edu/msw/wp-content/uploads/sites/236/2023/09/msw-079.pdf` (Waterman's page) | 2026-11-19, full text |
| Huang X, Hardison RC, Miller W (1990) CABIOS 6:373; Huang X, Miller W (1991) Adv Appl Math 12:337 | abstracts only (`digitalcommons.mtu.edu/michigantech-p/9773`, `/5622`); full texts: Elsevier / OUP HTTP 403 | 2026-11-19, abstracts |
| Flouri T, Kobert K, Rognes T, Stamatakis A (2015) bioRxiv 10.1101/031500 (border initialisation of Gotoh's recurrences) | WebSearch record (bioRxiv PDF blocked) | 2026-11-18, snippets |
| Wikipedia: Smith-Waterman algorithm | https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | 2026-03-06 |
| EMBOSS 6.6.0 `emboss/matcher.c` (`matcher_Sim`, `MATCHERDIAG`, `matcher_Diff`), installed `matcher` / `water` 6.6.0 | `raw.githubusercontent.com/kimrutherford/EMBOSS` mirror, local binaries | 2026-11-18 |
| Biopython 1.88 `Bio/Align/_pairwisealigner.c`, `PairwiseAligner(mode='local')`; parasail 1.3.4 `sw`, `sw_trace` | `raw.githubusercontent.com`, installed packages | 2026-11-18 |

---

## 2. Algorithm Definition (from the sources)

1. **Smith-Waterman 1981, equation (1)** (read in the paper): `H(i,j) = max{ H(i-1,j-1) + s(a_i,b_j), max_k { H(i-k,j) - W_k }, max_l { H(i,j-l) - W_l }, 0 }`, first row and column 0; "the pair of segments with maximum similarity is found by first locating the maximum element of H; the other matrix elements leading to this maximum value are then sequentially determined ... until an element of H equal to zero". The zero "is included to prevent calculated negative similarity". `LocalAlign` is the case `W_k = k·d` (linear gaps), `LocalAlignGeneralGap` the general `W_k` (Waterman-Smith-Beyer).
2. **Gotoh 1982** (read in the paper): for gap weights `w_k = u·k + v` the M²N-step algorithm of Waterman, Smith & Beyer (1976) is done in essentially MN steps with three recurrences (D, P, Q). `LocalAlignAffine` = the zero floor of (1) on Gotoh's three states (M, X, Y), a gap of length k scoring `GapOpen + k·GapExtend`; the border initialisation of the global case (Flouri et al. 2015) does not enter a local alignment (all border states unreachable).
3. **Waterman-Eggert 1987** (read in the paper): "Two alignments are said to intersect if they have a match or mismatch in common"; the simplest idea is "to recompute the matrix, simply not allowing matches or mismatches that were involved in alignments already output", the paper then recomputes only the cells that can change. `LocalAlignTopN` / `LocalAlignAffineTopN` implement exactly this: the forbidden set is the aligned residue pairs (diagonal columns) of the earlier alignments; EMBOSS `matcher`'s `MATCHERDIAG` implements the same definition (verified round by round on 1 695 of its alignments).
4. **Huang-Miller 1991 / Huang-Hardison-Miller 1990 (abstracts)**: local similarity in space proportional to M + N; "k best non-intersecting local alignments ... O(M + N + K) space". `LocalAlignLinearSpace` (locate interval, then global solve) and `LocalAlignTopNLinearSpace` (one linear-space run per reported alignment with the earlier pairs forbidden) follow this scheme.

The reported alignment ends at the first maximal cell in row-major order and the traceback (ties: diagonal > up > left, M > X > Y) stops at the first zero predecessor; this is the one fixed choice among co-optimal alignments (always a member of Biopython's co-optimal set).

---

## 3. Oracle numbers (all recorded in B09.md)

| Check | Result |
|-------|--------|
| F2 affine local vs Biopython 1.88 `mode='local'` and parasail 1.3.4 `sw_trace` | 3 192 cases (1 596 linear + 1 596 affine, 8 schemes): scores 100 % equal, alignments re-score, the end cell lies in Biopython's co-optimal set |
| Independent re-verification of F1-F3 (auditor harness) | 3 x 3 936 rows (random DNA, 2 / 4-letter alphabets, empty, single character, ties, gaps up to 40, lower case, 10 schemes) vs Biopython / parasail: 0 mismatches |
| F39 contract of all overloads | 4 090-row differential, 0 mismatches (empty = tagged-Local empty result; Biopython 1.88 / parasail / needle / scikit-bio reject empty input) |
| F41 `MatrixScoring` (BLOSUM62 / PAM250 / EDNAFULL, fractional gaps) | Biopython 1 986 / 1 986, parasail 616 / 616, EMBOSS needle / water 1 103 / 1 121 (the program's end-gap quirk) |
| F45 general gap costs | Biopython callable `gap_score` 960 / 960 rows (8 gap families) + exhaustive enumeration |
| F46 linear-space local | identical to the full matrix on 4 750 forced-recursion cases (all 16 end-gap policies for the family, local linear + affine) |
| F47 top-N (`TestData/Alignment/topn_oracle.tsv`, `scripts/alignment/generate_topn_oracle.py`) | 1 350 rows: independent full-recompute Waterman-Eggert identical (also 300 long pairs, 6 405 alignments); exhaustive enumeration of 198 + 600 tiny pairs; EMBOSS matcher 6.6.0 `-alternatives k`: 397 rows / 1 695 alignments verified optimal for the definition, score lists equal on 390 / 397 (7 tie divergences); Biopython co-optimal set identical on 1 345 rows |
| F77 top-N in linear space | 480 random pairs forced through the linear path = in-memory declumper; all 1 350 oracle rows equal; two 5 kb pairs above 25 M cells = in-memory declumper; memory 40-46 MB peak at 5-20 kb (matrix 95 MB - 1.5 GB) |

---

## 4. Documented corner cases

- No similarity (all scores <= 0): score 0, empty alignment, `Local`, positions -1; the same for an empty side (every overload).
- The traceback of a local alignment starts with a diagonal column (a local alignment never starts with a gap when gap scores are <= 0; the linear-space programs and the top-N rounds rely on this and require non-positive gap scores).
- Waterman-Eggert: two reported alignments may cross a matrix cell through a gap column but never align the same residue pair; ties between co-optimal alignments change later rounds (documented; every reported score is the true maximum of the remaining problem).
- `LocalAlignCoOptimal` follows Biopython's rule "zero-score extensions are not allowed".
- Wikipedia worked example: `TGTTACGG` / `GGTTGACTA`, +3 / -3 / -2: score 13, `GTT-AC` / `GTTGAC`.

## 5. Not reachable (proof)

Full texts of Huang-Hardison-Miller 1990 and Huang-Miller 1991 (Elsevier `sciencedirect.com/science/article/pii/019688589190017D` HTTP 403, OUP 403; repository records hold abstracts only) and Myers-Miller 1988 (OUP 403, author page 404); the linear-space local program is verified against the full-matrix aligner, EMBOSS `matcher` and the oracles above instead of the papers' text.

## Audit A2-28c (B09 F124): supermatcher

Sources actually opened: EMBOSS 6.6 `supermatcher.c`, `embword.c` (`embWordGetTable`, `embWordRabinKarpInit`, `embWordRabinKarpSearch`, `embWordMatchFirstMax`), `embaln.c` (`embAlignPathCalcSWFast`, `embAlignWalkSWMatrixFast`, `embAlignReportLocal`), `ajmem.h` (`AJCSET0` zeroes elements), `ajdefine.h`; the installed supermatcher 6.6.0. Oracle numbers: 260 / 260 runs (output files and error files).

## Audit A3-14b / A3-14c (B09 F143): water's border rule; why 7 of 397 matcher score lists differ

**water.** Sources actually opened: EMBOSS 6.6.0 `nucleus/embaln.c` `embAlignPathCalcSW` (text only). The running maximum `ret` is updated only inside `while (xpos != lenb) { ... while (ypos < lena) { ... if (mscore > ret) ret = mscore; ...` with `xpos`, `ypos` starting at 1, i.e. over the cells whose row AND column index are at least 1; the first row and column are filled by separate loops that never touch `ret`, and `ret` is taken before the floor at zero (all-negative matrices give a negative score). So an alignment that ENDS in the first row or column is never reported. `EmbossReferenceDp.Water` (score only, written from that text) gives the printed score of **639 of 639** oracle rows with dyadic penalties (of 860; the rows are compared exactly on all 860), and with the border cells added to the maximum it gives the textbook local optimum on all 250 sample rows. Of those 250 sample rows exactly **7** (ids 23, 51, 69, 240, 292, 305, 405) differ from the textbook optimum, each below it: the optimum ends in the first row / column (23, 305, 405) or the matrix has no positive cell (51, 69, 240, 292: the program reports the negative inner maximum).

**matcher `-alternatives` (F47).** F47's statement "= EMBOSS matcher" holds for round 1 (397 of 397 rows) and for the whole score list on **390 of 397** rows. The generator (`generate_topn_oracle.py`) proved that every matcher round is a legal Waterman-Eggert round (its score is the optimum with the pairs of matcher's OWN earlier rounds forbidden), so two lists can only differ through a **tie**: two different alignments of the same optimal score in an earlier round. `scripts/alignment/analyze_topn_matcher_ties.py` re-runs `matcher` 6.6.0 on every row and compares the aligned pairs of every round with the independent full-recompute reference (which `LocalAlignTopN` equals on every row): the 7 differing rows are lines 46, 487, 502, 787, 793, 865, 1348 of `topn_oracle.tsv`; in each the first alignment that differs has the SAME score as ours and comes no later than the first score difference (rounds 2, 2, 1, 2, 1, 2, 1; the score lists first differ in rounds 6, 4, 7, 3, 6, 3, 2). 195 of the 390 equal rows also break a tie differently without changing a score. Run output: "rows with a matcher list: 397; lists equal: 390 (of which 195 broke a tie differently without changing a score); lists differ: 7, explained by a tie broken differently in an earlier or the same round: 7". The tests (`SequenceAligner_TopN_Tests`, `SequenceAligner_TopNLinearSpace_Tests`) now assert 397 / 390 and the exact line list. Matching SIM's tie order needs the SIM program (`matcher_Sim`, `matcher_BigPass`, `matcher_Diff`, `matcher_Findmax`): the matcher port, audit A2-28d (DOABLE, carried).

## Audit A2-28d (B09 F145): `matcher` port, SIM's own tie order

Sources actually opened: EMBOSS 6.6.0 `emboss/matcher.c` in full (`raw.githubusercontent.com/kimrutherford/EMBOSS`, User-Agent `Seqeron-review/1.0`): `matcher_Sim`, `matcher_BigPass`, `matcher_Locate`, `matcher_SmallPass`, `matcher_Addnode`, `matcher_NoCross`, `matcher_Findmax`, `matcher_Diff`, `matcher_Calcons` and the macros `MATCHERORDER`, `MATCHERDIAG`, `MATCHERDEL/INS/REP`, `matchergap`; the installed `matcher` 6.6.0 as oracle.

`SequenceAligner.LocalAlignTopN` / `LocalAlignAffineTopN` (all overloads) take `TopNTieOrder tieOrder` (default `WatermanEggert` = unchanged behaviour). `TopNTieOrder.Sim` runs `SimMatcher` (`SequenceAligner.TopNSim.cs`), a line-by-line port of the program (nseq = 2): the K-vertex LIST keyed by the start point with its replacement rule and the static `low` / `most` pointers (`Addnode`, `Findmax`), the forward pass that records start points with the `MATCHERORDER` tie rule (`BigPass`), the recomputation of the affected area after each report (`Locate`, `NoCross`, `SmallPass`), the linear-space divide and conquer with `tb`/`te` end-gap states and the program's midpoint tie rule (`Diff`), and the used-pair rows (`MATCHERDIAG`; in `Diff` a forbidden diagonal leaves `c` at `CC[j] - qr`, ported as is). Integer scores: the scalar `ScoringMatrix` or the exact integer model of a `MatrixScoring` (so fractional gaps work, which the program cannot take). The list depends on K (the program's LIST size); when fewer than K alignments exist the program dies ("The number of alignments computed is too large"), the port returns the alignments found.

Oracle `matcher_sim_oracle.tsv` (`scripts/alignment/generate_matcher_sim_oracle.py`): every matcher row of `topn_oracle.tsv` (397, including the 7 that differ from Waterman-Eggert: lines 46, 487, 502, 787, 793, 865, 1348) plus 800 new tie-rich pairs (tandem repeats, repeated domains, alphabets of 2-4 letters, lengths up to 90, 7 gap schemes, K = 1-10; DNA EDNAFULL and protein EBLOSUM62): **1 197 rows, 1 138 complete lists equal alignment by alignment (score, 0-based coordinates, both rows), 59 rows where the program aborts and the port returns fewer than K alignments**. The tie rows are exact in Sim order and still differ in Waterman-Eggert order (asserted).

## Audit A2-28d (B09 F152): dotmatcher

Sources actually opened: EMBOSS 6.6.0 `emboss/dotmatcher.c` in full (the two sliding-window loops, `dotmatcher_pushpoint`, the `-stretch` / graph output) and the installed `dotmatcher` 6.6.0 (`-graph data` writes the plot as text: after the title lines the `Line x1 y1 x2 y2 colour 0` entries in drawing order, then `Rectangle` and the axes).

`DotMatcher.Compute` (Alignment, public) + MCP `dotmatcher`: for every diagonal the window sum of matrix scores is updated one pair per step; a maximal run of windows at or above the threshold is one segment (first residue of its first window to the last residue of its last window); loop 1 starts every diagonal on the first residue of the program's first sequence and at residue j of its second (j < length2 - window), loop 2 the other way round (i < length1 - window). **Measured: the program's "first" sequence is the SECOND input sequence** (the plot has the second sequence on the x axis; an earlier port with the C names failed 4 of 4 asymmetric probes, the swapped one reproduces all): a character that is not a letter of the matrix scores 0; sequences not longer than the window read past the string in the program and give no segments here. `MultipleAlignmentDisplay`'s private matrix-code class was replaced by `SubstitutionMatrix.TryGetScore` (the same rule: unknown characters score 0).

Oracle `dotmatcher_oracle.tsv` (`generate_dotmatcher_oracle.py`): **600 / 600 plots equal** (255 + 105 DNA with EDNAFULL / EDNASIMPLE, 240 protein with EBLOSUM62 / EPAM250 / EBLOSUM30; windows 3-20, thresholds a third to the whole of an identical window; 356 with segments, 2 000+ segments, 244 empty plots). Tests: `DotMatcher_Tests` (2), MCP `dotmatcher`. Still open in A2-28d: `wordmatch`, `est2genome`.
