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
