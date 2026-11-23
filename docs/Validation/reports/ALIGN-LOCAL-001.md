# Validation Report: ALIGN-LOCAL-001 — Local Alignment (Smith–Waterman)

- **Validated:** 2026-06-24   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` (canonical) and `LocalAlign(string, string, ScoringMatrix?)` (delegate), both routing to `LocalAlignCore` → `TracebackLocal` in `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

This is a re-validation in fresh context. The source file was last functionally touched for local alignment by the e19a8a02 ("O(n) traceback") perf work; the local core/traceback path reviewed below is logically identical to the previously-validated version. No defects found; no code changes made.

---

## B09 current state (2026-11-19) - supersedes the original verdict lines above and below

- **Stage A:** PASS-WITH-NOTES. **Stage B:** PASS-WITH-NOTES. **State:** FIXED (not CLEAN: the 2026-06 verdict below is superseded).
- **Fixes (B09.md, F-numbers):** F2 `LocalAlignAffine` (Smith-Waterman-Gotoh; `GapOpen` was ignored), F39 empty input = tagged-Local empty result on every overload, F41 `MatrixScoring`, F45 general gap costs, F46 exact linear-space local alignment (locate interval + global solve, SIM / EMBOSS matcher scheme), F47 Waterman-Eggert non-intersecting top-N alignments + capped co-optimal enumeration, F77 top-N in exact linear space above 25 M cells.
- **Oracle numbers:** F2 3 192 cases vs Biopython / parasail; F41 as GLOBAL; F46 forced recursion identical on 4 750 cases; F47 topn oracle 1 350 rows (independent full-recompute Waterman-Eggert, exhaustive enumeration of tiny pairs, EMBOSS matcher 6.6.0 on 397 rows / 1 695 alignments: score lists equal on 390, Biopython co-optimal sets on 1 345 rows); F77 480 random pairs + all 1 350 oracle rows forced through the linear path and two 5 kb pairs above the limit = the in-memory declumper.
- **Open items:** none doable. `LocalAlignCoOptimal` keeps the O(m·n) matrices (rejects > 25 M cells): an enumeration whose output can be exponential in the input has no linear-space form. The 2026-06 '7/7 tests' count is obsolete; see the TestSpec and `docs/Evidence/ALIGN-LOCAL-001-Evidence.md`.
- The unit's TestSpec (`tests/TestSpecs/ALIGN-LOCAL-001.md`), Evidence (`docs/Evidence/ALIGN-LOCAL-001-Evidence.md`) and algorithm document carry the same state; where the original text below says CLEAN / no defects / no code changes, that was the state of the first pass only.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "Smith–Waterman algorithm"** (https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm). Confirms: the local-alignment recurrence with the **zero floor**; first row/first column initialized to 0 (no end-gap penalty); the optimal score is the **maximum cell anywhere in the matrix** (not the corner); traceback starts at that max cell and **ends at the first cell with score 0**; linear (`W_k = k·W_1`) and affine (Gotoh) gap models; and the worked `TGTTACGG` / `GGTTGACTA` example (match +3, mismatch −3, gap 2 → score 13, alignment `GTT-AC` / `GTTGAC`). The page cites Smith, T.F. & Waterman, M.S. (1981), *J. Mol. Biol.* 147:195–197 as the primary reference.

### Formula check
Linear-gap S–W recurrence (Wikipedia §Linear gap penalty), matching TestSpec §1.2 point 5:

**H(i,j) = max( 0, H(i−1,j−1)+s(a_i,b_j), H(i−1,j) − W_1, H(i,j−1) − W_1 )**

- Zero floor (the `0` term) is the distinguishing feature vs Needleman–Wunsch and makes the alignment local.
- Initialization H(k,0)=H(0,l)=0.
- Optimum = max cell anywhere; traceback to first 0.

All match the spec (§1.2, §7) exactly.

### Edge-case semantics
- **No similarity:** every cell = max(0, neg, neg, neg) = 0 ⇒ max score 0, empty alignment (zero-floor property — sourced).
- **Identical sequences (match M>0, M>|W_1|):** diagonal dominates, H(i,i)=i·M, max = n·M, full alignment, no gaps (derivable).
- **Empty input:** a 0-length sequence leaves only the all-zero initialized row/column ⇒ max score 0 (sourced from scoring-matrix definition).

### Independent cross-check (hand computation)
Hand-recomputed the **full DP matrix** for the Wikipedia example A=`TGTTACGG`, B=`GGTTGACTA`, match +3, mismatch −3, linear gap W_1 = 2, recurrence `max(0, diag+s, up−2, left−2)` — reproduces the spec §7.1 matrix cell-for-cell:

```
       -   G   G   T   T   G   A   C   T   A
  -  [ 0   0   0   0   0   0   0   0   0   0 ]
  T  [ 0   0   0   3   3   1   0   0   3   1 ]
  G  [ 0   3   3   1   1   6   4   2   1   0 ]
  T  [ 0   1   1   6   4   4   3   1   5   3 ]
  T  [ 0   0   0   4   9   7   5   3   4   2 ]
  A  [ 0   0   0   2   7   6  10   8   6   7 ]
  C  [ 0   0   0   0   5   4   8  13  11   9 ]
  G  [ 0   3   3   1   3   8   6  11  10   8 ]
  G  [ 0   3   6   4   2   6   5   9   8   7 ]
```

**Maximum = 13 at H[6,7]** (row C = A[5], col C = B[6]; match → diag = H[5,6]+3 = 10+3 = 13). Traceback (diag → up → left priority):
H[6,7]=13 →diag H[5,6]=10 →diag H[4,5]=7 →**left** H[4,4]=9 →diag H[3,3]=6 →diag H[2,2]=3 →diag H[1,1]=0 STOP.
**Aligned region:** A[1..5]=`GTTAC` / B[1..6]=`GTTGAC` → `GTT-AC` / `GTTGAC`, score **13**; positions start1=1,end1=5, start2=1,end2=6. Matches Wikipedia's stated alignment exactly.

Swapped case (M2) hand-checked against spec §7.2: transposed matrix has max 13 at H^T[7,6], yielding `GTTGAC` / `GTT-AC` (gap moves to seq2 via the "up" branch).

### Findings / divergences
None. The description is biologically and mathematically correct and faithfully sourced.

---

## Stage B — Implementation

### Code path reviewed
`SequenceAligner.cs`: `LocalAlign(DnaSequence)` lines 287–296; `LocalAlign(string)` lines 301–313; `LocalAlignCore` lines 315–348; `TracebackLocal` lines 350–398. `ScoringMatrix` default `SimpleDna` at lines 22–26; gap penalty supplied via `GapExtend` (stored as a negative value, added to the predecessor).

### Formula realised correctly? (evidence)
- **Zero floor:** line 335 `score[i,j] = Math.Max(0, Math.Max(diag, Math.Max(up, left)))`. Present.
- **Recurrence terms:** diag = `score[i-1,j-1] + matchScore`, up = `score[i-1,j] + GapExtend`, left = `score[i,j-1] + GapExtend` (lines 331–333); `matchScore` = `Match` on equal chars else `Mismatch` (line 329). Linear gap = `GapExtend` per position; with the test's `GapExtend = -2` this equals Wikipedia's `W_1 = 2`. `GapOpen` is intentionally unused (no affine term). Correct.
- **Initialization:** `int[m+1, n+1]` defaults to 0 ⇒ row 0 / column 0 are 0 with no penalty loop. Correct.
- **Max-cell optimum:** lines 337–342 track `maxScore`/`maxI`/`maxJ` over every filled cell; traceback starts at `(maxI,maxJ)` (line 347). Does **not** take F(m,n). Correct.
- **Score reported:** `score[endI,endJ]` (line 392) = the max cell value. Correct.
- **Traceback-to-zero:** loop condition `while (i>0 && j>0 && score[i,j] > 0)` (line 359) terminates at the first 0 cell. Branch priority diag → up → left (lines 366–383) reproduces the spec's gap placement. Correct.
- **Positions:** 0-indexed origin/terminus of the traced region (lines 394–397), within bounds (INV-5).

### Cross-verification table recomputed vs code (tests run)
| Case | Input | Expected (source/hand) | Code result | Match |
|------|-------|------------------------|-------------|-------|
| M1 Wikipedia | TGTTACGG / GGTTGACTA, +3/−3/−2 | score 13, `GTT-AC`/`GTTGAC`, (1,5)/(1,6) | identical | ✅ |
| M2 swapped | GGTTGACTA / TGTTACGG | score 13, `GTTGAC`/`GTT-AC`, (1,6)/(1,5) | identical | ✅ |
| M3 string overload | same as M1 | score 13, `GTT-AC`/`GTTGAC` | identical | ✅ |
| M4 empty | "" / "ACGT" | `AlignmentResult.Empty` | identical | ✅ |
| M5 null | null DnaSequence | `ArgumentNullException` | thrown | ✅ |
| S1 identical | ACGTACGT ×2, +3 | score 24, full match | identical | ✅ |
| S2 dissimilar | AAAA / TTTT | score 0, empty alignment | identical | ✅ |

M1/M2 hand-recomputed in Stage A; all 7 confirmed by running `SequenceAligner_LocalAlign_Tests` (7/7 passed).

### Variant/delegate consistency
String overload applies `ToUpperInvariant()` then calls the same `LocalAlignCore`; returns `AlignmentResult.Empty` on empty/null-or-empty input (lines 306–312). DnaSequence overload null-guards then calls the identical core. M3 confirms parity.

### Test quality audit
Tests assert exact sourced values (score 13, exact aligned strings, exact 0-indexed positions), the zero-floor and AlignmentType invariants, gap-removal-to-substring (INV-4), and both gap-direction traceback branches (left in M1, up in M2). Deterministic; no tautology-only assertions. Edge cases (empty, null, all-match, all-mismatch) covered.

### Findings / defects
None.

---

## Verdict & follow-ups
- **Stage A: PASS** — recurrence (with zero floor), initialization, max-cell optimum, linear gap model, and traceback-to-zero all match Wikipedia / Smith–Waterman (1981); the full Wikipedia DP matrix was hand-recomputed to max score 13 with alignment `GTT-AC`/`GTTGAC`.
- **Stage B: PASS** — the implementation faithfully realises the validated recurrence: zero floor (line 335), optimum from the max cell not F(m,n) (lines 337–347), traceback terminates at 0 (line 359), linear gap via `GapExtend`. Hand-computed example reproduced by the code.
- **End-state: ✅ CLEAN.** No defects found; no code changes required. Build green; `SequenceAligner_LocalAlign_Tests` 7/7 pass.

## 2026-11 audit round 1 correction (A1-1, B09 F39)

The "M4 empty -> `AlignmentResult.Empty`" row and the string-overload remark above are **superseded**: an empty (or both-empty) side now gives the tagged-Local empty result `("", "", 0, Local, -1, -1, -1, -1)` on every overload (string/`DnaSequence`, linear/affine), and a null string throws `ArgumentNullException`. The Global-tagged `AlignmentResult.Empty` sentinel is no longer returned by any aligner.

## 2026-11 audit round 1 (A1-2, B09 F41)

Matrix-scored `LocalAlign`/`LocalAlignAffine` (water-compatible gap conventions): see ALIGN-GLOBAL-001 report and B09.md F41.

## 2026-11 audit round 1 (A1-5, B09 F45)

Arbitrary gap-cost functions are implemented (`LocalAlignGeneralGap`, Waterman-Smith-Beyer single-run model, Biopython callable `gap_score` oracle); the "not implemented: arbitrary W_k" statement is superseded (`Local_Alignment_Smith_Waterman.md` §5.3).

## 2026-11 audit round 1 (A1-6, B09 F46)

Exact linear-space alignment (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991 scheme): every pairwise aligner switches automatically above 25 000 000 cells (identical score, rows, coordinates, ties), explicit `*LinearSpace` methods for any size, MCP automatic. Oracles: Biopython 1.88 / parasail 1.3.4 (1 280 pairs), EMBOSS `stretcher` / `matcher` 6.6.0 (160 + 160 rows; large-pair scores equal up to 50 kb), full-matrix differential (4 750 forced-recursion cases + kb pairs + the automatic switch). 20 kb x 20 kb affine: 111 MB peak (4.8 GB full); 100 kb x 100 kb: 185 MB. Memory guards for the non-affine general-gap aligners (4 M cells) and the profile MSA DP (50 M cells). Algorithm document: `Linear_Space_Alignment.md`.


## 2026-11 audit round 1 (A1-7, B09 F47)

The "one optimal alignment only" limitation is lifted: Waterman-Eggert non-intersecting top-N (`LocalAlignTopN` / `LocalAlignAffineTopN`, the definition of EMBOSS `matcher -alternatives` / `lalign`, verified against matcher's source and 1 695 of its alignments) and the capped co-optimal enumeration (`LocalAlignCoOptimal`, Biopython `alignments` set) are implemented on the existing local cores (`Local_Alignment_Smith_Waterman.md` §5.5; MCP `local_align_top_n`, `local_align_co_optimal`).

## 2026-11 audit round 1 (A1-7c, B09 F77)

The O(m·n) memory limit of the top-N declumper is lifted: above 25 000 000 cells `LocalAlignTopN` / `LocalAlignAffineTopN` run the exact linear-space program (the aligned pairs of the earlier alignments are forbidden in the locate passes and interval solves of F46; SIM / EMBOSS `matcher` scheme), with identical alignments (480 random pairs, 1 350 oracle rows, two 5 kb pairs; 40-46 MB peak at 5-20 kb). `LocalAlignCoOptimal` keeps the matrices. Evidence file `docs/Evidence/ALIGN-LOCAL-001-Evidence.md` created (A1-8).

## Audit r2 addendum (A2-28, B09 F122): water

`SequenceAligner.WaterAlign` ports EMBOSS water (rows, starts and scores equal on 860 / 860 pairs).

## Audit r2 addendum (A2-28c, B09 F124): supermatcher

`SequenceAligner.Supermatcher` (Rabin-Karp word seeds + banded Smith-Waterman, the program's single-precision arithmetic); MCP `supermatcher_align`; 260 / 260 runs equal the installed supermatcher 6.6.0 (output files and error files).

## Audit A3-14b / A3-14c / A3-17 addendum (B09 F141 / F143): water border rule, matcher ties

MCP `water_align` (A3-17). The "at least 90 % reaches the textbook optimum" tolerance of water is replaced by an independent score-only program of `embAlignPathCalcSW` (639 / 639 dyadic oracle rows; the maximum is taken over the cells with both indices at least 1 only, so an alignment ending in the first row or column is never reported: 7 of the 250 sampled rows, ids 23, 51, 69, 240, 292, 305, 405). F47's "= EMBOSS matcher -alternatives" holds for round 1 (397 / 397) and for the whole score list on 390 / 397 rows; the 7 others (lines 46, 487, 502, 787, 793, 865, 1348 of `topn_oracle.tsv`) are explained by a tie broken differently in an earlier round (`scripts/alignment/analyze_topn_matcher_ties.py`); the exact tie order of SIM needs the matcher port (audit A2-28d).

## Audit A2-28d (B09 F145): matcher port

`TopNTieOrder.Sim` = EMBOSS `matcher` 6.6.0 (Huang-Miller SIM) ported from `matcher.c`; 1 138 / 1 138 alignment lists equal (incl. the 7 tie rows), 59 aborting rows return fewer than K; default order unchanged; MCP `local_align_top_n` `tieOrder`.

- **Audit A2-28d (F152):** EMBOSS `dotmatcher` ported (`DotMatcher.Compute`, MCP `dotmatcher`): 600 / 600 plots equal.

- **Audit A4-3 (F153):** the EMBOSS reader drops gap characters and every other non-letter (except `*`) from plain input sequences; one shared `EmbossSequenceReader` (Infrastructure) serves water, matcher, supermatcher, dotmatcher: 240 oracle rows each equal (`emboss_gap_oracle.tsv`).

- **Audit A4-1 (F157):** `DotMatcher.Compute`, matcher (`TopNTieOrder.Sim`) poll the ambient cancellation token (mid-run cancel tests); MCP tools take the request token.

- **Audit A4-2 (F158):** resource guards: `DotMatcher.Compute` (`DotMatcherMaxSegments`), `TopNTieOrder.Sim` (`TopNSimMaxVertices`, lazy vertex list: `maxAlignments = int.MaxValue` = every alignment).

- **Audit r5 / r6 (A5-7, A6-2, A6-3, A6-5; F172, F177, F178, F180):** `water` / `matcher` / `supermatcher` read the sequence-input qualifiers (180 oracle rows each); `MatcherProgram` accepts every integer `-gapopen` / `-gapextend` pair including `-gapopen` below `-gapextend` (SIM `q = gapopen - gapextend` < 0), resets a negative value with the program's warning (`EmbossProgramResult.StandardError`) and refuses a non-integer (872 rows); `Supermatcher` with `-wordlen` above every target prints `Error: no kmers found` and returns without the Rabin-Karp loops (`-wordlen 2147483647` took 30 s; 260 rows).
- **Audit A2-28d (F203, WP48):** EMBOSS `est2genome` ported (`EstEngine`, `EmbossEst2genomePrograms`, MCP `est2genome`): the -outfile text of 330 installed-program runs equal (all options, the linear-space recursion, ESTs on both strands, the crash rows); see docs/algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md.
- **Audit A5-8c3e (F241 - F243, WP54):** the feature library of EMBOSS 6.6.0 ported - vocabularies, the readers of eleven input formats with the library's file buffer and messages, the feature file of a sequence (`-ufo`, `-fformat`, `-fopenfile`: `diffseq`), the writers gff3 / gff2 / embl / genbank / debug, `featcopy` (MCP `featcopy`); 2 600 + 1 200 + 450 + 1 779 oracle runs equal.
- **Audit A5-8c3f (F244, F245, WP55):** the feature writers bed / refseq / refseqp / pir / swiss / dasgff / draw for a full table (`featcopy`: 2 750 / 2 750 oracle runs equal) and `featmerge` (MCP `featmerge`: 450 / 450).
- **Audit A5-8c3f (F246, WP55):** `feattext` (MCP `feattext`): 400 / 400 oracle runs equal.
- **Audit A5-8c3f (F247, WP55):** the list forms `@list` / `list::list` of the `features` input (`EmbossFeatureInput`; MCP `files` argument): 500 / 500 oracle runs equal.
- **Audit A5-8c3g (F248, WP56):** the sequence readers of `ajseqread.c` for 49 input formats (fasta, gcg, embl, swiss with its DE description, nbrf, genbank, refseqp, genpept, codata, strider, ig, staden, text, fastq, sam, gff2, gff3, gde, raw, clustal, msf, fitch, mase, biomart, treecon, acedb and their aliases; the feature tables of embl / genbank / swiss / nbrf / gff / refseqp inputs) and `seqret` (MCP `seqret`): 9 800 / 9 800 oracle runs equal.
- **Audit A5-8c3g (F249, WP56):** `maskfeat` (MCP `maskfeat`): 1 200 / 1 200 oracle runs equal.
- **Audit A5-8c3g (F250, WP56):** `descseq` (MCP `descseq`): 1 000 / 1 000 oracle runs equal.
- **Audit A5-8c3h (F251, WP57):** the sequence readers of phylip, phylipnon, nexus / paup, mega, selex, stockholm / pfam, jackknifer, hennig86, experiment, ace / consed, pdb, pdbseq, pdbnuc and pdbnucseq: 3 400 / 3 400 oracle runs equal (the runs on which the reference reads freed memory, overruns an array, loops or dies are refused).
- **Audit A5-8c3h (F252, WP57):** the sequence writers of gcg, nbrf, codata, strider, staden, text, ig, iguspto, gde, fitch, acedb, mase and the `-osformat` / `-ossingle` / `-osdirectory` / `-osextension` options of `seqret`: 3 800 / 3 800 oracle runs equal.
- **Audit A5-8c3h (F253, F254, WP57):** `union` and `seqretsplit` (MCP `union`, `seqretsplit`): 200 / 200 oracle runs equal each.
- **Audit A5-8c3h (F255, WP58):** the USA forms of the sequence input (ranges, list files, queries on files, `-sbegin` / `-send` / `-sreverse` on the output and its features, `-sid`, `-sdbname`, `-iquery`, `-ioffset`, `-squick`) of seqret, seqretsplit, union, maskfeat and descseq: 1 500 / 1 500 oracle runs equal (27 refused with proof); MCP `sequenceValue`.
- **Audit A5-8c3h (F256, WP58):** the sequence databases of an EMBOSS configuration file (`DBNAME`, method `direct`, `db:id` queries): 600 / 600 oracle runs equal (14 refused with proof).
- **Audit A5-8c3h (F257 - F263, WP58):** the feature programs `entret` (756 oracle runs), `extractfeat`, `coderet`, `showfeat`, `union -feature -source`, `featreport` and `twofeat` (1 000 each): equal to the installed 6.6.0 programs byte for byte (the runs on which the reference dies by a fault or never stops are refused with proof); MCP tools `entret`, `extractfeat`, `coderet`, `showfeat`, `featreport`, `twofeat` and the feature parameters of `union` (603 tools).
- **Audit A5-8c3h (F264, WP58):** the sequence writers ncbi, gifasta, nib / nibble, dasdna, das and asn1 (MCP `seqret` `outputFormat`, `union`, `extractfeat`, `coderet`): 260 / 260 oracle runs equal per format (1 820).
- **Audit A5-8c3h (F265 - F274, WP59):** every sequence writer of EMBOSS 6.6.0 (fastq family, sam, bam, scf, gff2 / gff3, embl / genbank / swiss families, experiment, debug and the set formats msf ... treecon) with the annotation the readers keep, the features written inside the sequence file, `-osdbname`, `-osdirectory` / `-osextension` / `-osname` / `-ossingle` of seqret, extractfeat and coderet: about 18 000 `seqret` runs and 2 200 `extractfeat` / `coderet` runs equal to the installed programs byte for byte (the runs on which the reference reads memory it never wrote or dies by a fault are refused with proof); the BGZF blocks of BAM by a port of zlib 1.3 deflate; MCP `seqret`, `union`, `extractfeat`, `coderet` (603 tools).
- **Audit A5-8c3h (F275 - F279, WP60):** the databases of the access methods `emboss` (index files of `dbxfasta` / `dbxflat` / `dbxcompress`) and `emblcd` (`dbifasta` / `dbiflat`), the `RES` / `SERVER` / `INCLUDE` / `IF` blocks of the configuration, `coderet -sreverse1` on GFF3, the output qualifiers of `union`, `-rname` / `-rextension` / `-rdirectory` and the features input qualifiers of `featreport` / `twofeat`: 1 200 + 1 200 + 400 database runs, 1 200 `union`, 1 000 + 1 000 report runs and 1 100 `coderet` runs equal to the installed programs byte for byte; registered: the methods gcg / embossgcg / blast / app, the DRCAT lookup (network: Docker package), the byte readers abi / bam / scf, das, ensembl, embl `CO` lines, the indexer programs.


WP61 (F280): 23 EMBOSS sequence programs (getorf, transeq, revseq, extractseq, cutseq, degapseq, trimseq, skipseq, nthseq, notseq, pasteseq, trimest, newseq, yank, seqcount, splitter, maskseq, maskambignuc, maskambigprot, trimspace, nospace, notab, noreturn) run from the command line (`EmbossProgramRunner`), each equal to the installed program on 320 runs (getorf 800); new `EmbossTrn`, `EmbossSeqType`, `EmbossSeqoutStream`, `EmbossAcdArgs`; see B09.md F280.

WP62 (F281): 11 more EMBOSS programs (shuffleseq, msbar, makenucseq, makeprotseq, nthseqset, listor, infoseq, wordcount, seqxref, showorf, prettyseq) run from the command line (`EmbossProgramRunner`), each equal to the installed program on 300 runs under a fixed clock; the clock / `timetoday` infrastructure (`EmbossRunEnvironment`), the generator of `ajmath.c` (`EmbossRandom`), the set input (`EmbossSeqset`), `ajSeqoutWriteSet`; see B09.md F281.

WP63 (F282): 5 more EMBOSS programs (tcode, checktrans, plotorf, wordfinder, marscan) run from the command line (`EmbossProgramRunner`), each equal to the installed program on 300 runs; the ACD adapters of the report layer (`EmbossReportArgs`), the XY graph (`EmbossGraphXy`: `-graph data` / `bedgraph` / `none`), the feature output (`EmbossFeatout`) and the alignment output (`EmbossAlignArgs`); library fixes found by the oracle: the default name of an `outfile` (the name of the first sequence, not the file), the blank under identical non-letters (`*`) in the alignment markup, `EMBOSS_001` names of unnamed rows in the sequence-style alignment formats, unknown `format::` prefixes of sequence outputs; 5 MCP tools (Analysis server).

- **Audit R14-1 (B09 F288):** the End coordinates of local / semi-global results follow the single 0-based inclusive convention of every producer (see ALIGN-GLOBAL-001); the EMBOSS formats of `local_align` / `semi_global_align` take the offsets from the rows (a local result is a water report; the fitting rows carry the whole reference, so it runs 1 .. length).

WP67 (F306): `patmatmotifs`, `patmatdb` and `prosextract` of EMBOSS 6.6.0 run from the command line (`EmbossProgramRunner`), each equal to the installed program on 300 runs and, for `patmatmotifs`, to the 85 runs of the user-run Docker package on the real PROSITE release 2026_04 (no run reads invalid memory); the regular expression engine of the library (`EmbossPcre`, PCRE 7 incl. its compile errors) equals the library on 41 670 rows; defects found by the oracles: the order of the matches (`ajListPush` front), the source and strand of the features, the message of a bad boolean (`Invalid Y/N value`), a NUL in a line of a sequence file, the stale bytes of `ajStrRemoveWhiteExcess`, the directory `PROSITE/` opened as a documentation file; 3 MCP tools (catalog 653 tools).

WP67 (F309): the USAs of databases that only the catalogue DRCAT defines (899 runs of `seqret` / `entret`) and `drfinddata`, `dbtell`, `showdb` of EMBOSS 6.6.0 (233 / 302 / 301 runs, standard output, error stream, exit status and files byte for byte; the user's run with network access integrated) equal the installed programs; not ported: user configuration files, `dbtell -server`, `drfinddata -odirectory` (registered A5-8c3h-db).
