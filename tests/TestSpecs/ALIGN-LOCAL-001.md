# Test Specification: ALIGN-LOCAL-001

**Test Unit ID:** ALIGN-LOCAL-001
**Area:** Alignment
**Algorithm:** Local Alignment (Smith–Waterman)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-28

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Smith–Waterman algorithm | https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | 2026-03-06 |

### 1.2 Key Evidence Points

1. Smith–Waterman performs local sequence alignment, finding the highest-scoring pair of subsequences between two sequences. (Smith–Waterman algorithm)
2. The key difference from Needleman–Wunsch is the **zero floor**: negative scores are set to 0, allowing alignment to restart at any position. (Smith–Waterman algorithm, §Comparison)
3. Initialization: first row and first column are set to 0 (no end-gap penalty). (Smith–Waterman algorithm, §Scoring matrix)
4. Traceback: begins at the cell with the **highest score** in the matrix and ends when a cell with score 0 is encountered. (Smith–Waterman algorithm, §Algorithm step 4)
5. Linear gap penalty recurrence: $H_{ij} = \max(0,\; H_{i-1,j-1} + s(a_i,b_j),\; H_{i-1,j} - W_1,\; H_{i,j-1} - W_1)$ (Smith–Waterman algorithm, §Linear gap penalty)
6. Wikipedia example: sequences `TGTTACGG` and `GGTTGACTA` with match +3, mismatch −3, linear gap penalty $W_1 = 2$ yields alignment `GTT-AC` / `GTTGAC` with score 13. (Smith–Waterman algorithm, §Example)

### 1.3 Documented Corner Cases

| Case | Evidence Source | Notes |
|------|-----------------|-------|
| No similarity between sequences | Smith–Waterman zero floor (§Algorithm) | When all comparisons yield negative scores, every $H_{ij} = \max(0, \text{neg}, \text{neg}, \text{neg}) = 0$. Score = 0, alignment is empty. |
| Identical sequences | Derivable from S-W recurrence | For identical sequences of length $n$ with match score $M > 0$ and $M > \|W_1\|$: diagonal dominates at every cell, $H_{i,i} = i \times M$. Maximum $= n \times M$, full alignment, no gaps. |
| Empty sequence input | Smith–Waterman definition (§Scoring matrix) | Matrix dimensions are $1 \times (m+1)$ or $(n+1) \times 1$. Only initialized row/column (all zeros) — no cells to fill. Max score = 0. |

### 1.4 Known Failure Modes / Pitfalls (from sources)

1. Multiple optimal local alignments can exist when several cells share the same maximum score. Implementation returns one. (Smith–Waterman algorithm, §Example)
2. Alignment quality depends on the scoring system and gap penalties. (Smith–Waterman algorithm, §Explanation)
3. Linear gap penalty may not model biological indels accurately; affine gaps (Gotoh) are often preferred. (Smith–Waterman algorithm, §Affine gap penalty)

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` | SequenceAligner | **Canonical** | Smith–Waterman local alignment |
| `LocalAlign(string, string, ScoringMatrix?)` | SequenceAligner | Delegate | String wrapper — delegates to same `LocalAlignCore` after `ToUpperInvariant()` |
| `LocalAlignAffine(DnaSequence, DnaSequence, ScoringMatrix?)` / `(string, string, ScoringMatrix?)` | SequenceAligner | **Canonical (affine)** | Smith-Waterman-Gotoh; gap of length k = GapOpen + k·GapExtend; shares `AffineAlignCore` with `GlobalAlignAffine` (2026-09 review) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Score ≥ 0 (zero floor property) | Yes | Smith–Waterman algorithm: $\max(0, \ldots)$ in recurrence |
| INV-3 | AlignmentType is Local | Yes | Implementation contract |
| INV-4 | Removing gaps from aligned sequences yields substrings of originals | Yes | Smith–Waterman traceback rules |
| INV-5 | StartPosition and EndPosition are within sequence bounds | Yes | Smith–Waterman traceback from matrix positions |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Wikipedia example: exact score and alignment (gap in seq1, "left" traceback branch) | seq1=`TGTTACGG`, seq2=`GGTTGACTA`, scoring: match +3, mismatch −3, gap −2 | Score=13, aligned1=`GTT-AC`, aligned2=`GTTGAC`, start1=1, end1=5, start2=1, end2=6, INV-1,3,4,5 hold | Smith–Waterman §Example, hand-verified scoring matrix |
| M2 | Swapped Wikipedia example: gap in seq2, "up" traceback branch | seq1=`GGTTGACTA`, seq2=`TGTTACGG`, same scoring | Score=13, aligned1=`GTTGAC`, aligned2=`GTT-AC`, start1=1, end1=6, start2=1, end2=5, INV-3,4,5 hold | Derived from Wikipedia via S-W transposition: $H^{AB}_{i,j} = H^{BA}_{j,i}$; hand-verified transposed matrix |
| M3 | String overload produces same Wikipedia-expected result | Same Wikipedia sequences as string | Score=13, aligned1=`GTT-AC`, aligned2=`GTTGAC` | Implementation contract — delegates to `LocalAlignCore`; validated against Wikipedia expected values |
| M4 | Empty input (either side or both, all overloads) returns the tagged-Local empty result `("", "", 0, Local, −1, −1, −1, −1)` — not the Global-tagged `AlignmentResult.Empty`; null → `ArgumentNullException` | seq1="", seq2="ACGT" (and mirrored, both empty) | Local empty result | S-W definition: 0-dimensional scoring matrix → maxScore=0 |
| M5 | Null DnaSequence throws `ArgumentNullException` | seq1=null, seq2=valid | ArgumentNullException | .NET convention: `ArgumentNullException.ThrowIfNull` |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | Identical sequences produce full-length alignment with exact score | seq1=seq2=`ACGTACGT`, scoring: match +3, mismatch −3, gap −2 | Score=24 (8×3), full alignment `ACGTACGT`/`ACGTACGT` | Derivable: identical seqs → diagonal dominates → score = n×match |
| S2 | Completely dissimilar sequences produce score exactly 0 | seq1=`AAAA`, seq2=`TTTT`, scoring: match +3, mismatch −3, gap −2 | Score=0, empty alignment | Derivable: no matches → all $H_{ij} = \max(0, \text{neg}, \text{neg}, \text{neg}) = 0$ |

### 4.2b Tie-breaking and affine-gap tests — added 2026-09 review

Reference values: Biopython 1.88 `PairwiseAligner(mode='local', open_gap_score=GapOpen+GapExtend, extend_gap_score=GapExtend)` and parasail 1.3.4 `sw_trace(open=-(GapOpen+GapExtend), extend=-GapExtend)`; scores identical on 3192 random/edge cases (linear + affine), every Seqeron alignment re-scores to its score, is among Biopython's co-optimal alignments, and its end cell is the row-major-first maximum.

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| T1 | Tied optima → first row-major end | ACGTTTTTACGT / ACGT, SimpleDna | 4, ACGT, (0,3)/(0,3) | Biopython first alignment; parasail end 3/3 |
| T2 | Traceback stops at first zero cell | ACGTTTTACGT / ACGTACGT, BlastDna (linear) | 10, TACGT, (6,10)/(3,7) | parasail; one of Biopython's 4 optima |
| A1 | Long gap charged one opening | ACGTACGTAAAAAAACGTACGT / ACGTACGTACGTACGT, HighIdentityDna | affine 64 (end 21/15); linear 74 | Biopython = parasail |
| A2 | Opening cost prevents bridging | ACGTACGTTTTTTACGT / ACGTACGTACGT, (1,−1,−5,−1) | 8, ACGTACGT, (0,7)/(0,7) | Biopython = parasail |
| A3 | Reference values | GGGACGTTTACGTCCC / TTACGTACGTAA BlastDna → 12 (7,12)/(0,5); ACGTTTTACGT / ACGTACGT BlastDna → 10 (6,10)/(3,7) | as listed | Biopython = parasail |
| A4 | GapOpen = 0 reduces to linear SW | Wikipedia example | 13, GTT-AC/GTTGAC | W_k = k·e |
| A5 | No positive region / empty input | AAAA/TTTT, ""/ACGT, ""/"" | Score 0, "", Local, coordinates −1 | same as typed linear path |
| A6 | Ties, overloads, nulls | AC/CA SimpleDna; DnaSequence vs lower-case string; null | A/A (0,0)/(1,1) = linear; equal; ArgumentNullException | API contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| C1 | Multiple optimal alignment detection | Sequences with tied max scores | Any optimal alignment is acceptable | Evidence: multiple optimal alignments can exist (§Example) |

---

## 5. Coverage Classification

| Area | Status | Tested by | Notes |
|------|--------|-----------|-------|
| Wikipedia example exact values (score, alignment, positions) | ✅ | M1 | Score=13, `GTT-AC`/`GTTGAC`, positions (1,5)/(1,6) |
| Traceback "left" branch (gap in seq1) | ✅ | M1 | `GTT-AC` has gap at position 4 in aligned seq1 |
| Traceback "up" branch (gap in seq2) | ✅ | M2 | `GTT-AC` has gap at position 4 in aligned seq2 |
| Traceback "diagonal" branch (match/mismatch) | ✅ | M1, M2 | Both tests traverse multiple diagonal steps |
| S-W transposition symmetry | ✅ | M1+M2 | Same score=13 for swapped sequences |
| Zero floor (INV-1) | ✅ | M1, S2 | Positive (13) and zero (0) boundary values |
| AlignmentType is Local (INV-3) | ✅ | M1, M2 | Checked in both dataset tests |
| Gaps removal → substring (INV-4) | ✅ | M1, M2 | Exact gap-removed strings + substring containment |
| Position bounds (INV-5) | ✅ | M1, M2 | Exact positions from hand-verified traceback |
| String overload parity | ✅ | M3 | Validated against Wikipedia expected values |
| Empty input handling | ✅ | M4, A5 | Returns the tagged-Local empty result (−1 coordinates) |
| Null argument handling | ✅ | M5 | `ArgumentNullException` |
| Identical sequences (all-match boundary) | ✅ | S1 | Score=24 (8×3), full alignment |
| Dissimilar sequences (all-mismatch boundary) | ✅ | S2 | Score=0, empty alignment |

---

## 6. Open Questions / Decisions

None.

---

## 7. Hand-Verified Scoring Matrices

### 7.1 Wikipedia Example (M1)

Sequences: A = `TGTTACGG`, B = `GGTTGACTA`
Scoring: match = +3, mismatch = −3, linear gap $W_1 = 2$

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

**Maximum**: 13 at H[6,7] (row=C, col=C)

**Traceback**: H[6,7]=13 →(diag) H[5,6]=10 →(diag) H[4,5]=7 →(**left**) H[4,4]=9 →(diag) H[3,3]=6 →(diag) H[2,2]=3 →(diag) H[1,1]=0 STOP

**Result**: `GTT-AC` / `GTTGAC`, score = 13

### 7.2 Swapped Wikipedia Example (M2)

Sequences: A = `GGTTGACTA`, B = `TGTTACGG`
Derived from §7.1 via S-W transposition: $H^T[i,j] = H[j,i]$

```
       -   T   G   T   T   A   C   G   G
  -  [ 0   0   0   0   0   0   0   0   0 ]
  G  [ 0   0   3   1   0   0   0   3   3 ]
  G  [ 0   0   3   1   0   0   0   3   6 ]
  T  [ 0   3   1   6   4   2   0   1   4 ]
  T  [ 0   3   1   4   9   7   5   3   2 ]
  G  [ 0   1   6   4   7   6   4   8   6 ]
  A  [ 0   0   4   3   5  10   8   6   5 ]
  C  [ 0   0   2   1   3   8  13  11   9 ]
  T  [ 0   3   1   5   4   6  11  10   8 ]
  A  [ 0   1   0   3   2   7   9   8   7 ]
```

**Maximum**: 13 at H^T[7,6] (row=C, col=C)

**Traceback**: H^T[7,6]=13 →(diag) H^T[6,5]=10 →(diag) H^T[5,4]=7 →(**up**) H^T[4,4]=9 →(diag) H^T[3,3]=6 →(diag) H^T[2,2]=3 →(diag) H^T[1,1]=0 STOP

**Result**: `GTTGAC` / `GTT-AC`, score = 13

## 2026-11 audit round 1 addition (A1-2, B09 F41)

Matrix-scored `LocalAlign`/`LocalAlignAffine` (BLOSUM62/PAM250/EDNAFULL, fractional gaps): see the table in ALIGN-GLOBAL-001.md (`SequenceAligner_MatrixScoring_Tests`; local rows = EMBOSS water 6.6.0 = Biopython `mode='local'` = parasail `sw`).

## 2026-11 audit round 1 addition (A1-5, B09 F45): general gap-cost functions

`LocalAlignGeneralGap` (Smith-Waterman with Waterman-Smith-Beyer gaps; Biopython 1.88 `PairwiseAligner(mode='local', gap_score=callable)`): 480 oracle rows + exhaustive substring-pair enumeration + affine/linear equality; see ALIGN-GLOBAL-001.md (A1-5 section, GG-1..GG-7).

## 2026-11 audit round 1 addition (A1-6, B09 F46): exact linear-space alignment

`LocalAlign` / `LocalAlignAffine` (resp. `SemiGlobalAlign` / `SemiGlobalAlignAffine`, all 16 `EndGaps` policies) switch automatically above `SequenceAligner.LinearSpaceCellThreshold` (25 M cells) to the exact linear-space program (locate the optimal interval / start vertex, then Hirschberg / Myers-Miller on the rectangle); explicit `LocalAlignLinearSpace` / `LocalAlignAffineLinearSpace`. Result identical to the full-matrix aligner (score, rows, coordinates, ties). Tests LS-1 … LS-12 are listed in `ALIGN-GLOBAL-001.md` (`SequenceAligner_LinearSpace_Tests`, oracle Biopython 1.88 / parasail 1.3.4 / EMBOSS `matcher` + `stretcher` 6.6.0).


## 2026-11 audit round 1 addition (A1-7, B09 F47): Waterman-Eggert top-N and co-optimal local alignments

`LocalAlignTopN` / `LocalAlignAffineTopN` (non-intersecting local alignments: no aligned residue pair shared; alignment 1 = `LocalAlign`; scores non-increasing; `maxAlignments`, `minScore`) and `LocalAlignCoOptimal` / `LocalAlignAffineCoOptimal` (all optimal local alignments, capped; Biopython `alignments` set), `AlignedPairs`. Tests (`SequenceAligner_TopN_Tests` x16, `Seqeron.Mcp.Alignment.Tests.TopNAlignTests` x4): oracle `TestData/Alignment/topn_oracle.tsv` (1 350 rows, `scripts/alignment/generate_topn_oracle.py`: independent full-recompute Waterman-Eggert, exhaustive enumeration, EMBOSS matcher 6.6.0 verified round by round, Biopython 1.88 co-optimal set) - see ALIGN-GLOBAL-001-Evidence.md; exhaustive-optimum check of every round on 600 tiny random pairs; first alignment = `LocalAlign`; prefix stability; hand-computed tandem repeat (3 x ACGTAC, 12/12/12 at 0/9/18), gap-placement set (GAAC/GAAAC: G-AAC, GA-AC, GAA-C = 19); `minScore`/`maxAlignments`; DNA/protein/matrix overloads; null/range/memory-guard validation.

## 2026-11 audit round 1 addition (A1-7c, B09 F77): Waterman-Eggert top-N in linear space

`LocalAlignTopN` / `LocalAlignAffineTopN` (all overloads) switch above `SequenceAligner.TopNMaxCells` (25 M cells) from the in-memory declumper to the exact linear-space program (`ForbiddenPairs` threaded through the locate passes and interval solves; explicit `LocalAlignTopNLinearSpace` / `LocalAlignAffineTopNLinearSpace`); `LocalAlignCoOptimal` still keeps the matrices (rejects > 25 M cells). Tests (`SequenceAligner_TopNLinearSpace_Tests` x10; `SequenceAligner_TopN_Tests` memory-guard test replaced by `TopN_AboveTheCellLimit_UsesTheLinearSpaceProgram_AndCoOptimalStillRejectsIt`; MCP `TopNAlignTests` +1): 480 random tie-rich pairs (scalar + matrix, linear + affine, base-case sizes 1 / 9 / 40 / 200) forced through the linear path = the in-memory declumper; all 1 350 rows of `topn_oracle.tsv` (+ the matcher score lists on >= 90 % of 397 rows); the automatic switch with a lowered threshold; two 5 kb pairs above the limit = the in-memory declumper (forced with `fullMatrixCells = long.MaxValue`), alignment 1 = `LocalAlign` / `LocalAlignAffine`, no aligned pair twice; a positive gap score above the limit is rejected with a clear error; empty / null / range validation. Sources and oracle numbers: `docs/Evidence/ALIGN-LOCAL-001-Evidence.md` (created in audit item A1-8).

## Audit A2-28c (B09 F124): supermatcher
`SequenceAligner_EmbossPrograms_Tests` (`emboss_multi_oracle.tsv`): 260 supermatcher runs against the installed EMBOSS 6.6.0 (1-3 queries x 1-4 targets, 4-600 residues, DNA / IUPAC / protein, every option: matrix, gapOpen, gapExtend, wordLength 3-15, width 1-100, minScore; whole output files in 12 -aformats and the error files, byte for byte), plus the no-shared-word error line.

## Audit A2-28d (B09 F145): `TopNTieOrder.Sim` (EMBOSS matcher's tie order)

`SequenceAligner_TopNSim_Tests` (`matcher_sim_oracle.tsv`, 1 197 rows of the installed `matcher` 6.6.0): (1) every list equals the program's alignment by alignment (1 138 rows) and the 59 aborting rows return fewer than K alignments; (2) the 7 tie rows of `topn_oracle.tsv` are exact in Sim order and differ in Waterman-Eggert order; (3) every Sim list is a legal Waterman-Eggert sequence (non-increasing scores, no aligned pair reused); (4) the first alignment has the optimal score, empty input gives nothing, `maxAlignments < 1` and an unknown order throw; (5) `minScore` cuts the list. MCP `TopNAlignTests.TopN_TieOrderSim_ReproducesEmbossMatcher_WhereWatermanEggertOrderDiffers`.

## Audit A2-28d (B09 F152): dotmatcher

`DotMatcher_Tests` (`dotmatcher_oracle.tsv`, 600 dotmatcher 6.6.0 plots): `Compute_EqualsDotmatcher_AllRows` (segments in drawing order, DNA / protein, six matrices); `DefaultMatrix_IsChosenByTheSequenceType_AndShortSequencesGiveNoSegments`. MCP `Dotmatcher_...`.

## Audit A4-3 (B09 F153): EMBOSS reader - gap characters and other non-letters

`EmbossGapCharacters_Oracle_Tests` (`emboss_gap_oracle.tsv`, `generate_emboss_gap_oracle.py`): water (240), matcher (240), supermatcher (240) and dotmatcher (240) - the installed EMBOSS 6.6.0 `water, matcher, supermatcher, dotmatcher` on inputs that contain `- . ~ ? # + _`, digits, blanks, punctuation, `*` and non-ASCII letters equal `WaterAlign`, `LocalAlignAffineTopN(.., TopNTieOrder.Sim)`, `Supermatcher`, `DotMatcher.Compute` on every row, and the result of the cleaned input equals the result of the input. The reader rule (one shared `EmbossSequenceReader.Clean` / `CleanGapped`, Infrastructure): only `A-Z a-z *` survive in a plain sequence (`ajseqread.c seqAppend` + `ajseqtype.c ajSeqTypeCheckIn` + `ajstr.c ajStrRemoveGap`); coordinates and lengths count the cleaned sequence. Also `Clean_KeepsLettersAndStar_AndDropsEveryOtherCharacter` (all 256 code points) and `IsNucleotide_IsTheAjSeqIsNucRule_OfEveryPort`.

## Audit A4-1 (B09 F157): cancellation

`SequenceAligner_Cancellation_Tests` (`PreCancelled_A4_1_...`, one `CancelledWhileRunning_...` test per class) and the MCP `CancellationTokenTests` / `EmbossProgramToolsTests`: `DotMatcher.Compute`, matcher (`TopNTieOrder.Sim`) poll the ambient token (a mid-run cancel stops with `OperationCanceledException` inside the named poll) and the MCP tool takes the request token.

## Audit A4-2 (B09 F158): resource guards

`SequenceAligner_ResourceLimits_Tests`: `DotMatcher.Compute` (`DotMatcherMaxSegments`), `TopNTieOrder.Sim` (`TopNSimMaxVertices`, lazy vertex list: `maxAlignments = int.MaxValue` = every alignment) - documented limit constants (`SequenceAligner.ResourceLimits.cs`), the exception of the `SequenceAligner.ResourceLimit` style raised before the allocation in 64-bit arithmetic, and a low-allocation assertion per limit.

## Audit A4-9 (B09 F164): dotmatcher edge cases

| ID | Test | Locks |
|---|---|---|
| A4-9 | `DotMatcher_Tests.Compute_EqualsDotmatcher_OnTheEdgeCaseRows` | 720 rows of the installed dotmatcher (`dotmatcher_edge_oracle.tsv`): lower / mixed case, ambiguity codes, `*`, gap characters, default matrix (the type of the first sequence; a nucleic first sequence with a non-nucleic second one is refused), window >= length (zero terminator beyond the end, up to 8 residues) |
| A4-9 | `DotMatcher_Tests.WindowEqualToASequence_ScoresTheDiagonalsThatFit_AndALongerWindowReadsZerosBeyondTheEnd` | a 12 / 10 nt pair: w = 10 -> 1:1:10:10, w = 11 -> 1:1:11:11, w = 12 -> none; the nucleic / protein refusal |


## WP43 additions (A6-3 / A6-5)
| ID | Test | Locks |
|---|---|---|
| A6-3 | `SequenceAligner_MatcherGap_Tests` | `MatcherProgram` equals matcher 6.6.0 on 872 `-gapopen` / `-gapextend` / `-alternatives` rows, in particular `-gapopen` below `-gapextend`; reset warnings (`StandardError`), refusal of non-integers |
| A6-5 | `SequenceAligner_SupermatcherWordLength_Tests` | supermatcher `-wordlen` above the longest target: `Error: no kmers found`, error file, alignments (260 rows); `int.MaxValue` returns at once |

## Audit A2-28d (B09 F203): est2genome

| Item | Tests | Oracle | Status |
|---|---|---|---|
| A2-28d | `EmbossEst2genome_Oracle_Tests.Est2genome_OutfileEqualsTheProgram_AllRows` | 330 `est2genome` 6.6.0 runs (`est2genome_oracle.jsonl.gz`, `generate_est2genome_oracle.py`): the -outfile text equal, 19 crash / error rows (-nobest with one strand mode: SIGSEGV), the space warning | as specified |
| A2-28d | `EmbossEst2genome_Oracle_Tests.Est2genome_WorkedExample_...`, MCP `EmbossProgramToolsTests.Est2genome_...`, `ToolRegistrationContractTests`, `ToolDocsContractTests` | the three-exon example of the program | as specified |

## Audit A5-8c3e (B09 F241 - F243): feature readers, UFO, featcopy

| Item | Tests | Oracle | Status |
|---|---|---|---|
| A5-8c3e | `EmbossFeatureVocab_Oracle_Tests.VocabularyTables_EqualTheLibrary_AllEntries` | 2 406 entries of the ten vocabularies (`emboss_feature_vocab_oracle.tsv.gz`, `emboss_feature_vocab_dump.c`) | done |
| A5-8c3e | `EmbossFeatureReaders_Oracle_Tests.FeatureReaders_EqualFeatcopyDebug_AllRows` | 2 600 `featcopy -offormat debug` runs of the 11 input formats (`featread_oracle.jsonl.gz`, `generate_featread_oracle.py`) | done |
| A5-8c3e | `EmbossFeatureReaders_Oracle_Tests.FeatureFileOfASequence_EqualsSeqretDebug_AllRows` | 1 200 `seqret -feature -ufo1` runs (`featread_ufo_oracle.jsonl.gz`) | done |
| A5-8c3e | `EmbossDiffseq_Oracle_Tests.Diffseq_FeatureFilesOfTheSequences_EqualTheProgram_AllRows` | 450 `diffseq -ufo1 / -ufo2` runs (`diffseq_ufo_oracle.jsonl.gz`, `generate_diffseq_ufo_oracle.py`) | done |
| A5-8c3e | `EmbossFeatureReaders_Oracle_Tests.Featcopy_WritesEveryFormat_EqualsTheProgram_AllRows`, MCP `EmbossProgramToolsTests.Featcopy_...`, `Diffseq_FeaturesFromSeparateFiles_...`, `ToolRegistrationContractTests`, `ToolDocsContractTests` | 1 800 `featcopy` output-format runs (`featwrite_oracle.jsonl.gz`) | done |

## Audit A5-8c3h (B09 F255 - F263, WP58): USA forms of the sequence input, sequence databases and the feature programs

| Item | Tests | Oracle | Status |
|---|---|---|---|
| A5-8c3h (2d) | `EmbossUsa_Oracle_Tests.UsaForms_EqualSeqret_AllRows` | 1 500 `seqret -auto [-feature] -sequence SPEC ...` runs (`usa_oracle.jsonl.gz`, `generate_usa_oracle.py`): `[b:e:r]` ranges, `:b:e` / `:b..e:r`, `format::`, `@list` / `list::` list files, `file:id` / `file:field:value` / `file{a\|b}` / `%offset` queries, `-sbegin1 -send1 -sreverse1 -sid1 -sdbname1 -squick1 -iquery1 -ioffset1`; sequences, features, messages byte for byte (27 refused with proof) | done |
| A5-8c3h (2d) | MCP `EmbossProgramToolsTests.Seqret_ReadsTheUsaForms_AsTheProgram`, `ToolDocsContractTests`, `ToolRegistrationContractTests` | ranges, strands, list files, queries and `-sid` on the texts of the installed program | done |
| A5-8c3h (3) | `EmbossUsa_Oracle_Tests.SequenceDatabases_EqualSeqret_AllRows` | 600 `EMBOSSRC=work seqret ...` runs (`usa_db_oracle.jsonl.gz`, `generate_db_oracle.py`): databases of a configuration file with the method `direct`, `db:id` / `db-field:value` / `db:{a\|b}` / wildcard queries (14 refused with proof) | done |
| A5-8c3h (4) | `EmbossEntret_Oracle_Tests.Entret_EqualsTheProgram_AllRows`, MCP `EmbossProgramToolsTests.Entret_ReturnsTheTextOfTheEntries_AsTheProgram` | 756 `entret -auto -sequence SPEC ...` runs (`entret_oracle.jsonl.gz`, `generate_entret_oracle.py`): the text of every entry as the readers of 63 input formats store it (F257; 747 equal, 9 refused) | done |
| A5-8c3h (4) | `EmbossExtractfeat_Oracle_Tests.Extractfeat_EqualsTheProgram_AllRows`, MCP `Extractfeat_ExtractsTheFeatures_AsTheProgram` | 1 000 `extractfeat` runs (`extractfeat_oracle.jsonl.gz`, `generate_extractfeat_oracle.py`): filters, before / after, join, describe, remote parts, 25 output formats (F258; 923 equal, 77 refused) | done |
| A5-8c3h (4) | `EmbossCoderet_Oracle_Tests.Coderet_EqualsTheProgram_AllRows`, MCP `Coderet_WritesTheCodingSequences_AsTheProgram` | 1 000 `coderet` runs (`coderet_oracle.jsonl.gz`, `generate_coderet_oracle.py`): CDS, mRNA, translations, non-coding regions, four outputs (F259; 948 equal, 52 refused) | done |
| A5-8c3h (4) | `EmbossShowfeat_Oracle_Tests.Showfeat_EqualsTheProgram_AllRows`, MCP `Showfeat_DrawsTheFeatures_AsTheProgram` | 1 000 `showfeat` runs (`showfeat_oracle.jsonl.gz`, `generate_showfeat_oracle.py`): every option (F260; 911 equal, 89 refused) | done |
| A5-8c3h (4) | `EmbossUnionFeature_Oracle_Tests.UnionFeature_EqualsTheProgram_AllRows`, MCP `Union_JoinsTheFeatures_AsTheProgram` | 1 000 `union -feature -source` runs (`union_feature_oracle.jsonl.gz`, `generate_union_feature_oracle.py`): joined features, feature files and their names and messages, ranges (F261; 957 equal, 43 refused) | done |
| A5-8c3h (4) | `EmbossFeatreport_Oracle_Tests.Featreport_EqualsTheProgram_AllRows`, MCP `Featreport_ReportsTheFeatureTable_AsTheProgram` | 1 000 `featreport` runs (`featreport_oracle.jsonl.gz`, `generate_featreport_oracle.py`): 7 feature file formats, all 24 report formats (F262; 999 equal, 1 refused) | done |
| A5-8c3h (4) | `EmbossTwofeat_Oracle_Tests.Twofeat_EqualsTheProgram_AllRows`, MCP `Twofeat_ReportsThePairs_AsTheProgram`, `ToolDocsContractTests`, `ToolRegistrationContractTests` | 1 000 `twofeat` runs (`twofeat_oracle.jsonl.gz`, `generate_twofeat_oracle.py`): every option and report format, the fatal error of the opposite sense (F263; 962 equal, 38 refused) | done |
| A5-8c3h (2c-rest) | `EmbossSeqout_Oracle_Tests.SequenceWriters_EqualSeqret_AllRows` (ncbi, gifasta, nib, nibble, dasdna, das, asn1), MCP `EmbossProgramToolsTests.Seqret_WritesTheNcbiDasAndAsn1Formats_AsTheProgram` | 1 820 `seqret -osformat FMT ...` runs (`seqout_oracle_ncbi.jsonl.gz`, `seqout_oracle_nibble.jsonl.gz`, `seqout_oracle_das.jsonl.gz`, `seqout_oracle_asn1.jsonl.gz`, `generate_seqout_oracle.py`): 260 per format, all equal (F264) | done |
| A5-8c3h (2c-rest) | `EmbossSeqout_Oracle_Tests.SequenceWriters_EqualSeqret_AllRows` (fastq, fastq-sanger / illumina / solexa, sam, gff2, gff3 / gff, msf, clustal, aln, selex, phylip, phylipnon, phylip3, hennig86, mega, meganon, nexus, paup, nexusnon, paupnon, jackknifer, jackknifernon, treecon, scf, experiment, debug, bam, embl, em, emblnew, emold, emblold, genbank, gb, ddbj, refseq, genpept, refseqp, swiss, sw, swissprot, uniprot, swissnew, swnew, swold, swissold, swissprotold, swissprotnew), MCP `Seqret_WritesTheEntryQualityAndSetFormats_AsTheProgram`, `Seqret_WritesTheOutputFormats_AsTheProgram` | `seqout_oracle_*.jsonl.gz` (`generate_seqout_oracle.py`, valgrind + `MALLOC_PERTURB_` repeated runs): 260 runs per quality / entry writer name (fastq*, sam, scf, experiment, debug, bam, embl family, genbank family, swiss family), 230 per set-format name, 4 400 runs with the features inside the sequence file, 900 `-osdbname`, 275 `-osdirectory` runs; every row equal or refused with proof (F265 - F272) | done |
| A5-8c3h (d) | `EmbossExtractfeat_Oracle_Tests`, `EmbossCoderet_Oracle_Tests`, `EmbossUnionFeature_Oracle_Tests`, MCP `Extractfeat_ExtractsTheFeatures_AsTheProgram`, `Coderet_WritesTheCodingSequences_AsTheProgram` | 1 100 `extractfeat` runs (168 with `-ossingle2`, 176 with `-osdirectory2`), 1 100 `coderet` runs (374 with `-ossingleN`, 395 with `-osdirectoryN`, `-osextensionN`, `-osnameN`), 1 000 `union -feature` runs, every sequence writer (F273) | done |
| A5-8c3h (2b-rest) | `EmbossNibbleReader_Tests` | `seqret -sformat1 nibble / nib` on FASTA text and on a nibble file: "Unable to read sequence" (F274) | done |
| A5-8c3h (c, d-rest) | `EmbossDbIndex_Oracle_Tests` (`MethodEmboss`, `MethodEmblcd`), `EmbossUsa_Oracle_Tests.ConfigurationBlocks`, `EmbossCoderet_Oracle_Tests`, `EmbossUnionFeature_Oracle_Tests`, `EmbossFeatreport_Oracle_Tests`, `EmbossTwofeat_Oracle_Tests`, MCP `Seqret_ReadsTheIndexOfTheMethodEmboss`, `Union_WritesTheFileOfOutputSingle`, `Featreport_NamesTheReportFile` | 1 200 `emboss` and 1 200 `emblcd` index runs of `seqret` / `entret`, 400 configuration-block runs, 1 100 `coderet` (with `-sreverse1`), 1 200 `union` (output qualifiers, experiment / debug / bam), 1 000 `featreport` and 1 000 `twofeat` (report file qualifiers, features qualifiers), 220 experiment writer runs with features (F275 - F279) | done; gcg / embossgcg / blast / app / DRCAT / abi, bam, scf, das, ensembl, `CO` lines registered |
