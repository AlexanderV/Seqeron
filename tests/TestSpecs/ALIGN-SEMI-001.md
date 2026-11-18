# Test Specification: ALIGN-SEMI-001

**Test Unit ID:** ALIGN-SEMI-001
**Area:** Alignment
**Algorithm:** Semi-Global Alignment (Fitting / Query-in-Reference)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-03-07

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Sequence alignment | https://en.wikipedia.org/wiki/Sequence_alignment | 2026-02-01 |
| Wikipedia: Needleman–Wunsch algorithm | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | 2026-02-01 |
| Rosalind: Finding a Motif with Modifications (SIMS) | https://rosalind.info/problems/sims/ | 2026-02-01 |
| Rosalind: Semiglobal Alignment (SMGB) | https://rosalind.info/problems/smgb/ | 2026-02-01 |
| Brudno et al. (2003) | doi:10.1093/bioinformatics/btg1005 | 2026-02-01 |

### 1.2 Key Evidence Points

1. Semi-global alignment (also "glocal") searches for the best partial alignment, allowing one or both ends to be unaligned without penalty. (Wikipedia: Sequence alignment)
2. Primary use case: aligning a short query to a longer reference where the query should be globally aligned but the reference may have free end gaps. (Wikipedia: Sequence alignment)
3. The implementation uses the **fitting alignment** variant, corresponding to Rosalind problem SIMS: "an alignment of a substring of s against all of t." (Rosalind: SIMS)
4. Algorithm is a modification of Needleman–Wunsch: first row = 0, first column = $d \cdot i$, traceback from $\max_j F_{m,j}$. (Wikipedia: NW algorithm; Rosalind: SIMS)
5. No zero floor in recurrence (unlike local/Smith–Waterman); scores can be negative. (Wikipedia: NW algorithm)
6. Complexity: $O(m \cdot n)$ time and space. (Wikipedia: NW algorithm)

### 1.3 Documented Corner Cases

| Case | Evidence Source | Expected Behavior |
|------|----------------|-------------------|
| Short query in long reference | Rosalind SIMS; Wikipedia (Sequence alignment) | Query fully aligned, reference has free leading/trailing gaps |
| Query equals reference | NW global = fitting when $m = n$ | Full alignment, score = $m \times \text{match}$ |
| All mismatches | NW recurrence (no zero floor) | Negative score |
| Mixed matches and mismatches | NW recurrence | Score = $\sum(\text{match scores}) + \sum(\text{mismatch scores})$ |
| Gap in optimal alignment | NW recurrence (up/left moves) | Score includes gap penalty |
| Null sequence | .NET API convention | `ArgumentNullException` |

### 1.4 Known Failure Modes / Pitfalls

1. Multiple equally optimal alignments may exist (as with all DP alignment algorithms).
2. Score must be taken from $F_{m, \text{maxJ}}$, not $F_{m,n}$ — these differ when query aligns in the middle of reference.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `SemiGlobalAlign(DnaSequence\|string, …, ScoringMatrix?)` | SequenceAligner | **Canonical** | Query-in-reference fitting alignment, linear gaps (GapExtend per position; GapOpen ignored); string overload added 2026-11 (null → `ArgumentNullException`, empty not special-cased) |
| `SemiGlobalAlignAffine(DnaSequence\|string, …, ScoringMatrix?)` | SequenceAligner | Variant | Same fitting variant, affine gaps (Gotoh; gap of length k = GapOpen + k·GapExtend) via the shared `AffineAlignCore` |

**Reference equivalents (2026-09 review, confirmed on 2428 random/edge cases):** Biopython 1.88
`PairwiseAligner(mode="global", end_insertion_score=0)` with target = sequence 1 (linear: open = extend = GapExtend;
affine: open_gap_score = GapOpen + GapExtend, extend_gap_score = GapExtend) and parasail 1.3.4 `sg_dx`
(s1 = sequence 1; affine open = −(GapOpen + GapExtend), extend = −GapExtend). Only end gaps *in sequence 1*
(overhanging reference residues) are free; a query overhang is charged (ACGTAA / CCACGT → 2, not 4).
Ties: smallest end column j of the last row (j = 0 allowed, i.e. query entirely against gaps), traceback
diagonal > up > left (affine: M > X > Y). Aligned strings span both full inputs (free reference flanks shown
against '-'). Coordinates (2026-11, A1-4): Start1/End1 = 0 and length − 1 (query aligned in full);
Start2/End2 = the **fitted interval** of the reference (first/last reference index consumed between the free end gaps; an
empty interval — empty query or query aligned entirely against gaps — has End2 = Start2 − 1). Biopython: leading/trailing
gap-in-sequence-1 columns of `alignment[0]` (the span of `aligned`); parasail: `end_ref` of `sg_dx_trace`. Verified on 4090
random/edge cases (1382 semi-global rows): reported interval ∈ the set of intervals of all Biopython co-optimal alignments
(`PairwiseAligner.align`, 823 fully enumerated; 32 capped), parasail score equal, interval global score = fitting score.

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | AlignmentType = SemiGlobal | Yes | Implementation contract |
| INV-2 | Aligned sequences have equal length | Yes | Alignment definition (all types) |
| INV-3 | `RemoveGaps(aligned1) == query` | Yes | Fitting alignment: query fully aligned (Rosalind SIMS) |
| INV-4 | Reference part of the scored region (columns between the first and last query residue) is a substring of the reference; `RemoveGaps(aligned2) == reference` (free flanks shown against gaps) | Yes | Fitting alignment: free reference end gaps |
| INV-5 | Score = $\max_j F_{m,j}$ | Yes | Fitting alignment traceback from max of last row |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Short query embedded in long reference | query="ATGC", ref="AAAATGCAAA" | Score=4, AlignmentType=SemiGlobal | Hand-computed DP; Rosalind SIMS |
| M2 | AlignmentType is SemiGlobal | query="ACGT", ref="ACGT" (default scoring) | `result.AlignmentType == SemiGlobal` | Implementation contract |
| M3 | Aligned sequences have equal length | query="ATGC", ref="AAAATGCAAA" | `len(aligned1) == len(aligned2)` | Alignment definition |
| M4 | Query fully represented + exact score | query="ATGC", ref="ATGCAAAA" | `RemoveGaps(aligned1) == "ATGC"`, Score=4 | INV-3; hand-computed DP |
| M5 | Reference is substring after gap removal | query="ATGC", ref="AAAATGCAAA" | `RemoveGaps(aligned2)` is substring of ref | INV-4 |
| M7 | Null seq1 throws | seq1=null | `ArgumentNullException` | .NET API convention |
| M8 | Null seq2 throws | seq2=null | `ArgumentNullException` | .NET API convention |

*M6 removed: was duplicate of M4 (same inputs, subset of assertions).*

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | Identical sequences | query=ref="ATGCATGC" | Score=8 (exact, 8×match) | NW: identical = fitting when m=n |
| S2 | Query at start of reference | query="ATG", ref="ATGCCCCC" | Score=3 (exact, 3×match) | Hand-computed DP |
| S3 | Query at end of reference | query="CCC", ref="ATGCCC" | Score=3 (exact, 3×match) | Hand-computed DP |
| S4 | Custom scoring in fitting context | query="ATGC", ref="AATGCCC", Match=5,Mm=-3,Gap=-2 | Score=20 (exact, 4×5) | Hand-computed DP |

### 4.3 Additional Regression / Coverage Tests

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| NEG | All mismatches → exact negative score | query="AAAA", ref="CCCC" | Score=-4; tie at j=0 → `AAAA----`/`----CCCC` (4×gap ties 4×mismatch) | NW recurrence: no zero floor; Biopython co-optimal set |
| MAX | Score is max of last row, not bottom-right | query="ATG", ref="ATGCCC" | Score=3 (not 0) | Fitting alignment: $\max_j F_{m,j}$ |
| OFS | Match with offset — exact score | query="ACG", ref="AACGG" | Score=3 | Hand-computed DP |
| INV | All invariants validated in one test | query="GCATGCG", ref="AAAGCATGCGAAA" | Score=7 (7×match), INV-1..5 | Hand-computed DP |
| MIX | Mixed matches and mismatches | query="AGT", ref="AAACTAAA" | Score=1 (2×match + 1×mismatch) | Hand-computed DP |
| GAP | Gap in optimal alignment | query="ACGT", ref="AGT" | Score=2 (3×match + 1×gap), gap in aligned ref | Hand-computed DP |
| REF | Biopython/parasail reference alignments (linear) | GATTACA/CCGATTTACAGG (2,−3,d=−2) → 12 `--GA-TTACA--`, Start2/End2 = 2/9; ACGTAA/CCACGT → 2 `--ACGTAA`/`CCACGT--`, 2/5; GGGACGT/ACGTCCC → 1, 0/3 | Exact strings + positions | Biopython 1.88, parasail sg_dx_trace |
| COORD-1 | Fitted interval (A1-4) vs Biopython `align` set / parasail `end_ref` | 8 unique-interval cases (linear + affine), e.g. ACGTAC/TTTTTTACGTACTTTTTT → 6..11; AGGGCT/TTAGGCTT (2,−3,−5,−2) → 2..6 | Exact | Biopython 1.88, parasail 1.3.4 |
| COORD-2 | Tied fittings report one co-optimal interval | AAAA/CCCC (1,−5,gap −1) → empty (0,−1) of {(0,−1)…(4,3)}; TTTTACGTA/GGGGACGTAGGGG → (0,8) of {(0..4,8)}; AC/CACACAC → (1,2) of {(1,2),(3,4),(5,6)} | Exact | Biopython 1.88 full enumeration |
| COORD-3 | Empty query / empty reference | ""/ACGT → coordinates (0,−1,0,−1); ACG/"" → (0,2,0,−1) | Exact | DP borders |
| COORD-4 | Brute-force definition oracle | 600 random pairs (linear + affine, 6 schemes): score = max over all reference substrings of the global alignment score; the reported interval attains it; equals the leading/trailing gap span of the query row | Exact | Independent of the DP |
| AFF-1 | Affine long gap | ACGTACGTTTTTTACGTACGT / GGGGACGTACGTACGTACGTGGGG, HighIdentityDna | 65 (linear 75) | Biopython = parasail |
| AFF-2 | Affine opening cost prevents gap | ACGTACGTACGT / TTACGTACGTTTTTACGTAA, (1,−1,−5,−1) | 6 ungapped (linear 8) | Biopython = parasail |
| AFF-3 | Affine reference values | ACGTAA/CCACGT → 0; AAAA/CCCC → −4 `AAAA`/`CCCC`; GATTACA/CCGATTTACAGG (2,−3,−5,−2) → 7 | Exact strings | Biopython = parasail |
| AFF-4 | GapOpen = 0 ⇒ affine = linear | GATTACA/CCGATTTACAGG (2,−3,0,−2) | 12 | Model reduction |
| AFF-5 | Empty inputs | ""/ACGT → 0 `----`/`ACGT`; ACG/"" → −5 (linear −3) | Exact | DP borders |

---

## 5. Coverage Classification

### 5.1 Discovery Summary

- **Canonical:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_SemiGlobalAlign_Tests.cs` — 17 tests
- **Property:** `tests/Seqeron/Seqeron.Genomics.Tests/Properties/AlignmentProperties.cs` — 1 semi-global test
- **Legacy:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAlignerTests.cs` — 0 (moved out, comments only)

### 5.2 Coverage Classification

#### Canonical (`SequenceAligner_SemiGlobalAlign_Tests.cs`) — 17 test methods

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `SemiGlobalAlign_ShortQueryInLongReference_FindsMatch` | M1 | ✅ |
| 2 | `SemiGlobalAlign_ValidInput_AlignmentTypeIsSemiGlobal` | M2 | ✅ |
| 3 | `SemiGlobalAlign_ValidInput_AlignedSequencesHaveEqualLength` | M3 | ✅ |
| 4 | `SemiGlobalAlign_QueryFullyRepresented` | M4 | ✅ |
| 5 | `SemiGlobalAlign_ShortQueryInReference_ReferenceIsSubstring` | M5 | ✅ |
| 6 | `SemiGlobalAlign_NullSequence1_ThrowsArgumentNullException` | M7 | ✅ |
| 7 | `SemiGlobalAlign_NullSequence2_ThrowsArgumentNullException` | M8 | ✅ |
| 8 | `SemiGlobalAlign_IdenticalSequences_FullAlignmentMaxScore` | S1 | ✅ |
| 9 | `SemiGlobalAlign_QueryAtStart_ExactScore` | S2 | ✅ |
| 10 | `SemiGlobalAlign_QueryAtEnd_ExactScore` | S3 | ✅ |
| 11 | `SemiGlobalAlign_CustomScoring_ExactScore` | S4 | ✅ |
| 12 | `SemiGlobalAlign_ValidatesAllInvariants` | INV | ✅ |
| 13 | `SemiGlobalAlign_AllMismatches_NegativeExactScore` | NEG | ✅ |
| 14 | `SemiGlobalAlign_ScoreIsMaxOfLastRow_NotBottomRight` | MAX | ✅ |
| 15 | `SemiGlobalAlign_MatchWithOffset_ExactScore` | OFS | ✅ |
| 16 | `SemiGlobalAlign_MixedMatchMismatch_ExactScore` | MIX | ✅ |
| 17 | `SemiGlobalAlign_GapInAlignment_ExactScore` | GAP | ✅ |

#### Property (`AlignmentProperties.cs`) — 1 semi-global test

| # | Test Method | Status | Notes |
|---|-------------|--------|-------|
| 1 | `SemiGlobalAlign_AlignedSequences_HaveEqualLength` | ✅ | Structural invariant, different inputs |

#### Classification Summary

- ✅ Covered: 17 canonical + 1 property = 18 total
- ❌ Missing: 0
- ⚠ Weak: 0
- 🔁 Duplicate: 0

### 5.3 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M6 | 🔁 Duplicate | Removed (subset of M4, same inputs) | ✅ Done |
| 2 | S4 | ⚠ Weak | Rewritten: fitting context with non-identical seqs | ✅ Done |
| 3 | INV | ⚠ Weak | Fixed erroneous comment ("6×match+1×mm" → "7×match") | ✅ Done |
| 4 | NEG | ⚠ Weak | Strengthened: `Is.LessThan(0)` → `Is.EqualTo(-4)` | ✅ Done |
| 5 | MIX | ❌ Missing | Added: mixed match/mismatch test, Score=1 | ✅ Done |
| 6 | GAP | ❌ Missing | Added: gap in alignment test, Score=2 | ✅ Done |

**Total items:** 6
**✅ Done:** 6 | **⛔ Blocked:** 0 | **Remaining:** 0

---

## 6. Design Decisions

| Decision | Rationale | Source |
|----------|-----------|--------|
| Query-in-reference (fitting) variant | Corresponds to Rosalind SIMS; most common bioinformatics use case | Rosalind SIMS; Wikipedia |
| Linear gap cost (GapExtend only) in `SemiGlobalAlign` | GapOpen unused; affine gaps are provided by `SemiGlobalAlignAffine` | NW linear gap model; Gotoh (1982) |
| Overlap and all-ends-free variants (A1-3, F44) | Implemented as one `EndGaps` option on the same programs (no second DP); fitting = `FreeStart2 \| FreeEnd2` stays the default | Biopython `left/right_insertion/deletion_score`; parasail `sg_*`; EMBOSS needle `-endweight` |

---

## 7. Validation Checklist

- [x] All MUST tests implemented and passing (M1–M5, M7–M8)
- [x] All SHOULD tests implemented and passing (S1–S4)
- [x] Regression tests for critical bug fixes (NEG, MAX, OFS, INV, MIX, GAP)
- [x] Tests are deterministic
- [x] Invariants verified with `Assert.Multiple` (INV test)
- [x] No duplicate tests across files (M6 removed)
- [x] Evidence documented for all test rationale — no ASSUMPTION tags
- [x] All scores hand-computed from DP matrices, not taken from implementation output
- [x] External sources: Wikipedia, Rosalind SIMS/SMGB
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
- [x] Mixed match/mismatch path tested (MIX)
- [x] Gap penalty path tested (GAP)
- [x] Negative score exact value tested (NEG)

## 2026-11 audit round 1 addition (A1-2, B09 F41)

Matrix-scored `SemiGlobalAlign`/`SemiGlobalAlignAffine`: see the table in ALIGN-GLOBAL-001.md (semi rows = Biopython `end_insertion_score = 0` = parasail `sg_dx`; fitted-interval re-scoring and coordinates checked in `Oracle_Alignment_RescoresToScore_AndReproducesInputs`).

## 2026-11 audit round 1 addition (A1-3, B09 F44): ends-free / overlap variants (`EndGaps`)

**Scope.** `SemiGlobalAlign(…, EndGaps)` / `SemiGlobalAlignAffine(…, EndGaps)` for `string`, `DnaSequence` (+ `ScoringMatrix`) and `string`, `ProteinSequence`, `DnaSequence` (+ `MatrixScoring`): each of the four end gaps free or penalised (16 policies: global, fitting in both directions, overlap in both orientations, all free = EMBOSS needle default, single flags), linear and affine. Coordinates = consumed intervals of both sequences (AlignmentResult docs). Tests: `SequenceAligner_EndGaps_Tests` (14) + MCP `EndGapsAlignTests` (8 cases).

**Reference equivalents.** Biopython 1.88 `left_deletion_score` = `FreeStart1`, `right_deletion_score` = `FreeEnd1`, `left_insertion_score` = `FreeStart2`, `right_insertion_score` = `FreeEnd2` (sequence 1 = target); parasail `_qb`/`_qe`/`_db`/`_de`; EMBOSS needle default `-endweight N` = `All`.

| ID | Test | Input / data | Expected | Evidence |
|----|------|--------------|----------|----------|
| EG-1 | `Oracle_Score_EqualsBiopython_AllRows` | 3 590 rows (2 240 scalar, 1 350 BLOSUM62/EDNAFULL, fractional gaps), all 16 policies, linear + affine | score = Biopython `align().score` (1e-6) | endgaps_oracle.tsv |
| EG-2 | `Oracle_Score_EqualsParasail_ExactVariants_AndNeverAboveOptimumForTheOthers` | 1 680 integer rows | equal for `nw`, `sg_qb`, `sg_db`, `sg_qx`, `sg_dx`, `sg_qb_db`; parasail never above the optimum for the variants with a free end flag (23 lower rows: parasail skips the border end cells (m,0)/(0,n)) | parasail 1.3.4 |
| EG-3 | `Oracle_Emboss_NeedleDefaultEndweight_IsEndGapsAll` | 287 needle default rows (`-endweight N`) | >= 95 % equal (277 of 287), never above the optimum (embaln.c opens a last-row/column end gap only from the match matrix) | EMBOSS 6.6.0 |
| EG-4 | `Oracle_Alignment_IsOneOfBiopythonCoOptimalAlignments` | 2 050 rows with <= 40 Biopython co-optimal alignments | returned rows in the (gap-order-normalised) set | Biopython `align` |
| EG-5 | `Oracle_Alignment_RescoresToScore_AndCoordinatesAreTheConsumedIntervals` | all rows | rows reproduce inputs; independent re-scoring with the policy = score; flank columns match the coordinates; no overhang on a penalised end | independent rescorer |
| EG-6 | `BruteForce_EveryAlignment_MaximumEqualsDynamicProgram_All16Policies` | 120 pairs (<= 5 residues), every alignment enumerated, 16 policies x linear/affine | maximum = DP score | exhaustive definition |
| EG-7 | `Fitting_IsTheClassicSemiGlobal_AndNone_IsGlobal_RowForRow`, `MatrixScoring_FittingAndNone_Equal…` | 300 x 6 schemes (+ DNA overloads, 200 protein pairs) | `Fitting` = `SemiGlobalAlign`, `None` = `GlobalAlign(Affine)` incl. rows, coordinates, tags | additivity |
| EG-8 | `Overlap_SuffixOfSeq1_PrefixOfSeq2_RowsAndCoordinates` | GGGGACGTAC/ACGTACTTTT (1,−1,−1) | `GGGGACGTAC----`/`----ACGTACTTTT`, 6, (4,9,0,5); swapped + other orientation (0,5,4,9) | Biopython (unique alignment) |
| EG-9 | `FittingSequence2_IsFittingOfTheSwappedPair`, `Score_IsMonotoneInTheFreeFlags` | 200 / 120 pairs | symmetry of scores; freeing a flag never lowers the score | model |
| EG-10 | `ZeroLengthOverlap_ScoresZero_BiopythonAlignScore`, `EmptyAndNull_Inputs_FollowTheAlignerContract` | AC/T overlap → `AC-`/`--T` 0, (2,1,0,−1); empty sides; null | exact; documents Biopython's `score()` vs `align()` disagreement (25 oracle rows) | Biopython 1.88 |

## 2026-11 audit round 1 addition (A1-5, B09 F45): general gap-cost functions

`SemiGlobalAlignGeneralGap(…, Func<int,double> gapScore, …, EndGaps)` (ends-free with Waterman-Smith-Beyer gaps; Biopython ignores end-gap properties for a callable `gap_score`, so the policy is locked by the independent reference and an exhaustive enumeration): see the table in ALIGN-GLOBAL-001.md (A1-5 section) and `SequenceAligner_GeneralGap_Tests`.
