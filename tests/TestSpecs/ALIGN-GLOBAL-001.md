# Test Specification: ALIGN-GLOBAL-001

**Test Unit ID:** ALIGN-GLOBAL-001
**Area:** Alignment
**Algorithm:** Global Alignment (Needleman–Wunsch)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-03-06

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Needleman–Wunsch algorithm | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | 2026-03-06 |
| Wikipedia: Sequence alignment (global alignment) | https://en.wikipedia.org/wiki/Sequence_alignment | 2026-02-01 |

### 1.2 Key Evidence Points

1. **Algorithm**: Needleman–Wunsch uses dynamic programming to compute the optimal end-to-end alignment score for two sequences. (Wikipedia: Needleman–Wunsch algorithm)
2. **Scoring**: The alignment score is the sum of per-position scores for match, mismatch, or gap. (Wikipedia: "Choosing a scoring system")
3. **Linear gap penalty**: The standard NW uses a single gap penalty `d` applied uniformly per gap position. Initialization: `F(i,0) = d·i`, `F(0,j) = d·j`. Recurrence: `F(i,j) = max(F(i−1,j−1)+S(Aᵢ,Bⱼ), F(i−1,j)+d, F(i,j−1)+d)`. (Wikipedia: "Advanced presentation of algorithm")
4. **Traceback**: Follows the max choices from `F(n,m)` back to `F(0,0)`, inserting gaps on vertical/horizontal moves. (Wikipedia: "Tracing arrows back to origin")
5. **Example dataset**: GCATGCG vs GATTACA with match=+1, mismatch/indel=−1 yields optimal score 0. One optimal alignment: `GCATG-CG` / `G-ATTACA`. (Wikipedia: "Choosing a scoring system")
6. **Multiple optimal alignments**: When several traceback paths share the same score, each represents an equally valid alignment; which one is returned depends on implementation choices. (Wikipedia: "Tracing arrows back to origin")

### 1.3 Documented Corner Cases (from sources)

1. **Completely different sequences**: All mismatches and gaps; score is negative with negative penalties. (Wikipedia: Sequence alignment)
2. **Identical sequences**: Perfect diagonal alignment, score = n × match. (Trivially follows from the recurrence.)
3. **Scoring system dependence**: Different match/mismatch/gap values yield different optimal alignments. (Wikipedia: "Scoring systems")

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GlobalAlign(seq1, seq2, scoring)` | SequenceAligner | **Canonical** | Needleman–Wunsch with linear gap penalty |
| `GlobalAlign(string, string, scoring)` | SequenceAligner | Delegate | Wrapper: normalizes to uppercase; null → `ArgumentNullException`; empty side → all-gap NW alignment (same as every other overload) |
| `GlobalAlign(string, string, scoring, CancellationToken, IProgress?)` | SequenceAligner | Delegate | Cancellation wrapper; shares the same linear NW core (2026-09 review) |
| `GlobalAlign(DnaSequence, DnaSequence, scoring, CancellationToken, IProgress?)` | SequenceAligner | Delegate | Same core as the typed overload (empty typed input → all-gap NW border alignment, like the non-cancellation typed overload) |
| `GlobalAlignAffine(seq1, seq2, scoring)` / `(string, string, scoring)` | SequenceAligner | **Canonical (affine)** | Gotoh (1982) three-state DP, Flouri et al. (2015) border initialization; gap of length k = GapOpen + k·GapExtend |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Aligned sequences have equal length | Yes | Wikipedia: Sequence alignment (global alignment pads both sequences) |
| INV-2 | Removing gaps from aligned sequences yields the original sequences | Yes | Wikipedia: Needleman–Wunsch traceback rules |
| INV-3 | Alignment score equals sum of match/mismatch/gap scores across aligned positions | Yes | Wikipedia: "Choosing a scoring system" — score is the sum of per-position pairings |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Wikipedia example: score is 0 | seq1=GCATGCG, seq2=GATTACA, match=+1, mismatch=−1, gap=−1 | Score = 0 | Wikipedia: "Choosing a scoring system" — explicitly shows score 0 |
| M2 | Wikipedia example: all invariants hold | same as M1 | INV-1, INV-2, INV-3 hold | Wikipedia: algorithm definition + scoring rules |
| M3 | Unequal lengths (short vs long): correct linear border score | seq1=T, seq2=ACGT, gap=−1 | Score = −2 | Wikipedia: F(0,j) = d·j from pseudocode |
| M4 | Unequal lengths (long vs short): correct linear border score | seq1=ACGT, seq2=T, gap=−1 | Score = −2 | Wikipedia: F(i,0) = d·i from pseudocode |
| M5 | Completely different equal-length sequences | seq1=AAAA, seq2=TTTT, match=+1, mismatch=−1, gap=−1 | Score = −4 (4 mismatches) | Wikipedia: Sequence alignment corner case |
| M6 | Identical sequences: perfect alignment | seq1=ACGTACGT, seq2=ACGTACGT, match=+1, gap=−1 | Score = 8, no gaps | Follows from recurrence — diagonal dominates |
| M7 | Single deletion: gap penalty applied correctly | seq1=ACGT, seq2=AGT, match=+1, gap=−1 | Score = 2 | Linear gap model: 3 matches + 1 gap |
| M8 | Different scoring matrix yields correct result | seq1=ACGT, seq2=ACGT, match=+5, gap=−1 | Score = 20 | Wikipedia: "Basic scoring schemes" — arbitrary values |
| M9 | String overload matches DnaSequence overload | same as M1 | Same score and alignment | API contract |
| M10 | Single deletion: statistics exact values | seq1=ACGT, seq2=AGT, match=+1, gap=−1 | Matches=3, Mismatches=0, Gaps=1, Identity=75.0%, GapPercent=25.0% | Deterministic traceback: ACGT/A-GT |
| M11 | Score symmetry: reversed inputs produce same score | seq1=GCATGCG, seq2=GATTACA | score(A,B) = score(B,A) | NW recurrence symmetry (S is symmetric, d is uniform) |

### 4.1b Affine-gap (Gotoh) Tests — added 2026-09 review

Reference values: Biopython 1.88 `PairwiseAligner(mode='global', open_gap_score=GapOpen+GapExtend, extend_gap_score=GapExtend)` and parasail 1.3.4 `nw_trace(open=-(GapOpen+GapExtend), extend=-GapExtend)` — identical on 2008 random + edge cases.

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| G1 | BlastDna preset = NCBI blastn defaults | — | (2, −3, −5, −2) | BLAST+ User Manual (NBK279684), blastn options |
| G2 | One long internal gap charged one opening | ACGTTTTACGT / ACGTACGT, BlastDna | 5 | Biopython = parasail |
| G3 | Reference scores | TTTTTTTTTTACGT/ACGT BlastDna → −17; GCATGCG/GATTACA BlastDna → −6; ACGTAAAAAAAAACGT/ACGTCGT (5,−4,−10,−1) → 16; AC/GC (1,−10,−1,−1) → −3 | as listed | Biopython = parasail |
| G4 | Flouri initialization trap | A / CAC, (1,−10,−5,−1) | −11, `-A-`/`CAC` (Gotoh's original border gives −9) | Flouri et al. 2015; Biopython = parasail |
| G5 | GapOpen = 0 reduces to linear NW | Wikipedia example | 0 | w(k) = k·e |
| G6 | Empty side | "" / ACGT, BlastDna | −13 (o + n·e); both empty → score 0, empty rows, coordinates 0/0/−1/−1 | Gotoh border |
| G7 | Overload agreement | DnaSequence vs lower-case string | equal | API contract |

### 4.2 API Contract Tests

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| A1 | Empty side = all-gap Needleman-Wunsch alignment, every overload (string, DnaSequence, both cancellable) | ""/GATTACA → `-------`/`GATTACA` −7; GATTACA/"" → −7; ""/"" → score 0, "", "", coordinates 0/0/−1/−1; ""/A → −1; CCAC/"" with gap −2 → −8 | Exact | NW border F(0,j) = d·j (Needleman & Wunsch 1970; Durbin et al. 1998 §2.3). Biopython 1.88 / parasail 1.3.4 / EMBOSS needle 6.6 / scikit-bio 0.7 reject zero-length input (no oracle); non-empty cases of the same code path are oracle-verified (4090-row differential run, 0 score mismatches) |
| A2 | Null DnaSequence / null string throws ArgumentNullException (all string, DnaSequence, cancellable, affine overloads) | seq1=null, seq2=GATTACA | ArgumentNullException | API contract — null guard (replaces the former silent `AlignmentResult.Empty` of the string overloads) |

---

## 5. Test File

**Canonical file:** `Seqeron.Genomics.Tests/SequenceAligner_GlobalAlign_Tests.cs`

| Test Method | Spec ID |
|-------------|---------|
| `WikipediaExample_OptimalScore_IsZero` | M1 |
| `WikipediaExample_AllInvariants_Hold` | M2 |
| `UnequalLengths_ShortVsLong_CorrectLinearBorderScore` | M3 |
| `UnequalLengths_LongVsShort_CorrectLinearBorderScore` | M4 |
| `CompletelyDifferent_EqualLength_AllMismatchScore` | M5 |
| `IdenticalSequences_PerfectScore_NoGaps` | M6 |
| `SingleDeletion_CorrectGapPenalty` | M7 |
| `DifferentScoringMatrix_HighMatchReward_CorrectScore` | M8 |
| `StringOverload_ProducesSameResultAsDnaSequenceOverload` | M9 |
| `SingleDeletion_Statistics_ExactValues` | M10 |
| `ScoreSymmetry_ReversedInputs_SameScore` | M11 |
| `EmptyInput_AllOverloads_AreTheAllGapNeedlemanWunschAlignment` ×4, `EmptyInput_LinearScore_IsGapTimesLength` ×3 | A1 |
| `NullDnaSequence_ThrowsArgumentNullException`, `NullString_ThrowsArgumentNullException_OnEveryStringOverload` | A2 |

---

## 6. Deviations and Assumptions

**None.** All tests and implementation are grounded in the authoritative sources listed in §1.1.

- The implementation uses the standard Needleman–Wunsch linear gap penalty model exactly as described in the Wikipedia pseudocode.
- `ScoringMatrix.GapExtend` serves as the linear gap penalty `d`. `ScoringMatrix.GapOpen` is not used by `GlobalAlign`; the affine model (gap of length k = GapOpen + k·GapExtend, Gotoh 1982 / BLAST convention — GapOpen excludes the first extension) is `GlobalAlignAffine`.
- When multiple optimal alignments exist, the implementation returns one deterministically; this is explicitly supported by Wikipedia ("more than one choice may have the same value, leading to alternative optimal alignments").
- Empty-input handling and null-argument validation are API contract behaviors, not algorithm-level specifications.

---

## 7. Coverage Classification

### 7.1 Coverage Classification Summary

| Coverage Area | Status | Tests | Notes |
|---------------|--------|-------|-------|
| Score correctness (known examples) | ✅ Covered | M1, M3–M8 | 7 tests with exact scores derived from NW theory |
| Structural invariants (INV-1,2,3) | ✅ Covered | M2, M3, M4, M7 | Verified across Wikipedia + border + deletion cases |
| Statistics exact values | ✅ Covered | M10 | Exact match/mismatch/gap counts + identity/similarity |
| Score symmetry | ✅ Covered | M11 | score(A,B) == score(B,A) |
| API contract (string overload) | ✅ Covered | M9 | Same result as DnaSequence overload |
| API contract (empty/null) | ✅ Covered | A1, A2 | Guard behavior |
| CancellationToken path | ✅ Covered | PerformanceExtensionsTests (strengthened) + AlignmentProperties (equivalence) | Score + invariants + equivalence with standard overload |
| Performance regression | ✅ Covered | PerformanceRegressionTests | 500bp completion time |
| Snapshot regression | ✅ Covered | AlignmentSnapshotTests | Full output golden master |
| Property: Similarity+GapPercent=100 | ✅ Covered | AlignmentProperties | Structural percentage invariants |
| Property: Matches+Mismatches+Gaps=Length | ✅ Covered | AlignmentProperties | Statistics decomposition |

### 7.2 Actions Taken

| Action | Test | Rationale |
|--------|------|-----------|
| 🔁 **Removed** | `AlignmentProperties.GlobalAlign_AlignedSequences_HaveEqualLength` | Duplicate of INV-1 tested in M2, M3, M4, M7 |
| 🔁 **Removed** | `AlignmentProperties.GlobalAlign_IdenticalSequences_MaxScore` | Duplicate of M6 (`Score > 0` → M6 asserts `Score == 8`) |
| ⚠ **Strengthened** | `AlignmentProperties.Statistics_Identity_InRange` → `Statistics_PercentageFields_SatisfyInvariants` | Was `Identity ∈ [0,100]` only; now verifies 5 structural invariants across all percentage fields |
| ⚠ **Strengthened** | `PerformanceExtensionsTests.GlobalAlign_WithCancellation_CompletesNormally` | Was `result ≠ null`; now verifies Score=8, no gaps, AlignmentType, both aligned sequences |
| ❌ **Added** | `SequenceAligner_GlobalAlign_Tests.SingleDeletion_Statistics_ExactValues` (M10) | Missing: no test verified exact statistics counts |
| ❌ **Added** | `SequenceAligner_GlobalAlign_Tests.ScoreSymmetry_ReversedInputs_SameScore` (M11) | Missing: NW symmetry property not tested |
| ❌ **Added** | `AlignmentProperties.GlobalAlign_ScoreSymmetry_ReversedInputs` | Score symmetry as structural property |
| ❌ **Added** | `AlignmentProperties.GlobalAlign_CancellationOverload_SameResultAsStandard` | CancellationToken path equivalence with standard overload |

## 2026-11 audit round 1 addition (A1-2, B09 F41): substitution-matrix scoring

Tests: `SequenceAligner_MatrixScoring_Tests` (14) and `Seqeron.Mcp.Alignment.Tests.MatrixAlignTests`. Oracle table `TestData/Alignment/matrix_align_oracle.tsv` (1 986 rows: 331 random protein / DNA / IUPAC-DNA pairs x global|local|semi x linear|affine, integer and fractional gap scores; Biopython 1.88, EMBOSS needle `-endweight Y` / water 6.6.0, parasail 1.3.4; generator `scripts/alignment/generate_matrix_align_oracle.py`).

| Case | Test | Locks |
|------|------|-------|
| Score = Biopython (all rows) and parasail (integer-gap rows) | `Oracle_Score_EqualsBiopythonAndParasail_AllRows` | 1 986 / 616 rows |
| EMBOSS needle/water agreement (>= 98 %; needle never above the optimum) | `Oracle_Emboss_NeedleAndWater_AgreeExceptKnownEndGapQuirk` | 1 103 of 1 121 equal; embaln.c end-gap quirk documented |
| Returned alignment re-scores to its score, reproduces inputs, coordinate contracts | `Oracle_Alignment_RescoresToScore_AndReproducesInputs` | independent scorer |
| HBA/HBB with EMBOSS defaults | `HbaVsHbb_*` | needle 292.5, water 293.5 |
| Gap-open convention (EMBOSS gapopen = gap of length 1) | `SingleGap_ConventionIsEmbossGapopenForLengthOne` | 10 / 19.5 / 9.5 |
| One DP: integer matrix + integer gaps = scalar aligners row for row | `IntegerMatrix_IntegerGaps_IsIdenticalToScalarAligners_RowForRow` | 6 x 400 random pairs |
| Overloads, case-insensitivity, validation, empty inputs, GapOpen 0 = linear | remaining tests | string = ProteinSequence = DnaSequence |

## 2026-11 audit round 1 addition (A1-5, B09 F45): general gap-cost functions (Waterman-Smith-Beyer)

**Scope.** `GlobalAlignGeneralGap` / `LocalAlignGeneralGap` / `SemiGlobalAlignGeneralGap(…, EndGaps)` with `Func<int,double> gapScore` (any w(k): concave, convex, capped, step, two-slope, non-monotone), scalar `ScoringMatrix` or `SubstitutionMatrix` scoring, `string`/`DnaSequence`/`ProteinSequence` overloads. Model: every maximal gap run scores w(k) once (Biopython's callable `gap_score`; never split). O(mn(m+n)); linear/affine functions delegate to the O(mn) programs. Tests: `SequenceAligner_GeneralGap_Tests` (10) + MCP `GeneralGapAlignTests`.

| ID | Test | Data | Expected | Evidence |
|----|------|------|----------|----------|
| GG-1 | `Oracle_Score_GeneralProgramAndPublicApi_EqualBiopythonAndReference` | 1 440 rows: 960 global/local (8 families x 60 pairs x 2 modes; scalar / BLOSUM62 / EDNAFULL) + 480 ends-free | general program and public API = Biopython 1.88 `PairwiseAligner(gap_score=callable)` (global/local) = independent three-state reference (all rows), 1e-8 | general_gap_oracle.tsv |
| GG-2 | `Oracle_Alignment_RescoresToScore_AndCoordinatesAreConsistent` | all rows, API and general program | rows reproduce inputs; independent run-charging re-scoring = score; local interval / global coordinates / tags | independent rescorer |
| GG-3 | `AffineAndLinearFunctions_UseTheAffineAndLinearAligners_RowForRow` | 150 pairs x 4 schemes, global/local/semi + 4 end policies | rows = `GlobalAlignAffine` / `LocalAlignAffine` / `SemiGlobalAlignAffine(EndGaps)` / `GlobalAlign` (GapOpen 0); general program same score | fast path |
| GG-4 | `FractionalAffineFunction_WithBlosum62_EqualsTheMatrixAffineAligner_AndTheGeneralProgram` | w(k) = -(9.5 + 0.5 k), 60 protein pairs | = `GlobalAlignAffine(MatrixScoring.Emboss)` incl. exact fractional score | EMBOSS gap numbers |
| GG-5 | `BruteForce_EveryAlignment_MaximumEqualsTheGeneralProgram` | 60 pairs (<= 5), 6 gap functions (incl. non-monotone), 16 policies; local on 15 pairs vs all substring-pair alignments | max = program | exhaustive definition |
| GG-6 | `ConvexAndConcaveGaps_BiopythonValues` | AAAATTTT/AAAACCCCTTTT (2,−3): w = -k^2 → 5 (a run of 4 would give 0); w = -(2 + 3 ln k) → 9.841116916640328, `AAAA----TTTT`/`AAAACCCCTTTT` (unique) | exact | Biopython 1.88 callable gap_score |
| GG-7 | `GapFunction_IsEvaluatedOncePerLength_AndMustBeValid`, `EmptyInputs_FollowTheAlignerContract`, `DnaAndProteinOverloads_EqualTheStringOverloads` | NaN/+inf rejected, -inf forbids a length, null guards, residue outside the matrix alphabet, empty contract | as stated | API contract |

## 2026-11 audit round 1 addition (A1-6, B09 F46): exact linear-space alignment

**Scope.** `SequenceAligner.LinearSpaceCellThreshold` (25 M cells: automatic switch of every pairwise aligner) and the explicit `GlobalAlignLinearSpace`, `GlobalAlignAffineLinearSpace`, `LocalAlignLinearSpace`, `LocalAlignAffineLinearSpace`, `SemiGlobalAlignLinearSpace(…, EndGaps)`, `SemiGlobalAlignAffineLinearSpace(…, EndGaps)`; algorithm document `docs/algorithms/Alignment/Linear_Space_Alignment.md`. Contract: score, rows, coordinates, type and tie behaviour IDENTICAL to the full-matrix aligner; memory `O(m + n)`.

Tests: `SequenceAligner_LinearSpace_Tests` (20) and `Seqeron.Mcp.Alignment.Tests.LinearSpaceAlignTests` (3). Oracle `TestData/Alignment/linear_space_oracle.tsv` (`scripts/alignment/generate_linear_space_oracle.py`): 1 280 Biopython 1.88 / parasail 1.3.4 rows (global, local, fitting and nine other end-gap policies, linear + affine, ties, long indels), 160 EMBOSS 6.6.0 `stretcher` rows (Myers-Miller) and 160 `matcher` rows (Huang-Miller SIM).

| ID | Test | Data | Expected | Evidence |
|----|------|------|----------|----------|
| LS-1 | `LinearSpace_EqualsFullMatrix_RowForRow_AllFamiliesAndPolicies_ForcedRecursion` | 45 pairs (alphabets A / AC / ACGT, related + unrelated, long indels), linear + affine, local + 16 `EndGaps` policies, base cells 1-70 | `AlignmentResult` equal (score, rows, coordinates, type) | full-matrix core (`forceFullMatrix`) |
| LS-2 | `LinearSpace_EqualsFullMatrix_ManyRandomPairs_Global_Local_Fitting_Affine_Linear` | 700 pairs x {global, local, fitting, all-free} = 2 800 | same | full-matrix core |
| LS-3 | `LinearSpace_MatrixScoring_EqualsFullMatrix_Protein_Dna_FractionalGaps` | 60 BLOSUM62 / EDNAFULL pairs x 7 families, gaps (-9.5,-0.5) … (0,0) | same as the public full-matrix `MatrixScoring` aligners | F41 aligners |
| LS-4 | `LinearSpace_Overloads_DnaProteinString_Agree` | DnaSequence / ProteinSequence / string overloads | equal | API |
| LS-5 | `LinearSpace_KilobasePairs_EqualFullMatrix_AtDefaultBase` | 2.6-3 kb related pairs with two planted indels, affine/linear, global/fitting/all-free/overlap/local, default base and 50 000-cell base | rows equal | full-matrix core |
| LS-6 | `AutomaticSwitch_AboveThreshold_ReturnsTheFullMatrixAlignment`, `AutomaticSwitch_MatrixScoring_ProteinAbove25MCells` | 5 300 x ~5 200 (27.6 M cells) all six aligners; protein BLOSUM62 5 200 x 5 200 | public aligner = full-matrix core row for row | threshold contract |
| LS-7 | `LinearSpace_Memory_IsLinear_NotQuadratic`, `AutomaticSwitch_Cancellation_IsHonoured` | 5 600 x 5 600 affine (matrices 376 MB); 8 000 x 8 000 with a 30 ms token | allocation < 1/8 of the matrices; `OperationCanceledException` | O(m+n) memory |
| LS-8 | `GlobalAffine_LongGap_KnownValue_ThroughRecursion`, `LocalAffine_ReferenceValue_Biopython` | ACGTTTTACGT/ACGTACGT (BLAST DNA) = 5; F2 local value 64 | exact | Biopython / parasail (F1, F2) |
| LS-9 | `Oracle_Biopython_Parasail_ScoresAndCoOptimalRows_ThroughRecursion` | 1 280 rows, base cells 1-60 | score = Biopython (and parasail where a variant exists), independent re-scoring = score, alignment ∈ Biopython co-optimal set (≤ 40), = full-matrix rows | Biopython 1.88, parasail 1.3.4 |
| LS-10 | `Oracle_EmbossStretcher_MyersMiller_GlobalScoreAndRows`, `Oracle_EmbossMatcher_HuangMillerSim_LocalScoreAndRows` | 160 + 160 rows (EDNAFULL / EBLOSUM62, gaps 16/4, 10/0.5, 12/2, 8/8, 20/1) | score equal; EMBOSS rows re-score to the score; rows (+ matcher coordinates) equal where Biopython finds a unique optimum (49 / 93 rows) | EMBOSS 6.6.0 |
| LS-11 | `EmptyInputs_AndTinySequences_FollowTheFullMatrixContract`, `Validation_NullInputs_And_PositiveLocalGapScores`, `LowerCaseInput_IsUpperCased_LikeEveryAligner` | empty / 1-char inputs x all policies; null; positive local gap scores; matrix alphabet | = full-matrix contract; `ArgumentNullException` / `ArgumentException` | API contract |
| LS-12 | MCP `GlobalAlignAffine_AboveThreshold_…`, `LocalAlignAffine_AboveThreshold_…`, `GeneralGapAlign_NonAffineGapFunction_AboveCellLimit_…` | 5 300 x ~5 300 tools; 2 100 x 2 100 log gap table | = explicit linear-space result; clear `ArgumentException` above `GeneralGapMaxCells` | tool contract |

---

## 2026-11 audit round 1 (A1-7b, B09 F59): co-optimal global / ends-free alignments

`GlobalAlignCoOptimal` / `GlobalAlignAffineCoOptimal` / `SemiGlobalAlignCoOptimal` / `SemiGlobalAlignAffineCoOptimal` (`ScoringMatrix` string + `DnaSequence`, `MatrixScoring` string + `ProteinSequence` + `DnaSequence`; `maxAlignments` cap; MCP `global_align_co_optimal`). Test class `SequenceAligner_GlobalCoOptimal_Tests`, oracle `TestData/Alignment/global_cooptimal_oracle.tsv` (1 920 Biopython 1.88 `PairwiseAligner(mode='global').align` rows, all 16 end-gap policies, linear/affine, scalar + BLOSUM62 fractional gaps; `scripts/alignment/generate_global_cooptimal_oracle.py`).

| ID | Test | Input | Expected | Source |
|----|------|-------|----------|--------|
| GC-1 | `CoOptimalSet_EqualsBiopythonSet_OnEveryOracleRow` | 1 859 rows with <= 150 alignments (935 with several) | set of (row1, row2) equal, no duplicates | Biopython 1.88 |
| GC-2 | `EveryAlignment_RescoresToTheBiopythonOptimum_WithConsumedIntervals` | every alignment of every row | independent re-score = Biopython optimum, coordinates = consumed intervals of the rows, rows spell the inputs, type Global / SemiGlobal | definition |
| GC-3 | `FirstAlignment_IsTheSingleAligners_OnEveryOracleRow` | maxAlignments = 1 | = `GlobalAlign` / `Affine` / `SemiGlobalAlign*` (zero gap costs with a free flank: all fields but the coordinates) | canonical single aligner |
| GC-4 | `Cap_OnRowsWithMoreThanTheOracleLimit_ReturnsDistinctOptimalPrefix` | 61 rows with > 150 alignments, cap 50 | 50 distinct alignments, all optimal | definition |
| GC-5 | `CoOptimalSet_EqualsExhaustiveEnumeration_AllSixteenPolicies_LinearAndAffine` | pairs <= 5 residues, every column sequence enumerated | set equal, 16 policies, linear + affine | definition (exhaustive) |
| GC-6 | `GlobalExample_GaacVsGaaac_ThreeGapPlacements`, `BiopythonTutorialExample_GaactVsGat_TwoAlignments` | hand / Biopython tutorial | 3 / 2 alignments | Biopython tutorial |
| GC-7 | `StringDnaProteinAndMatrixOverloads_ReturnTheSameAlignments`, `EndGapsNone_IsTheGlobalProblem_AndDefaultPolicyOfTheGlobalNames`, `EmptySide_HasTheSingleAllGapAlignment`, `ZeroGapCosts_FreeFlankIsListedOnce` | overload equivalence, `None` = global, empty side, zero costs | equal / single alignment | contract |
| GC-8 | `Validation_NullCapEndGapsAndMemoryGuard`, `ManyTies_CapStopsTheEnumeration_Quickly` | null, cap < 1, bad `EndGaps`, > 25 M cells, 60 x 60 all-ties | exceptions; cap honoured quickly | contract |
| GC-9 (MCP) | `GlobalCoOptimalAndExactMsaTests.GlobalCoOptimal_*` | gap placements, affine, `endGaps`, BLOSUM62 fractional, validation | delegates to the core; first = `global_align` | library call |

---

## 2026-11 audit round 2 (A2-37, A2-43; B09 F80, F85): one scoring contract and `null` for the scoring argument

Scalar `ScoringMatrix` and `MatrixScoring` share ONE contract (`PairScoring`, checked by every aligner): `GapOpen` and `GapExtend` are penalties (<= 0, greater than `int.MinValue`; a positive gap score is an `ArgumentException` "must be a finite score <= 0", previously silently accepted by the scalar path: `GlobalAlign("ACGT","AC", new ScoringMatrix(1,-1,5,5))` = 30), and the score of any path, `max(|Match|, |Mismatch|, |GapOpen| + |GapExtend|) * (m + n + 2)`, must stay below 2^28 (`ArgumentException` "Score range too large for the integer dynamic program"; previously a silently wrapped `Score`: `GlobalAlign("A"x50, "A"x50, new ScoringMatrix(1e9,-1,-2,-1))` = 1 999 999 904 for the true 5e10). The multiple-alignment sum-of-pairs totals accumulate in 64 bits and are refused when they do not fit the reported 32-bit total.

| ID | Test (`SequenceAligner_ScoringContract_Tests`) | Input | Expected | Source |
|----|------|-------|----------|--------|
| SC-1 | `ScalarOverflow_IsRejectedByEveryScalarAligner_NotWrapped` | the 5e10 repro (match 1e9) through 14 aligners (global / local / semi-global, linear / affine, linear space, top-N, co-optimal) | `ArgumentException` "Score range too large" | exact integer arithmetic (50 x 1e9) |
| SC-2 | `ScalarLargeScores_BelowTheBound_AreExact` | match 1e6, 50 identities | 50 000 000 for global / local / affine / semi-global | closed form n x match |
| SC-3 | `ScalarRangeBound_IsTheSameAsTheMatrixPathBound` | match 1e6: m = n = 133 accepted, 134 rejected | the 2^28 bound | definition |
| SC-4 | `IntMinValueGaps_AreRejected_WithArgumentException_NotIndexOutOfRange` | gap = int.MinValue (negation overflow) | `ArgumentException` | contract |
| SC-5 | `PositiveGapScore_IsRejected_ByScalarAndMatrixModelAlike` | scalar (5,5), (0,1), (1,0), (-1,2) and `MatrixScoring` (1,-1) | same message; zero gaps stay valid | contract |
| SC-6 | `IntMinValueGaps_ThroughMultipleAlign_SurfaceAsArgumentException_ForEveryK` | k = 2, 3, 5, affine | `ArgumentException` (was `AggregateException` for k = 3, 5) | A2-42 |
| SC-7 | `SumOfPairs_Accumulates64Bit_AndRefusesATotalBeyondInt` | 3 rows x 2 columns, match 2e9 (total 1.2e10) | `ArgumentException`; `MatrixSumOfPairsScore` = 24.0 (BLOSUM62 A/A = 4 x 6) | exact arithmetic, BLOSUM62 |
| OV-1 | `SequenceAligner_InputContract_Tests.NullScoring_Compiles_AndSelectsTheDefault_OnEveryPairwiseFamily` | `GlobalAlign(a, b, null)`, `LocalAlignTopN(a, b, null, 3, 1)` ... 17 calls (CS0121 before `[OverloadResolutionPriority(1)]`) | same result as the call without the argument | compile-level + value |
| OV-2 | `NullScoring_TheTypedMatrixOverloadsStillBind`, `EveryScoringMatrixOverloadWithAMatrixScoringTwin_CarriesTheResolutionPriority` | reflection over all public overloads | every `ScoringMatrix?` overload with a `MatrixScoring` twin carries priority >= 1 | contract |
| MCP-1 | `EmptyInputContractTests` (Mcp.Alignment.Tests) | `global_align("", "ACGT")` ... | the library's all-gap alignment / empty lists; null still rejected | library contract (F39) |

## Audit A2-29 (B09 F123): ksw2-style dual affine gaps

`SequenceAligner_DualAffine_Tests`: 1 600 oracle cases (global / local x full / banded; Biopython callable gap score + independent WSB), re-scoring and coordinates of every alignment, single-regime equality with `GlobalAlignAffine` / `LocalAlignAffine`, long-gap regime, contract; MCP `DualAffineAlign`.

## Audit A2-28c / A2-29b (B09 F124 / F125): needleall, ksw2 extension alignment
`SequenceAligner_EmbossPrograms_Tests`: 260 needleall runs (whole output files, error files, `-endweight` / `-endopen` / `-endextend`, `-minscore`, `-nobrief`). `Ksw2Aligner_Tests` (`ksw2_exact_oracle.tsv`, 1 320 cases, 330 per variant extz / extz2 / extd / extd2): max, mqe, mte and score equal an independent Gotoh program over the band, every CIGAR is re-scored along its own path to score / max / mqe (corner, best cell, query end with the end bonus), `reach_end`; identical sequences, N penalty, Z-drop stop, end bonus. The Z-drop and approximate-mode cases wait for the real program (handoff A2-29b-ksw2).

## Audit A4-3 (B09 F153): EMBOSS reader - gap characters and other non-letters

`EmbossGapCharacters_Oracle_Tests` (`emboss_gap_oracle.tsv`, `generate_emboss_gap_oracle.py`): needle (240 rows, e.g. `AC-GTTAGC.TT*` vs `ACGTTAGCTT` = score 58, `Length: 11`; 43 before), stretcher (240) and needleall (240) - the installed EMBOSS 6.6.0 `needle, stretcher, needleall` on inputs that contain `- . ~ ? # + _`, digits, blanks, punctuation, `*` and non-ASCII letters equal `NeedleAlign`, `StretcherAlign`, `NeedleAll` on every row, and the result of the cleaned input equals the result of the input. The reader rule (one shared `EmbossSequenceReader.Clean` / `CleanGapped`, Infrastructure): only `A-Z a-z *` survive in a plain sequence (`ajseqread.c seqAppend` + `ajseqtype.c ajSeqTypeCheckIn` + `ajstr.c ajStrRemoveGap`); coordinates and lengths count the cleaned sequence. Also `Clean_KeepsLettersAndStar_AndDropsEveryOtherCharacter` (all 256 code points) and `IsNucleotide_IsTheAjSeqIsNucRule_OfEveryPort`.

## Audit R14-1 (B09 F288): one End convention for every `AlignmentResult` producer

`AlignmentResult_EndConvention_Tests` (Genomics.Tests): `EndPosition1` / `EndPosition2` are 0-based INCLUSIVE (length - 1 for a global result, the last consumed index of a local / fitting result, `Start - 1` for an empty consumed interval, -1 / -1 for an empty local result) in every producer: `GlobalAlign`, `GlobalAlignAffine`, `GlobalAlignLinearSpace`, `GlobalAlignGeneralGap`, `GlobalAlignCoOptimal`, `LocalAlign` (+ affine, linear space, general gap, top-N, SIM top-N), `SemiGlobalAlign` (+ affine, ends-free), `NeedleAlign`, `WaterAlign`, `StretcherAlign`, `GlobalAlignDualAffine`, `LocalAlignDualAffine` over 7 fixed + 60 seeded pairs; the pair "TTACGTTT" / "ACGT" gives end 7,3 (global, needle, stretcher, dual) and 2..5 / 0..3 (local, water, dual local). `EmbossAlignedSequence.Flanks` (offset / offend from the rows: a full row has none, a local row `Start` and the rest, an empty local result the whole sequence as the end flank) and the unchanged offsets of needle / water / stretcher. MCP `AlignmentResultFormatOffsetsTests`: the markx10 / pair / sam of `local_align` / `semi_global_align` / `global_align` / `needle_align` / `stretcher_align` / `dual_affine_align` equal the numbers printed by the installed EMBOSS 6.6.0 `water` (">>>a, 8 nt vs b, 4 nt", al_start 3 / al_stop 6, "a 3 ACGT 6", sam "4M"); the fitting rows run 1 .. 8.

## Audit R14-15 (B09 F294): explicit cancellation token on every pairwise aligner

`SequenceAligner_CancellationOverloads_Tests`: reflection - every public static pairwise aligner of `SequenceAligner` (AlignmentResult / IReadOnlyList<AlignmentResult>, first parameters `sequence1`, `sequence2`) has an overload with the same parameters and a trailing `CancellationToken` (generated by `scripts/alignment/generate_cancellable_overloads.py`, 163 overloads); 25 families (GlobalAlign, GlobalAlignAffine, LocalAlign, LocalAlignAffine, SemiGlobalAlign with ScoringMatrix / EndGaps / GapScoreProfile, LinearSpace, GeneralGap, CoOptimal, TopN (Waterman-Eggert and SIM), TopNLinearSpace, Needle, Water, Stretcher, dual affine): the result with a live token equals the original, a pre-cancelled token throws `OperationCanceledException`; mid-run cancellation of LocalAlignAffine / GlobalAlignAffine / LocalAlignTopN / NeedleAlign / WaterAlign stops at the named poll.
