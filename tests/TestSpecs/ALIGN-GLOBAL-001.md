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
| `GlobalAlign(string, string, scoring)` | SequenceAligner | Delegate | Wrapper: normalizes to uppercase, returns `AlignmentResult.Empty` for empty input |
| `GlobalAlign(string, string, scoring, CancellationToken, IProgress?)` | SequenceAligner | Delegate | Cancellation wrapper (tested as smoke elsewhere) |

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

### 4.2 API Contract Tests

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| A1 | Empty string input returns AlignmentResult.Empty | seq1="", seq2=GATTACA | Empty result | API contract — empty input guard |
| A2 | Null DnaSequence throws ArgumentNullException | seq1=null, seq2=GATTACA | ArgumentNullException | API contract — null guard |

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
| `EmptyStringInput_ReturnsEmptyResult` | A1 |
| `NullDnaSequence_ThrowsArgumentNullException` | A2 |

---

## 6. Deviations and Assumptions

**None.** All tests and implementation are grounded in the authoritative sources listed in §1.1.

- The implementation uses the standard Needleman–Wunsch linear gap penalty model exactly as described in the Wikipedia pseudocode.
- `ScoringMatrix.GapExtend` serves as the linear gap penalty `d`. `ScoringMatrix.GapOpen` is not used by `GlobalAlign` (it exists in the record for other alignment types).
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
