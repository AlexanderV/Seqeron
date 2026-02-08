# Test Specification: ALIGN-LOCAL-001

**Test Unit ID:** ALIGN-LOCAL-001
**Area:** Alignment
**Algorithm:** Local Alignment (Smith–Waterman)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-01

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Smith–Waterman algorithm | https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | 2026-02-01 |
| Wikipedia: Sequence alignment | https://en.wikipedia.org/wiki/Sequence_alignment | 2026-02-01 |

### 1.2 Key Evidence Points

1. Smith–Waterman performs local sequence alignment, finding the highest-scoring pair of subsequences between two sequences. (Smith–Waterman algorithm)
2. The key difference from Needleman–Wunsch is the **zero floor**: negative scores are set to 0, allowing alignment to restart at any position. (Smith–Waterman algorithm)
3. Initialization: first row and first column are set to 0 (no end-gap penalty). (Smith–Waterman algorithm)
4. Traceback: begins at the cell with the **highest score** in the matrix and ends when a cell with score 0 is encountered. (Smith–Waterman algorithm)
5. Recurrence with linear gap penalty: $H_{i,j} = \max(0, H_{i-1,j-1} + s(a_i,b_j), H_{i-1,j} - W_1, H_{i,j-1} - W_1)$ (Smith–Waterman algorithm)
6. Example from Wikipedia: sequences `TGTTACGG` and `GGTTGACTA` with match +3, mismatch -3, linear gap penalty $W_1 = 2$ yields alignment `GTT-AC` / `GTTGAC`. (Smith–Waterman algorithm)

### 1.3 Documented Corner Cases

| Case | Evidence Source | Notes |
|------|-----------------|-------|
| No similarity between sequences | Smith–Waterman algorithm (zero floor) | When all comparisons yield negative scores, result should be score 0 or empty alignment |
| Identical sequences | ASSUMPTION | Should align fully with maximum score |
| Empty sequence input | ASSUMPTION | Implementation returns `AlignmentResult.Empty` |

### 1.4 Known Failure Modes / Pitfalls (from sources)

1. Multiple optimal local alignments can exist when several cells share the same maximum score. (Smith–Waterman algorithm)
2. Alignment quality depends on the scoring system and gap penalties. (Smith–Waterman algorithm)
3. Linear gap penalty may not model biological indels accurately; affine gaps (Gotoh) are often preferred. (Smith–Waterman algorithm)

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` | SequenceAligner | **Canonical** | Smith–Waterman local alignment |
| `LocalAlign(string, string, ScoringMatrix?)` | SequenceAligner | Delegate | String wrapper over core implementation |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Score ≥ 0 (zero floor property) | Yes | Smith–Waterman algorithm (negative scores → 0) |
| INV-2 | Aligned region is a contiguous subsequence of the original | Yes | Smith–Waterman traceback from max to 0 |
| INV-3 | AlignmentType is Local | Yes | Implementation contract |
| INV-4 | Removing gaps from aligned sequences yields substrings of originals | Yes | Smith–Waterman traceback rules |
| INV-5 | StartPosition and EndPosition are within sequence bounds | Yes | Smith–Waterman traceback from matrix positions |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Wikipedia example: local alignment finds correct subsequences | seq1 = `TGTTACGG`, seq2 = `GGTTGACTA`, scoring = match +3, mismatch −3, gap −2 | Aligned region contains best local match, INV-1..INV-5 hold | Smith–Waterman algorithm example |
| M2 | Score is non-negative (zero floor) | Any input | result.Score ≥ 0 | Smith–Waterman algorithm |
| M3 | AlignmentType is Local | Any valid input | result.AlignmentType == Local | Implementation contract |
| M4 | Aligned positions are valid | Any valid input | 0 ≤ StartPosition ≤ EndPosition ≤ sequence length | Smith–Waterman traceback |
| M5 | Removing gaps from aligned sequences yields originals' substrings | Wikipedia example | Gaps removed = substring of original | Smith–Waterman traceback rules |
| M6 | String overload returns same result as DnaSequence overload | Same sequences | Same alignment | ASSUMPTION (wrapper parity) |
| M7 | Empty string input returns AlignmentResult.Empty | seq1 = "", seq2 = "ACGT" | Empty result | ASSUMPTION (implementation) |
| M8 | Null DnaSequence throws ArgumentNullException | seq1 = null, seq2 = valid | ArgumentNullException | ASSUMPTION (implementation) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| S1 | Identical sequences produce full-length alignment | seq1 = seq2 = "ACGTACGT" | Full match, high score | ASSUMPTION |
| S2 | Completely dissimilar sequences produce score 0 or minimal | seq1 = "AAAA", seq2 = "TTTT" | Score ≥ 0, possibly empty alignment | Zero floor property |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| C1 | Multiple optimal alignment detection | Sequences with tied max scores | Any optimal alignment is acceptable | Evidence: multiple optimal alignments can exist |

---

## 5. Audit of Existing Tests

### 5.1 Test Discovery Summary

Existing tests found in [Seqeron.Genomics.Tests/SequenceAlignerTests.cs](../tests/Seqeron/Seqeron.Genomics.Tests/SequenceAlignerTests.cs) `#region Local Alignment Tests`:

| Test Method | Lines | Status |
|-------------|-------|--------|
| `LocalAlign_FindsBestSubsequence()` | 338-346 | Weak: checks AlignmentType and Contains("TGC") |
| `LocalAlign_NoSimilarity_ReturnsEmptyOrLow()` | 348-356 | **Weak**: assertion `result.Score >= 0` is always true |
| `LocalAlign_IdenticalSequences_FullMatch()` | 358-367 | Adequate: checks aligned sequence |
| `LocalAlign_ReturnsPositions()` | 369-378 | Adequate: checks position validity |
| `LocalAlign_StringOverload_Works()` | 380-389 | Smoke test for wrapper |
| `LocalAlign_NullSequence_ThrowsException()` | ~null handling region | Edge case: null throws |

### 5.2 Coverage Classification (Pre-Consolidation)

| Area | Status | Notes |
|------|--------|-------|
| Zero floor invariant (INV-1) | ✅ Covered | Strengthened with proper assertions |
| Wikipedia example dataset | ✅ Covered | Dataset-based verification added |
| AlignmentType is Local (INV-3) | ✅ Covered | `FindsBestSubsequence` tests this |
| Position validity (INV-5) | ✅ Covered | `ReturnsPositions` tests bounds |
| Empty input handling | ✅ Covered | Added |
| Null argument handling | ✅ Covered | Present in null handling region |
| Gaps removal invariant (INV-4) | ✅ Covered | Added |

### 5.3 Consolidation Plan

- **Canonical file:** new [Seqeron.Genomics.Tests/SequenceAligner_LocalAlign_Tests.cs](../tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_LocalAlign_Tests.cs)
  - Evidence-based invariant tests using Wikipedia example
  - All MUST tests (M1-M8)
  - SHOULD tests (S1-S2)
- **Remove** local alignment tests from shared [Seqeron.Genomics.Tests/SequenceAlignerTests.cs](../tests/Seqeron/Seqeron.Genomics.Tests/SequenceAlignerTests.cs)
- **Remove** weak `LocalAlign_NoSimilarity_ReturnsEmptyOrLow` test (meaningless assertion)

### 5.4 Final State After Consolidation

| File | Tests | Role |
|------|-------|------|
| [Seqeron.Genomics.Tests/SequenceAligner_LocalAlign_Tests.cs](../tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_LocalAlign_Tests.cs) | 10 | Canonical + SHOULD |

---

## 6. Open Questions / Decisions

None.
