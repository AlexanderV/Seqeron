# Test Specification: ALIGN-SEMI-001

**Test Unit ID:** ALIGN-SEMI-001  
**Area:** Alignment  
**Algorithm:** Semi-Global Alignment (Ends-Free / Glocal)  
**Status:** ⏳ In Progress  
**Owner:** Algorithm QA Architect  
**Last Updated:** 2026-02-01  

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Sequence alignment | https://en.wikipedia.org/wiki/Sequence_alignment#Global_and_local_alignments | 2026-02-01 |
| Wikipedia: Needleman–Wunsch algorithm | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | 2026-02-01 |
| Brudno et al. (2003) | doi:10.1093/bioinformatics/btg1005 | 2026-02-01 |

### 1.2 Key Evidence Points

1. Semi-global alignment (also "glocal") searches for the best partial alignment, allowing one or both ends to be unaligned without penalty. (Wikipedia: Sequence alignment)
2. Primary use case: aligning a short query to a longer reference where the query should be globally aligned but the reference may have free end gaps. (Wikipedia: Sequence alignment)
3. Useful for overlap alignment in sequence assembly. (Wikipedia: Sequence alignment)
4. Algorithm is a modification of Needleman–Wunsch with modified initialization and traceback. (ASSUMPTION from algorithmic theory)
5. Complexity: O(n × m) time and space. (Wikipedia: Needleman–Wunsch algorithm)

### 1.3 Documented Corner Cases

| Case | Evidence Source | Expected Behavior |
|------|-----------------|-------------------|
| Short query in long reference | Wikipedia (Sequence alignment) | Query fully aligned, reference has leading/trailing gaps |
| Query equals reference | ASSUMPTION | Full alignment like global |
| Query overlaps reference end | Wikipedia (Sequence alignment) | Overlap detected without trailing gap penalty |
| Null sequence | ASSUMPTION | ArgumentNullException |

### 1.4 Known Failure Modes / Pitfalls

1. Multiple equally optimal alignments may exist (like all DP alignment algorithms).
2. Implementation variant (query-in-reference) may not be suitable for all semi-global use cases.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` | SequenceAligner | **Canonical** | Query-in-reference semi-global alignment |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | AlignmentType is SemiGlobal | Yes | Implementation contract |
| INV-2 | Aligned sequences have equal length | Yes | Alignment definition |
| INV-3 | Removing gaps from aligned seq1 yields original seq1 | Yes | Query fully aligned |
| INV-4 | Removing gaps from aligned seq2 yields substring of original seq2 | Yes | Reference has free end gaps |
| INV-5 | Score reflects matches/mismatches/gaps correctly | Yes | Scoring consistency |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Short query embedded in long reference | query = "ATGC", ref = "AAAATGCAAA" | Query fully aligned, AlignmentType = SemiGlobal | Wikipedia (Sequence alignment) |
| M2 | AlignmentType is SemiGlobal | Any valid input | result.AlignmentType == SemiGlobal | Implementation contract |
| M3 | Aligned sequences have equal length | Any valid input | len(aligned1) == len(aligned2) | Alignment definition |
| M4 | Query fully represented | query = "ATGC", ref = "ATGCAAAA" | RemoveGaps(aligned1) == "ATGC" | INV-3 (query fully aligned) |
| M5 | Reference is substring after gap removal | query = "ATGC", ref = "AAAATGCAAA" | RemoveGaps(aligned2) is substring of ref | INV-4 |
| M6 | Score is non-negative for matching sequences | query = "ATGC", ref = "ATGCAAAA" | Score >= 0 | ASSUMPTION |
| M7 | Null sequence throws ArgumentNullException | seq1 = null | ArgumentNullException | ASSUMPTION (implementation) |
| M8 | Null sequence (seq2) throws ArgumentNullException | seq2 = null | ArgumentNullException | ASSUMPTION (implementation) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| S1 | Identical sequences produce full alignment | query = ref = "ATGCATGC" | Full match, score = len × matchScore | ASSUMPTION |
| S2 | Query at start of reference | query = "ATG", ref = "ATGCCCCC" | Trailing gaps in reference, full query match | Semi-global overlap |
| S3 | Query at end of reference | query = "CCC", ref = "ATGCCC" | Leading gaps in reference portion | Semi-global overlap |
| S4 | Custom scoring matrix | Custom scores | Score computed correctly | Scoring contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| C1 | Long query, short reference | query longer than ref | Alignment computed (query may have gaps) | Edge case |
| C2 | Single base sequences | query = "A", ref = "A" | Match score = 1 | Minimal input |

---

## 5. Audit of Existing Tests

### 5.1 Test Discovery Summary

| File | Tests Found | Coverage |
|------|-------------|----------|
| SequenceAlignerTests.cs | 3 tests (semi-global region) | Weak: basic smoke tests only |

### 5.2 Existing Test Classification

| Test Method | Classification | Notes |
|-------------|----------------|-------|
| SemiGlobalAlign_ShortInLong_FindsMatch | **Weak** | Only checks AlignmentType, no invariants |
| SemiGlobalAlign_FreeEndGaps_NoGapPenalty | **Weak** | Score >= 0 check, basic seq1 reconstruction |
| SemiGlobalAlign_NullSequence_ThrowsException | Covered | Null handling |

### 5.3 Consolidation Plan

1. **Move existing tests**: Extract semi-global tests from SequenceAlignerTests.cs to new canonical file
2. **Strengthen tests**: Add invariant validation (INV-1 through INV-5)
3. **Add missing Must tests**: M3, M4, M5, M6, M8
4. **Add Should tests**: S1, S2, S3, S4
5. **Remove weak tests**: Replace with invariant-focused versions

---

## 6. Open Questions / Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| Which semi-global variant is implemented? | Query-in-reference | Implementation analysis: first row = 0, first col = gap penalties |
| Should we test other variants? | No | Out of scope for this test unit |

---

## 7. Validation Checklist

- [ ] All MUST tests implemented and passing
- [ ] Tests are deterministic
- [ ] Invariants verified with Assert.Multiple where appropriate
- [ ] No duplicate tests across files
- [ ] Evidence documented for all test rationale
