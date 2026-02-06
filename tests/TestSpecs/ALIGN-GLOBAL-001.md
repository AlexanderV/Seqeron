# Test Specification: ALIGN-GLOBAL-001

**Test Unit ID:** ALIGN-GLOBAL-001
**Area:** Alignment
**Algorithm:** Global Alignment (Needleman–Wunsch)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-01

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Needleman–Wunsch algorithm | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | 2026-02-01 |
| Wikipedia: Sequence alignment (global alignment) | https://en.wikipedia.org/wiki/Sequence_alignment | 2026-02-01 |

### 1.2 Key Evidence Points

1. Global alignment (Needleman–Wunsch) uses dynamic programming to compute the optimal end-to-end alignment score for two sequences given a scoring system for matches, mismatches, and gaps (indels). (Needleman–Wunsch algorithm)
2. The scoring function is the sum over aligned positions: match, mismatch, or indel (gap), with initialization of the first row/column based on gap penalties. (Needleman–Wunsch algorithm)
3. Traceback rules align a character to a gap when the optimal path comes from a horizontal or vertical move. (Needleman–Wunsch algorithm)
4. Global alignment attempts to align the entire length of the sequences; gaps are inserted so aligned sequences are the same length. (Sequence alignment, global vs local)
5. Example sequences shown in the Needleman–Wunsch page: GCATGCG and GATTACA with a simple scoring scheme (match = +1, mismatch/indel = −1). (Needleman–Wunsch algorithm)

### 1.3 Documented Corner Cases

No authoritative sources explicitly specify empty-sequence handling for Needleman–Wunsch in this repository’s API. Empty-input behavior is therefore treated as **ASSUMPTION** based on current implementation.

### 1.4 Known Failure Modes / Pitfalls (from sources)

1. Alignment quality depends on the scoring system and gap penalties; different choices can yield different optimal alignments. (Needleman–Wunsch algorithm)
2. Multiple optimal alignments can exist when several traceback paths share the same score. (Needleman–Wunsch algorithm)

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GlobalAlign(seq1, seq2, scoring)` | SequenceAligner | **Canonical** | Needleman–Wunsch global alignment |
| `GlobalAlign(string, string, scoring)` | SequenceAligner | Delegate | Wrapper over core implementation |
| `GlobalAlign(string, string, scoring, CancellationToken, IProgress<double>?)` | SequenceAligner | Delegate | Cancellation wrapper (tested as smoke elsewhere) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Aligned sequences have equal length | Yes | Sequence alignment (global alignment) |
| INV-2 | Removing gaps from aligned sequences yields the original sequences | Yes | Needleman–Wunsch traceback rules |
| INV-3 | Alignment score equals sum of match/mismatch/indel scores across aligned positions | Yes | Needleman–Wunsch scoring definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Wikipedia example: global alignment reconstructs inputs and score matches scoring | seq1 = GCATGCG, seq2 = GATTACA, scoring = match +1, mismatch −1, gap −1 | Invariants INV-1..INV-3 hold | Needleman–Wunsch algorithm example and scoring rules |
| M2 | String overload uses same alignment result as DnaSequence overload | same as M1 | Same aligned sequences and score | **ASSUMPTION** (wrapper behavior from implementation) |
| M3 | Empty string input returns AlignmentResult.Empty | seq1 = "", seq2 = GATTACA | Empty result | **ASSUMPTION** (implementation behavior) |
| M4 | Null DnaSequence throws ArgumentNullException | seq1 = null, seq2 = GATTACA | ArgumentNullException | **ASSUMPTION** (implementation behavior) |

### 4.2 SHOULD Tests (Important edge cases)

None evidence-backed beyond MUST at this time.

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| C1 | Alternative optimal alignment paths | sequences with multiple optimal alignments | Any optimal alignment is acceptable | Evidence: multiple optimal alignments can exist; requires explicit dataset |

---

## 5. Audit of Existing Tests

### 5.1 Test Discovery Summary

- Canonical tests found in [Seqeron.Genomics.Tests/SequenceAlignerTests.cs](Seqeron.Genomics.Tests/SequenceAlignerTests.cs) (global alignment region).
- Wrapper smoke test for cancellation overload found in [Seqeron.Genomics.Tests/PerformanceExtensionsTests.cs](Seqeron.Genomics.Tests/PerformanceExtensionsTests.cs).

### 5.2 Coverage Classification (Pre-Consolidation)

| Area | Status | Notes |
|------|--------|-------|
| Global alignment invariants | Weak | Existing tests use minimal assertions (score > 0, gap presence) without evidence grounding |
| Wikipedia example dataset | Missing | No dataset-based verification |
| Empty input handling | Covered | Present but not evidence-based |
| Null argument handling | Covered | Present but not evidence-based |

### 5.3 Consolidation Plan

- **Canonical file:** new [Seqeron.Genomics.Tests/SequenceAligner_GlobalAlign_Tests.cs](Seqeron.Genomics.Tests/SequenceAligner_GlobalAlign_Tests.cs)
  - Deep, evidence-based invariants using Wikipedia example and scoring rules.
- **Wrapper/delegate smoke tests:**
  - Keep cancellation overload smoke test in [Seqeron.Genomics.Tests/PerformanceExtensionsTests.cs](Seqeron.Genomics.Tests/PerformanceExtensionsTests.cs).
  - Include single wrapper test for string overload in canonical file.
- **Remove** weak or duplicate global alignment tests from shared [Seqeron.Genomics.Tests/SequenceAlignerTests.cs](Seqeron.Genomics.Tests/SequenceAlignerTests.cs).

### 5.4 Final State After Consolidation

| File | Tests | Role |
|------|-------|------|
| [Seqeron.Genomics.Tests/SequenceAligner_GlobalAlign_Tests.cs](Seqeron.Genomics.Tests/SequenceAligner_GlobalAlign_Tests.cs) | 4 | Canonical + wrapper smoke |
| [Seqeron.Genomics.Tests/PerformanceExtensionsTests.cs](Seqeron.Genomics.Tests/PerformanceExtensionsTests.cs) | 1 (existing) | Wrapper smoke (cancellation overload) |

---

## 6. Open Questions / Decisions

None.
