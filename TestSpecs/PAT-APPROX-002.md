# Test Specification: PAT-APPROX-002

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | PAT-APPROX-002 |
| **Title** | Approximate Matching (Edit Distance) |
| **Area** | Pattern Matching |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-01-22 |

---

## Canonical Methods Under Test

| Method | Class | Type | Complexity |
|--------|-------|------|------------|
| `EditDistance(string s1, string s2)` | ApproximateMatcher | Canonical | O(m × n) |
| `FindWithEdits(string sequence, string pattern, int maxEdits)` | ApproximateMatcher | Canonical | O(n × m²) |
| `FindWithEdits(DnaSequence sequence, string pattern, int maxEdits)` | ApproximateMatcher | Wrapper | Delegates to string version |

---

## Evidence Summary

### Sources Consulted
1. **Wikipedia - Levenshtein Distance**: Definition, mathematical formula, canonical examples
2. **Wikipedia - Edit Distance**: Properties, metric axioms, algorithm types
3. **Rosetta Code - Levenshtein Distance**: Test vectors, cross-language validation
4. **Navarro (2001)**: "A Guided Tour to Approximate String Matching" - theoretical foundation

### Canonical Test Vectors (from Sources)

| String 1 | String 2 | Expected Distance | Source |
|----------|----------|-------------------|--------|
| "kitten" | "sitting" | 3 | Wikipedia, Rosetta Code |
| "rosettacode" | "raisethysword" | 8 | Rosetta Code |
| "saturday" | "sunday" | 3 | Rosetta Code |
| "" | "abc" | 3 | Wikipedia (definition) |
| "abc" | "" | 3 | Wikipedia (definition) |
| "flaw" | "lawn" | 2 | Wikipedia |
| "stop" | "tops" | 2 | Rosetta Code |
| "sleep" | "fleeting" | 5 | Rosetta Code |

### Invariants (from Sources)
1. **Symmetry:** EditDistance(a, b) == EditDistance(b, a)
2. **Identity:** EditDistance(a, a) == 0
3. **Empty string:** EditDistance("", s) == length(s)
4. **Triangle inequality:** EditDistance(a, c) ≤ EditDistance(a, b) + EditDistance(b, c)
5. **Bounds:** |len(a) - len(b)| ≤ EditDistance(a, b) ≤ max(len(a), len(b))

---

## Test Classification

### MUST Tests (Evidence-Backed)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| M01 | EditDistance_IdenticalStrings_ReturnsZero | Identity property | Wikipedia |
| M02 | EditDistance_EmptyAndNonEmpty_ReturnsLength | Base case definition | Wikipedia |
| M03 | EditDistance_KittenSitting_ReturnsThree | Canonical example | Wikipedia, Rosetta Code |
| M04 | EditDistance_RosettacodeRaisethysword_ReturnsEight | Canonical example | Rosetta Code |
| M05 | EditDistance_Symmetry_CommutativeProperty | Metric property | Wikipedia |
| M06 | EditDistance_SingleSubstitution_ReturnsOne | Substitution operation | Definition |
| M07 | EditDistance_SingleInsertion_ReturnsOne | Insertion operation | Definition |
| M08 | EditDistance_SingleDeletion_ReturnsOne | Deletion operation | Definition |
| M09 | EditDistance_NullInput_ThrowsArgumentNullException | Error handling | Implementation contract |
| M10 | EditDistance_CaseInsensitive_IgnoresCase | Implementation behavior | Implementation |
| M11 | FindWithEdits_ExactMatch_Found | maxEdits=0 behavior | Definition |
| M12 | FindWithEdits_NegativeMaxEdits_ThrowsException | Error handling | Implementation contract |
| M13 | EditDistance_FlawLawn_ReturnsTwo | Levenshtein < Hamming example | Wikipedia |

### SHOULD Tests (Good Practice)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| S01 | EditDistance_SaturdaySunday_ReturnsThree | Additional canonical case | Rosetta Code |
| S02 | EditDistance_StopTops_ReturnsTwo | Transposition-like case | Rosetta Code |
| S03 | EditDistance_BothEmpty_ReturnsZero | Edge case | Definition |
| S04 | FindWithEdits_WithSubstitution_Found | Substitution matching | Definition |
| S05 | FindWithEdits_WithInsertion_Found | Insertion matching | Definition |
| S06 | FindWithEdits_EmptyInputs_ReturnsEmpty | Edge case | Implementation |
| S07 | EditDistance_TriangleInequality_Holds | Metric property | Wikipedia |

### COULD Tests (Comprehensive)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| C01 | FindWithEdits_DnaSequenceOverload_DelegatesToStringVersion | Wrapper verification | Implementation |
| C02 | EditDistance_LongerStrings_CorrectDistance | Performance case | ASSUMPTION |

---

## Consolidation Plan

### Current State
- Existing tests in: `ApproximateMatcherTests.cs` (Edit Distance region, lines 19-96)
- 7 existing tests for EditDistance
- 3 existing tests for FindWithEdits

### Target State
- Create dedicated file: `ApproximateMatcher_EditDistance_Tests.cs`
- Deep, evidence-based tests for canonical algorithms
- Remove tests from `ApproximateMatcherTests.cs` (Edit Distance region)

### Migration
| Source | Target | Action |
|--------|--------|--------|
| ApproximateMatcherTests.cs (Edit Distance region) | ApproximateMatcher_EditDistance_Tests.cs | Move + enhance |
| ApproximateMatcherTests.cs (Find With Edits region) | ApproximateMatcher_EditDistance_Tests.cs | Move + enhance |

---

## Audit of Existing Tests

| Existing Test | Coverage | Assessment | Action |
|---------------|----------|------------|--------|
| EditDistance_IdenticalStrings_ReturnsZero | M01 | Covered | Keep (enhance) |
| EditDistance_OneSubstitution_ReturnsOne | M06 | Covered | Keep |
| EditDistance_OneInsertion_ReturnsOne | M07 | Covered | Keep |
| EditDistance_OneDeletion_ReturnsOne | M08 | Covered | Keep |
| EditDistance_EmptyAndNonEmpty_ReturnsLength | M02 | Covered | Keep |
| EditDistance_ComplexCase_CalculatesCorrectly | Partial | Weak (not canonical) | Replace with canonical vectors |
| EditDistance_CaseInsensitive | M10 | Covered | Keep |
| FindWithEdits_ExactMatch_Found | M11 | Covered | Keep |
| FindWithEdits_WithInsertion_Found | S05 | Weak | Enhance |
| FindWithEdits_WithDeletion_Found | S04 | Weak | Enhance |

### Missing Tests
- M03, M04, M05, M09, M12, M13 (canonical vectors, symmetry, error handling)
- S01, S02, S03, S07 (additional canonical cases, triangle inequality)

---

## Open Questions / Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| Should we test all Rosetta Code vectors? | Yes, for M03, M04, S01, S02 | Canonical cross-language validation |
| Keep complex GATTACA test? | No | Replace with Wikipedia-sourced "flaw"→"lawn" example |

---

## ASSUMPTIONS

| ID | Assumption | Justification |
|----|------------|---------------|
| A01 | Case-insensitive comparison is correct behavior | Implementation choice, documented |

---

## Validation Criteria

- [ ] All MUST tests pass
- [ ] All SHOULD tests pass  
- [ ] Zero warnings in test file
- [ ] Tests are deterministic
- [ ] Tests follow NUnit conventions
- [ ] Naming follows `Method_Scenario_ExpectedResult` pattern
