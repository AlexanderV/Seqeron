# Test Specification: PAT-APPROX-002

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | PAT-APPROX-002 |
| **Title** | Approximate Matching (Edit Distance) |
| **Area** | Pattern Matching |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-03-01 |

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
| "saturday" | "sunday" | 3 | Wikipedia (matrix example), Rosetta Code |
| "" | "abc" | 3 | Wikipedia (definition) |
| "abc" | "" | 3 | Wikipedia (definition) |
| "flaw" | "lawn" | 2 | Wikipedia (bounds section) |
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
| M10 | EditDistance_CaseSensitive_DistinguishesCase | Standard definition: characters are distinct symbols | Wikipedia (definition: `head(a) = head(b)`) |
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
| S08 | EditDistance_Bounds_WithinExpectedRange | Distance bounds property | Wikipedia |
| S09 | FindWithEdits_WithDeletion_Found | Deletion matching | Definition |
| S10 | EditDistance_SleepFleeting_ReturnsFive | Canonical example | Rosetta Code |

### COULD Tests (Comprehensive)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| C01 | FindWithEdits_DnaSequenceOverload_DelegatesToStringVersion | Wrapper verification | Implementation |

---

## Coverage Classification

| ID | Test Name | Status |
|----|-----------|--------|
| M01 | EditDistance_IdenticalStrings_ReturnsZero | ✅ Covered |
| M02 | EditDistance_EmptyAndNonEmpty_ReturnsLength | ✅ Covered |
| M03 | EditDistance_KittenSitting_ReturnsThree | ✅ Covered |
| M04 | EditDistance_RosettacodeRaisethysword_ReturnsEight | ✅ Covered |
| M05 | EditDistance_Symmetry_CommutativeProperty | ✅ Covered |
| M06 | EditDistance_SingleSubstitution_ReturnsOne | ✅ Covered |
| M07 | EditDistance_SingleInsertion_ReturnsOne | ✅ Covered |
| M08 | EditDistance_SingleDeletion_ReturnsOne | ✅ Covered |
| M09 | EditDistance_NullInput_ThrowsArgumentNullException | ✅ Covered |
| M10 | EditDistance_CaseSensitive_DistinguishesCase | ✅ Covered |
| M11 | FindWithEdits_ExactMatch_Found | ✅ Covered |
| M12 | FindWithEdits_NegativeMaxEdits_ThrowsException | ✅ Covered |
| M13 | EditDistance_FlawLawn_ReturnsTwo | ✅ Covered |
| S01 | EditDistance_SaturdaySunday_ReturnsThree | ✅ Covered |
| S02 | EditDistance_StopTops_ReturnsTwo | ✅ Covered |
| S03 | EditDistance_BothEmpty_ReturnsZero | ✅ Covered |
| S04 | FindWithEdits_WithSubstitution_Found | ✅ Covered |
| S05 | FindWithEdits_WithInsertion_Found | ✅ Covered |
| S06 | FindWithEdits_EmptyInputs_ReturnsEmpty | ✅ Covered |
| S07 | EditDistance_TriangleInequality_Holds | ✅ Covered |
| S08 | EditDistance_Bounds_WithinExpectedRange | ✅ Covered |
| S09 | FindWithEdits_WithDeletion_Found | ✅ Covered |
| S10 | EditDistance_SleepFleeting_ReturnsFive | ✅ Covered |
| C01 | FindWithEdits_DnaSequenceOverload_DelegatesToStringVersion | ✅ Covered |

---

## Validation Criteria

- [x] All MUST tests pass (13/13)
- [x] All SHOULD tests pass (10/10)
- [x] Zero warnings in test file
- [x] Tests are deterministic
- [x] Tests follow NUnit conventions
- [x] Naming follows `Method_Scenario_ExpectedResult` pattern
