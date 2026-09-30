# Test Specification: PAT-APPROX-002

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | PAT-APPROX-002 |
| **Title** | Approximate Matching (Edit Distance) |
| **Area** | Pattern Matching |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-09-30 |

---

## Canonical Methods Under Test

| Method | Class | Type | Complexity |
|--------|-------|------|------------|
| `EditDistance(string s1, string s2)` | ApproximateMatcher | Canonical | O(m × n) |
| `FindWithEdits(string sequence, string pattern, int maxEdits)` | ApproximateMatcher | Canonical | O(n × m × (m+k)) |
| `FindEditEndPositions(string sequence, string pattern, int maxEdits)` | ApproximateMatcher | Canonical (Sellers 1980) | O(n × m) |
| `FindWithEdits(DnaSequence sequence, string pattern, int maxEdits)` | ApproximateMatcher | Wrapper | Delegates to string version |
| `GetEditAlignment(string query, string target)` | ApproximateMatcher | Canonical (traceback, edlib CIGAR) | O(m × n) |
| `GetEditAlignmentLinearSpace(string query, string target)` | ApproximateMatcher | Canonical (Hirschberg 1975, linear space) | O(m × n) time, O(m + n) space |
| `OptimalStringAlignmentDistance(string s1, string s2)` | ApproximateMatcher | Canonical (restricted Damerau) | O(m × n) |
| `DamerauLevenshteinDistance(string s1, string s2)` | ApproximateMatcher | Canonical (Lowrance–Wagner 1975) | O(m × n) |

`EditDistance` and `FindEditEndPositions` run on the Myers (1999) bit-parallel engine (Hyyrö 2003 global form, ⌈m/64⌉ words, edlib `calculateBlock`); the Wagner–Fischer references `EditDistanceDp` / `FindEditEndPositionsDp` (internal) are the test oracles.

---

## Evidence Summary

### Sources Consulted
1. **Wikipedia - Levenshtein Distance**: Definition, mathematical formula, canonical examples
2. **Wikipedia - Edit Distance**: Properties, metric axioms, algorithm types
3. **Rosetta Code - Levenshtein Distance**: Test vectors, cross-language validation
4. **Navarro (2001)**: "A Guided Tour to Approximate String Matching" - theoretical foundation; §5.1 Sellers DP and `survey`/`surgery` example
5. **edlib 1.3 (infix/HW mode)**, rapidfuzz 3.14 Levenshtein: reference implementations used for the review-2026-09 cross-check
6. **Myers (1999) J. ACM 46(3):395; Hyyrö (2003); edlib source** (`edlib.cpp`: `calculateBlock`, `obtainAlignmentTraceback`, `edlibAlignmentToCigar`; `edlib.h`: EDLIB_EDOP_*, EDLIB_CIGAR_*) — bit-parallel engine and CIGAR convention (B05 follow-up, 2026-09-30)
7. **Lowrance & Wagner (1975) J. ACM 22(2):177; Damerau (1964); Boytsov (2011)**; rapidfuzz `OSA` / `DamerauLevenshtein`, jellyfish — Damerau variants

### B05 follow-up reference cross-check (2026-09-30)

| Check | Cases | Result |
|-------|-------|--------|
| `EditDistance` (Myers) vs rapidfuzz Levenshtein | 3160 random pairs, lengths 0–300, alphabets {AC, ACGT, ACGTN, a–z, ACαβ} | 3160/3160 equal |
| `GetEditAlignment` vs edlib NW `task='path'` | 2000 random pairs (≤25, {AC}/{ACGT}) | distance + CIGAR replay 2000/2000; CIGAR identical 116; identical on all 115 pairs with a unique optimal path |
| edlib tie-break reproduction (Python traceback I → D → diagonal) | 2000 pairs | 2000/2000 identical to edlib ⇒ differences are co-optimal ties only |
| `FindEditEndPositions` (Myers) vs Python Sellers DP; edlib HW best distance + end locations | 600 (text ≤300, pattern ≤150) | 600/600 |
| `FindWithEdits` windows, CIGAR replay, MismatchPositions, MismatchType | 500 cases, 13019 hits | 0 errors; 7461 CIGARs identical to edlib NW |
| `GetEditAlignmentLinearSpace` (Hirschberg) vs edlib NW `editDistance`, `GetEditAlignment`, `EditDistance` | 3501 pairs (3000 random ≤120 incl. non-ASCII + mutated copies, 500 exhaustive-small, one 3000×3300) | distance 3501/3501 equal; CIGAR replay 3501/3501; path identical to `GetEditAlignment` 1757/3501 (co-optimal ties ⇒ separate method) |
| OSA / DL vs rapidfuzz, DL vs jellyfish | 4507 pairs | 0 differences (CA/ABC: OSA 3, DL 2; pyxDamerauLevenshtein gives 3 = OSA) |

**Tie-break (documented, deterministic):** traceback from (m, n) takes the diagonal when optimal, then `I`, then `D` (edlib: `I`, `D`, diagonal). Diagonal-first makes a substitution-only optimum (equal-length window with ed = Hamming) come back as the Hamming path, so `MismatchPositions` of `Substitution` hits equal the Hamming mismatch indices.

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
| C02 | FindEditEndPositions_NavarroSurveySurgery_ReturnsSellersEnds | Sellers worked example (ends 4,5,6 at d=2) | Navarro 2001 §5.1; edlib |
| C03 | FindEditEndPositions_TtacInGattaca_MatchesEdlib | Sellers end positions | edlib HW |
| C04 | FindEditEndPositions_AcgaK2_MatchesEdlib | Sellers end positions, case-insensitive text | edlib HW |
| C05 | FindEditEndPositions_Guards | Negative k / empty input | Contract |
| C06 | FindWithEdits_NavarroSurveySurgery_ReturnsAllWindows | Window enumeration | Navarro 2001 §5.1; rapidfuzz |
| C07 | FindWithEdits_EndPositionSet_EqualsSellers | Window ends ≡ Sellers ends (300 seeded cases) | Sellers 1980 |
| C08 | GetEditAlignment_UniqueOptimalPath_EqualsEdlibCigar (4 cases) | CIGAR = edlib on unique optimal paths | edlib NW |
| C09 | GetEditAlignment_CoOptimalPaths_UsesDiagonalFirstTieBreak (3 cases) | Documented tie-break; edlib path same cost | edlib NW |
| C10 | GetEditAlignment_AlignedStringsAndSubstitutions / _EmptyInputs | Gapped strings, guards | edlib convention |
| C11 | FindWithEdits_SurveySurgery_AlignmentsEqualEdlib / _TtacInGattaca_Cigars | Per-hit CIGAR + MismatchPositions | edlib NW per window |
| C12 | FindWithEdits_HugeMaxEdits_DoesNotOverflow | maxEdits = int.MaxValue reports every window | k-differences definition |
| C13 | EditDistance_MultiBlock_EqualsRapidfuzz, _Myers_EqualsDp_Exhaustive, _NonAsciiSymbols_EqualDp, FindEditEndPositions_Myers_EqualsDp | Myers engine == DP / rapidfuzz, m > 64 | Myers 1999; rapidfuzz |
| C14 | DamerauVariants_EqualRapidfuzz (11 cases), DamerauVariants_NullInput_Throws | OSA / DL reference values | rapidfuzz, jellyfish |
| P01 | Properties/EditAlignmentProperties P1–P5 | Myers == DP; CIGAR replay cost = distance; hit alignments; DL ≤ OSA ≤ Lev; DL triangle | Myers 1999; Lowrance–Wagner 1975 |
| P02 | Metamorphic/PatternApproxB05MetamorphicTests.EditDistance_ReversingBothStrings_PreservesDistance | Reversal invariance | Definition |
| P03 | Fuzzing/PatternApproxEditFuzzTests | Word-boundary lengths 63–129, non-ASCII, int.MaxValue | Myers 1999 |

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
| C02 | FindEditEndPositions_NavarroSurveySurgery_ReturnsSellersEnds | ✅ Covered |
| C03 | FindEditEndPositions_TtacInGattaca_MatchesEdlib | ✅ Covered |
| C04 | FindEditEndPositions_AcgaK2_MatchesEdlib | ✅ Covered |
| C05 | FindEditEndPositions_Guards | ✅ Covered |
| C06 | FindWithEdits_NavarroSurveySurgery_ReturnsAllWindows | ✅ Covered |
| C07 | FindWithEdits_EndPositionSet_EqualsSellers | ✅ Covered |
| C08–C14 | ApproximateMatcher_EditAlignment_Tests | ✅ Covered |
| P01–P03 | Property / Metamorphic / Fuzzing (B05 follow-up) | ✅ Covered |

---

## Validation Criteria

- [x] All MUST tests pass (13/13)
- [x] All SHOULD tests pass (10/10)
- [x] Zero warnings in test file
- [x] Tests are deterministic
- [x] Tests follow NUnit conventions
- [x] Naming follows `Method_Scenario_ExpectedResult` pattern

### Linear-space alignment tests (B05 follow-up, 2026-09-30)

| ID | Test | Evidence |
|----|------|----------|
| LS-M1 | 14 literal pairs: distance = edlib NW, CIGAR replays with that cost, distance = `GetEditAlignment` (`ApproximateMatcher_EditAlignmentLinearSpace_Tests`) | edlib 1.3 |
| LS-M2 | co-optimal path may differ from `GetEditAlignment` (locked example, distance 5) | documented decision |
| LS-S1 | aligned strings / SubstitutionPositions / HasIndels / STANDARD CIGAR; null guards; 3000×3300 input | contract |
| LS-P6 | `Properties/EditAlignmentProperties.P6`: 1500 random pairs, distance = full traceback = Myers, replay valid | Hirschberg 1975 optimality |
| LS-MR | `Metamorphic/PatternApproxB05MetamorphicTests`: swap (I↔D script valid for swapped pair), reversal preserve distance | symmetry of ed |
| LS-F | `Fuzzing/PatternApproxEditFuzzTests`: extreme shapes (0×n, 1×700, 900×2, 1000×1000), non-ASCII/lone surrogates | robustness |
