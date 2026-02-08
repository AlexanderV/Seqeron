# Test Specification: PAT-APPROX-001

**Test Unit ID:** PAT-APPROX-001
**Area:** Pattern Matching
**Algorithm:** Approximate Matching (Hamming Distance)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-01-22

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Hamming distance | https://en.wikipedia.org/wiki/Hamming_distance | 2026-01-22 |
| Wikipedia: Approximate string matching | https://en.wikipedia.org/wiki/Approximate_string_matching | 2026-01-22 |
| Rosalind: Counting Point Mutations (HAMM) | https://rosalind.info/problems/hamm/ | 2026-01-22 |
| Hamming, R.W. (1950) | Error detecting and error correcting codes. Bell System Technical Journal | Reference |
| Gusfield (1997) | Algorithms on Strings, Trees and Sequences | Reference |
| Navarro (2001) | A guided tour to approximate string matching. ACM Computing Surveys | Reference |

### 1.2 Algorithm Description

#### Hamming Distance (Wikipedia)

> The Hamming distance between two equal-length strings of symbols is the number
> of positions at which the corresponding symbols are different.

**Definition:** For two strings s and t of equal length n:
$$d_H(s, t) = \sum_{i=0}^{n-1} \mathbf{1}_{s[i] \neq t[i]}$$

**Properties (Wikipedia):**
- **Non-negativity:** d_H(s, t) ≥ 0
- **Identity:** d_H(s, t) = 0 if and only if s = t
- **Symmetry:** d_H(s, t) = d_H(t, s)
- **Triangle inequality:** d_H(s, u) ≤ d_H(s, t) + d_H(t, u)

**Constraint:** Strings must have equal length; undefined for unequal lengths.

#### Approximate Pattern Matching with k Mismatches (Wikipedia/Navarro)

> The specialized problem of pattern matching with k mismatches is defined by
> disallowing insertions and deletions in the matching of the pattern to the text.
> Thus, the Hamming distance of the pattern to the corresponding segment of the
> text has to be at most k.

**Complexity:** O(n × m) brute force, O(n√k log k) with advanced algorithms.

### 1.3 Reference Examples from Evidence

#### Wikipedia Hamming Distance Examples

| String 1 | String 2 | Distance | Source |
|----------|----------|----------|--------|
| "karolin" | "kathrin" | 3 | Wikipedia |
| "karolin" | "kerstin" | 3 | Wikipedia |
| "kathrin" | "kerstin" | 4 | Wikipedia |
| 0000 | 1111 | 4 | Wikipedia |

#### Rosalind HAMM Problem

**Input:**
- s = "GAGCCTACTAACGGGAT"
- t = "CATCGTAATGACGGCCT"

**Output:** 7

This is a canonical bioinformatics test case for Hamming distance between DNA sequences.

### 1.4 Edge Cases from Evidence

| Edge Case | Expected Behavior | Source |
|-----------|-------------------|--------|
| Identical strings | Distance = 0 | Wikipedia (identity property) |
| Completely different | Distance = length | Wikipedia |
| Unequal lengths | Exception/Error | Wikipedia (constraint) |
| Empty strings (equal) | Distance = 0 | Mathematical definition |
| Case variation | Implementation-defined (typically case-insensitive for DNA) | Implementation |
| maxMismatches = 0 | Exact matching only | Definition |
| maxMismatches ≥ pattern length | All positions match | Definition |
| Null input | ArgumentNullException | Implementation contract |
| Negative maxMismatches | ArgumentOutOfRangeException | Implementation contract |

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `HammingDistance(string, string)` | ApproximateMatcher | **Canonical** | Distance calculation |
| `HammingDistance(ReadOnlySpan<char>, ReadOnlySpan<char>)` | SequenceExtensions | **Canonical** | Span API |
| `FindWithMismatches(string, string, int)` | ApproximateMatcher | **Canonical** | Pattern matching |
| `FindWithMismatches(DnaSequence, string, int)` | ApproximateMatcher | Overload | DnaSequence support |
| `FindWithMismatches(..., CancellationToken)` | ApproximateMatcher | Variant | Cancellable |

---

## 3. Invariants

| ID | Invariant | Verifiable |
|----|-----------|------------|
| INV-1 | HammingDistance(s, t) ≥ 0 | Yes |
| INV-2 | HammingDistance(s, s) = 0 for any string s | Yes |
| INV-3 | HammingDistance(s, t) = HammingDistance(t, s) (symmetry) | Yes |
| INV-4 | HammingDistance(s, u) ≤ HammingDistance(s, t) + HammingDistance(t, u) (triangle) | Yes |
| INV-5 | HammingDistance requires equal lengths | Yes (exception) |
| INV-6 | FindWithMismatches(seq, pat, 0) ⊆ ExactMatches | Yes |
| INV-7 | For all results r: HammingDistance(seq[r.Position..], pat) ≤ maxMismatches | Yes |
| INV-8 | For all results r: r.Distance = HammingDistance(r.MatchedSequence, pattern) | Yes |
| INV-9 | Result positions are in range [0, seq.Length - pat.Length] | Yes |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

#### HammingDistance Tests

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Identical strings | "ACGT", "ACGT" | 0 | Wikipedia (identity) |
| M2 | One difference | "ACGT", "ACGG" | 1 | Definition |
| M3 | All different | "AAAA", "TTTT" | 4 | Definition |
| M4 | Case insensitive | "acgt", "ACGT" | 0 | Implementation |
| M5 | Unequal lengths throws | "ACGT", "ACG" | ArgumentException | Wikipedia (constraint) |
| M6 | Null s1 throws | null, "ACGT" | ArgumentNullException | Contract |
| M7 | Null s2 throws | "ACGT", null | ArgumentNullException | Contract |
| M8 | Empty strings | "", "" | 0 | Mathematical definition |
| M9 | Symmetry property | "ABCD", "AXYZ" | Same both ways | Wikipedia (symmetry) |
| M10 | Rosalind HAMM | "GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT" | 7 | Rosalind |
| M11 | Wikipedia example 1 | "karolin", "kathrin" | 3 | Wikipedia |
| M12 | Wikipedia example 2 | "karolin", "kerstin" | 3 | Wikipedia |

#### FindWithMismatches Tests

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M13 | Zero mismatches = exact | "ACGTACGT", "ACGT", 0 | [0, 4] | Definition |
| M14 | One mismatch found | "ACGTACGT", "ACGG", 1 | Positions with dist ≤ 1 | Definition |
| M15 | Too many mismatches not found | "ACGT", "TGCA", 2 | [] (dist=4 > 2) | Definition |
| M16 | Empty pattern returns empty | "ACGT", "", 1 | [] | Implementation |
| M17 | Empty sequence returns empty | "", "AC", 1 | [] | Standard |
| M18 | Pattern longer than sequence | "ACG", "ACGT", 1 | [] | Standard |
| M19 | Negative mismatches throws | "ACGT", "AC", -1 | ArgumentOutOfRangeException | Contract |
| M20 | maxMismatches = pattern length | "XXXX", "ACGT", 4 | [0] | All mismatches allowed |
| M21 | Invariant: result distance ≤ max | Any valid input | All r.Distance ≤ max | INV-7 |
| M22 | Invariant: result positions valid | Any valid input | All r.Position valid | INV-9 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | Single character pattern | "AAATAAAA", "A", 0 | All A positions | Standard |
| S2 | Pattern = sequence | "ACGT", "ACGT", 0 | [0] | Standard |
| S3 | Pattern = sequence with mismatches | "ACGT", "AXXX", 3 | [0] | Standard |
| S4 | Overlapping matches with mismatches | "AAAA", "AA", 1 | [0, 1, 2] | Standard |
| S5 | Large maxMismatches matches all | "ACGT", "XY", 10 | [0, 1, 2] (all positions) | Definition |
| S6 | DnaSequence overload works | DnaSequence("ACGT"), "AC", 0 | [0] | API |
| S7 | Span HammingDistance matches string | Same input | Same result | API consistency |
| S8 | MismatchPositions are correct | "ACGT", "AXGX", 2 | Positions [1, 3] | Implementation |
| S9 | Cancellation token respected | Long sequence | Completes or cancels | API |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| C1 | Triangle inequality | Three strings | d(a,c) ≤ d(a,b) + d(b,c) | Wikipedia |
| C2 | Large sequence performance | 10000+ chars | Completes in time | Performance |
| C3 | SNP detection use case | Genome-like sequence | Correct SNP locations | Bioinformatics |

---

## 5. ASSUMPTIONS

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Case-insensitive comparison for DNA | DNA sequences conventionally use case-insensitive matching |
| A2 | Empty pattern returns empty collection | Implementation choice; differs from some libraries |

---

## 6. Test Consolidation Plan

### Canonical Test File

- **Location:** `Seqeron.Genomics.Tests/ApproximateMatcher_HammingDistance_Tests.cs`
- **Contains:** All HammingDistance and FindWithMismatches MUST/SHOULD tests for PAT-APPROX-001
- **Status:** New file to separate Hamming-specific tests from Edit Distance tests

### Existing Tests to Migrate/Consolidate

From `ApproximateMatcherTests.cs`:
- HammingDistance tests (region) → Move to new canonical file
- FindWithMismatches tests (region) → Move to new canonical file
- Keep Edit Distance tests in ApproximateMatcherTests.cs (scope: PAT-APPROX-002)

### Smoke Tests (Wrapper/Span API)

- **Location:** `Seqeron.Genomics.Tests/PerformanceExtensionsTests.cs`
- **Status:** Contains `HammingDistance_Span_ReturnsCorrectDistance` - retain as smoke test

### Actions

1. Create new canonical file `ApproximateMatcher_HammingDistance_Tests.cs`
2. Migrate Hamming-related tests from `ApproximateMatcherTests.cs`
3. Add missing evidence-based tests (Rosalind HAMM, Wikipedia examples)
4. Add invariant tests (symmetry, triangle inequality)
5. Keep minimal smoke test in PerformanceExtensionsTests.cs

---

## 7. Audit Results

### Coverage Analysis of Existing Tests

| Category | Status | Notes |
|----------|--------|-------|
| Identical strings | ✅ Covered | `HammingDistance_IdenticalStrings_ReturnsZero` |
| One difference | ✅ Covered | `HammingDistance_OneDifference_ReturnsOne` |
| All different | ✅ Covered | `HammingDistance_AllDifferent_ReturnsLength` |
| Case insensitive | ✅ Covered | `HammingDistance_CaseInsensitive` |
| Unequal lengths throws | ✅ Covered | `HammingDistance_DifferentLengths_ThrowsException` |
| Empty strings | ✅ Covered | Added |
| Null input | ✅ Covered | Added |
| Rosalind HAMM example | ✅ Covered | Added |
| Wikipedia examples | ✅ Covered | Added |
| Symmetry invariant | ✅ Covered | Added |
| Zero mismatches = exact | ✅ Covered | `FindWithMismatches_ExactMatch_FoundWithZeroMismatches` |
| One mismatch found | ✅ Covered | `FindWithMismatches_OneMismatch_Found` |
| Too many mismatches | ✅ Covered | `FindWithMismatches_TooManyMismatches_NotFound` |
| Empty pattern | ✅ Covered | `FindWithMismatches_EmptyPattern_ReturnsEmpty` |
| Pattern too long | ✅ Covered | `FindWithMismatches_PatternLongerThanSequence_ReturnsEmpty` |
| Negative mismatches | ✅ Covered | `FindWithMismatches_NegativeMismatches_ThrowsException` |
| DnaSequence overload | ✅ Covered | `FindWithMismatches_DnaSequence_Works` |
| MismatchPositions correct | ✅ Covered | Explicit verification added |
| Result invariants | ✅ Covered | Invariant tests added |
| Span API | ✅ Covered | In PerformanceExtensionsTests.cs (smoke) |

### Weak/Redundant Tests

- `FindWithMismatches_MultipleMismatches_AllReturned` - Weak assertion, tests for exact match existence only
- Real-world tests (`SNP_Detection`, `PrimerBinding`) - Good but could use more precise assertions

### Duplicates Found

- None across files

---

## 8. Open Questions / Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| Separate Hamming from Edit tests? | Yes | Different Test Units (PAT-APPROX-001 vs PAT-APPROX-002) |
| Test triangle inequality? | Yes (COULD) | Mathematical property from Wikipedia |
| Include real-world genomics tests? | Keep existing | Good integration examples |

---
