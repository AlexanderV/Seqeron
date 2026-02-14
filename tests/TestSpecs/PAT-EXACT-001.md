# Test Specification: PAT-EXACT-001

**Test Unit ID:** PAT-EXACT-001
**Area:** Pattern Matching
**Algorithm:** Exact Pattern Search (Suffix Tree)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Suffix tree | https://en.wikipedia.org/wiki/Suffix_tree | 2026-01-22 |
| Gusfield (1997) | Algorithms on Strings, Trees and Sequences | Reference |
| CP-Algorithms: Ukkonen's Algorithm | https://cp-algorithms.com/string/suffix-tree-ukkonen.html | 2026-01-22 |
| Rosalind: Finding a Motif in DNA | https://rosalind.info/problems/subs/ | 2026-01-22 |

### 1.2 Algorithm Description (Wikipedia/Gusfield)

**Exact Pattern Matching using Suffix Trees:**

Given a text T of length n and a pattern P of length m:
1. Build suffix tree for T in O(n) time
2. Search for pattern P by traversing tree edges matching P characters
3. All occurrences are found by enumerating leaves in the subtree below match point

**Complexity:**
- Search: O(m + z) where m = pattern length, z = number of occurrences
- Contains check: O(m)
- Count: O(m) when leaf counts are pre-computed

### 1.3 Edge Cases from Evidence

| Edge Case | Expected Behavior | Source |
|-----------|-------------------|--------|
| Empty pattern | Returns all positions (0..n-1) | Formal language theory: ε is a substring of every string at every position |
| Pattern not found | Returns empty collection | Standard |
| Pattern = entire text | Returns [0] | Gusfield |
| Overlapping occurrences | All positions returned | Rosalind example: "ATAT" in "GATATATGCATATACTT" → 2,4,10 |
| Pattern longer than text | Returns empty | Standard |
| Single character pattern | All occurrences of that character | Standard |
| Null pattern | ArgumentNullException | Implementation |

### 1.4 Rosalind Test Case (SUBS Problem)

From Rosalind bioinformatics platform:
- **Input:** s = "GATATATGCATATACTT", t = "ATAT"
- **Output:** 2, 4, 10 (1-indexed; our 0-indexed: 1, 3, 9)
- **Note:** Positions are overlapping occurrences

### 1.5 Known Test Strings from Literature

| Text | Pattern | Occurrences (0-indexed) | Source |
|------|---------|-------------------------|--------|
| "banana" | "ana" | [1, 3] | Wikipedia suffix tree |
| "banana" | "a" | [1, 3, 5] | Standard |
| "banana" | "na" | [2, 4] | Standard |
| "mississippi" | "issi" | [1, 4] | Gusfield |
| "mississippi" | "i" | [1, 4, 7, 10] | Standard |

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindAllOccurrences(string)` | SuffixTree | **Canonical** | Core implementation |
| `FindAllOccurrences(ReadOnlySpan<char>)` | SuffixTree | **Canonical** | Span overload |
| `Contains(string)` | SuffixTree | **Canonical** | Existence check |
| `Contains(ReadOnlySpan<char>)` | SuffixTree | **Canonical** | Span overload |
| `CountOccurrences(string)` | SuffixTree | **Canonical** | Count via LeafCount |
| `CountOccurrences(ReadOnlySpan<char>)` | SuffixTree | **Canonical** | Span overload |
| `FindMotif(DnaSequence, string)` | GenomicAnalyzer | Wrapper | Delegates to SuffixTree |
| `FindExactMotif(DnaSequence, string)` | MotifFinder | Wrapper | Delegates to SuffixTree |

---

## 3. Invariants

| ID | Invariant | Verifiable |
|----|-----------|------------|
| INV-1 | All returned positions are valid: 0 ≤ pos ≤ text.Length - pattern.Length | Yes |
| INV-2 | text[pos..pos+pattern.Length] == pattern for all returned positions | Yes |
| INV-3 | CountOccurrences == FindAllOccurrences.Count | Yes |
| INV-4 | Contains == (CountOccurrences > 0) | Yes |
| INV-5 | All substrings of text are found | Yes (exhaustive test) |
| INV-6 | Patterns not in text return empty | Yes |
| INV-7 | Empty text → no matches (except empty pattern) | Yes |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Null pattern throws | `null` | ArgumentNullException | Implementation contract |
| M2 | Empty pattern returns all positions | `""` | [0..n-1] for n-length text | Formal language theory (ε ⊆ every string) |
| M3 | Empty tree returns empty | tree(""), pattern("a") | [] | Standard |
| M4 | Single occurrence at start | tree("hello world"), "hello" | [0] | Standard |
| M5 | Single occurrence at end | tree("hello world"), "world" | [6] | Standard |
| M6 | Single occurrence in middle | tree("hello world"), "lo wo" | [3] | Standard |
| M7 | Pattern not found | tree("hello world"), "xyz" | [] | Standard |
| M8 | Pattern longer than text | tree("abc"), "abcdef" | [] | Standard |
| M9 | Pattern = full text | tree("abcdef"), "abcdef" | [0] | Gusfield |
| M10 | Multiple non-overlapping | tree("abcabc"), "abc" | [0, 3] | Standard |
| M11 | Overlapping occurrences | tree("aaaa"), "aa" | [0, 1, 2] | Standard |
| M12 | Banana classic test | tree("banana"), "ana" | [1, 3] | Wikipedia |
| M13 | Mississippi test | tree("mississippi"), "issi" | [1, 4] | Gusfield |
| M14 | Single character multiple | tree("abracadabra"), "a" | [0, 3, 5, 7, 10] | Standard |
| M15 | Rosalind SUBS test (adapted) | tree("GATATATGCATATACTT"), "ATAT" | [1, 3, 9] | Rosalind |
| M16 | Contains matches FindAll | any | Contains == !FindAll.IsEmpty | INV-4 |
| M17 | Count matches FindAll.Count | any | Count == FindAll.Count | INV-3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | Case sensitivity | tree("AbCd"), "abcd" | [] | Implementation (case-sensitive) |
| S2 | Special characters: spaces | tree("hello world"), "o w" | [4] | Standard |
| S3 | Special characters: newlines | tree("line1\nline2"), "1\nl" | [4] | Standard |
| S4 | Large text performance | 1000+ chars | completes quickly | Performance |
| S5 | Span overload matches string | any | same results | API consistency |
| S6 | All suffixes contained | any | Contains(suffix) == true | Suffix tree property |
| S7 | All substrings contained | any | Contains(substring) == true | Suffix tree property |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| C1 | DNA sequences | ATGC patterns | correct positions | Bioinformatics use case |
| C2 | Unicode safety | non-ASCII | defined behavior | Robustness |
| C3 | Very long patterns | 100+ chars | correct positions | Edge case |

---

## 5. Wrapper/Delegate Tests

### GenomicAnalyzer.FindMotif and MotifFinder.FindExactMotif

These are wrappers that delegate to `SuffixTree.FindAllOccurrences`. Minimal smoke tests only.

| ID | Test Case | Notes |
|----|-----------|-------|
| W1 | Single occurrence smoke | Verify delegation works |
| W2 | Case normalization | Wrappers normalize to uppercase |
| W3 | Empty motif returns empty | Different behavior from SuffixTree |

---

## 6. Coverage Classification

### FindAllOccurrencesTests.cs (19 methods / 21 test runs)

| # | Test Method | Spec | Classification |
|---|-------------|------|----------------|
| 1 | `FindAll_NullPattern_ThrowsArgumentNullException` | M1 | ✅ Covered |
| 2 | `FindAll_EmptyPattern_ReturnsAllPositions` | M2 | ✅ Strengthened: exact positions [0,1,2] |
| 3 | `FindAll_EmptyTree_ReturnsEmpty` | M3 | ✅ Covered |
| 4 | `FindAll_SingleOccurrence_AtPosition` ×3 | M4,M5,M6 | ✅ Parametrized: start/end/middle |
| 5 | `FindAll_MultipleOccurrences_ReturnsAllPositions` | M10 | ✅ Covered |
| 6 | `FindAll_OverlappingOccurrences_FindsAll` | M11 | ✅ Covered |
| 7 | `FindAll_Banana_FindsAllOccurrences` | M12 | ✅ Covered |
| 8 | `FindAll_Mississippi_IssiPattern` | M13 | ✅ Covered |
| 9 | `FindAll_RosalindSubs_FindsOverlappingDnaMotif` | M15 | ✅ Covered |
| 10 | `FindAll_NonExistent_ReturnsEmpty` | M7 | ✅ Covered |
| 11 | `FindAll_PatternLongerThanText_ReturnsEmpty` | M8 | ✅ Added |
| 12 | `FindAll_FullString_ReturnsZero` | M9 | ✅ Covered |
| 13 | `FindAll_SingleCharacter_FindsAllOccurrences` | M14 | ✅ Covered |
| 14 | `FindAll_AllSubstrings_MatchLinearSearch` | S7/INV-2,5 | ✅ Covered |
| 15 | `FindAll_WithSpaces_Works` | S2 | ✅ Covered |
| 16 | `FindAll_WithNewlines_Works` | S3 | ✅ Covered |
| 17 | `FindAll_SpanOverload_MatchesStringOverload` | S5 | ✅ Covered |
| 18 | `FindAll_SpanFromSlice_Works` | S5 | ✅ Covered |
| 19 | `FindAll_EmptySpan_ReturnsAllPositions` | M2/S5 | ✅ Strengthened: exact positions [0,1,2] |

### ContainsTests.cs (12 methods / 12 test runs)

| # | Test Method | Spec | Classification |
|---|-------------|------|----------------|
| 1 | `Contains_NullPattern_ThrowsArgumentNullException` | M1 | ✅ Covered |
| 2 | `Contains_EmptyPattern_ReturnsTrue` | M2 | ✅ Covered |
| 3 | `Contains_EmptyTree_AnyPattern_ReturnsFalse` | M3/INV-7 | ✅ Covered |
| 4 | `Contains_FullString_ReturnsTrue` | M9 | ✅ Covered |
| 5 | `Contains_NonExistentPatterns_ReturnsFalse` | M7/M8 | ✅ Covered |
| 6 | `Contains_AllSubstrings_ReturnsTrue` | S6,S7/INV-5 | ✅ Covered |
| 7 | `Contains_OverlappingPatterns_Works` | M11 | ✅ Covered |
| 8 | `Contains_RepeatingCharacter_Works` | — | ✅ Covered |
| 9 | `Contains_SpanOverload_MatchesStringOverload` | S5 | ✅ Covered |
| 10 | `Contains_SpanFromSlice_Works` | S5 | ✅ Covered |
| 11 | `Contains_SpanFromCharArray_Works` | S5 | ✅ Covered |
| 12 | `Contains_IsCaseSensitive` | S1 | ✅ Covered |

### CountOccurrencesTests.cs (14 methods / 14 test runs)

| # | Test Method | Spec | Classification |
|---|-------------|------|----------------|
| 1 | `Count_NullPattern_ThrowsArgumentNullException` | M1 | ✅ Covered |
| 2 | `Count_EmptyPattern_ReturnsTextLength` | M2 | ✅ Covered |
| 3 | `Count_EmptyTree_ReturnsZero` | M3 | ✅ Covered |
| 4 | `Count_SingleOccurrence_ReturnsOne` | M4/M5 | ✅ Covered |
| 5 | `Count_MultipleOccurrences_ReturnsCorrectCount` | M10 | ✅ Covered |
| 6 | `Count_OverlappingPatterns_CountsOverlaps` | M11 | ✅ Covered |
| 7 | `Count_Banana_CorrectCounts` | M12 | ✅ Covered |
| 8 | `Count_NonExistent_ReturnsZero` | M7 | ✅ Covered |
| 9 | `Count_PatternLongerThanText_ReturnsZero` | M8 | ✅ Added |
| 10 | `Count_MatchesFindAllCount` | M17/INV-3 | ✅ Covered |
| 11 | `Count_ManyOccurrences_Works` | S4 | ✅ Covered |
| 12 | `Count_SpanOverload_MatchesStringOverload` | S5 | ✅ Covered |
| 13 | `Count_SpanFromSlice_Works` | S5 | ✅ Covered |
| 14 | `Count_EmptySpan_ReturnsTextLength` | M2/S5 | ✅ Covered |

### Changes Applied

| Action | Details |
|--------|---------|
| ⚠ Strengthened (2) | `FindAll_EmptyPattern`, `FindAll_EmptySpan`: exact positions [0,1,2] instead of count-only |
| 🔁 Merged (3→1) | `FindAll_SingleOccurrence` + `FindAll_AtBeginning` + `FindAll_AtEnd` → parametrized `FindAll_SingleOccurrence_AtPosition` with M4(hello→0), M5(world→6), M6(lo wo→3) |
| 🔁 Removed (1) | `FindAll_IsLazyEnumerated`: IReadOnlyList is not lazy; test asserted nothing meaningful |
| 🔁 Removed (4) | EdgeCaseTests duplicates: `Contains_PatternLongerThanText`, `FindAll_PatternLongerThanText`, `Contains_PatternEqualsText`, `FindAll_PatternEqualsText` (covered by canonical files) |
| ❌ Added (2) | `FindAll_PatternLongerThanText_ReturnsEmpty` (M8), `Count_PatternLongerThanText_ReturnsZero` (M8) |

**Total: 45 canonical methods / 47 test runs** (was 50 methods / 50 test runs)

---

## 7. Design Decisions

| Decision | Value | Source |
|----------|-------|--------|
| Empty pattern behavior | Returns all positions [0..n-1] | Formal language theory: the empty string ε is a substring of every string at every position |
| Case sensitivity | Case-sensitive matching | Suffix tree operates on raw characters; wrappers normalize to uppercase |
| Span overloads | Identical semantics to string overloads | Zero-allocation equivalent, not a different operation |
