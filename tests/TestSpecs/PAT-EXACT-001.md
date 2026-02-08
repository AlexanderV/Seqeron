# Test Specification: PAT-EXACT-001

**Test Unit ID:** PAT-EXACT-001
**Area:** Pattern Matching
**Algorithm:** Exact Pattern Search (Suffix Tree)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-01-22

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
| Empty pattern | Returns all positions (0..n-1) OR empty, implementation-defined | Wikipedia: substring definition |
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
| M2 | Empty pattern returns all positions | `""` | [0..n-1] for n-length text | Regex convention |
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

## 6. Test Consolidation Plan

### Canonical Test File
- **Location:** `SuffixTree.Tests/Search/FindAllOccurrencesTests.cs`
- **Contains:** All MUST and SHOULD tests for FindAllOccurrences
- **Status:** Existing file with good coverage; needs Rosalind test case

### Related Canonical Files
- **ContainsTests.cs:** Tests for Contains method (complete)
- **CountOccurrencesTests.cs:** Tests for CountOccurrences method (complete)

### Wrapper Smoke Tests
- **Location:** `Seqeron.Genomics.Tests/MotifFinderTests.cs`
- **Status:** Contains smoke tests; adequate for wrappers

### Consolidation Actions
1. ✅ FindAllOccurrencesTests.cs - comprehensive, add Rosalind test
2. ✅ ContainsTests.cs - comprehensive, no changes needed
3. ✅ CountOccurrencesTests.cs - comprehensive, no changes needed
4. ✅ MotifFinderTests.cs - smoke tests exist, no changes needed

---

## 7. Audit Results

### Coverage Analysis

| Category | Status | Notes |
|----------|--------|-------|
| Null input | ✅ Covered | Tests exist |
| Empty pattern | ✅ Covered | Tests exist |
| Empty tree | ✅ Covered | Tests exist |
| Single occurrence | ✅ Covered | Start, end, middle |
| Multiple non-overlapping | ✅ Covered | Tests exist |
| Overlapping occurrences | ✅ Covered | "aaaa"/"aa" test |
| Pattern not found | ✅ Covered | Tests exist |
| Full text match | ✅ Covered | Tests exist |
| Banana test (Wikipedia) | ✅ Covered | Tests exist |
| Mississippi test (Gusfield) | ✅ Covered | "issi" pattern added |
| Rosalind SUBS test | ✅ Covered | Added |
| Exhaustive substring test | ✅ Covered | Tests exist |
| Span overload consistency | ✅ Covered | Tests exist |
| Contains/Count consistency | ✅ Covered | Tests exist |

### Duplicates Found
- None significant; tests are well-organized

### Missing Tests
~~All gaps closed:~~
1. ~~Rosalind SUBS problem example~~ ✅ Added
2. ~~Mississippi "issi" pattern~~ ✅ Added

---

## 8. Open Questions / Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| Empty pattern behavior | Returns all positions | Consistent with regex "" behavior |
| Case sensitivity | Case-sensitive | Canonical SuffixTree behavior; wrappers normalize |
| Span empty pattern | Returns empty (differs from string) | Implementation note documented |

---

## 9. Assumptions

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Tests use "banana" as canonical test string | Wikipedia suffix tree article uses it |
| A2 | Overlapping matches are all reported | Rosalind and standard string search semantics |
