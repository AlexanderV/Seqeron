# Test Specification: GENOMIC-REPEAT-001

**Test Unit ID:** GENOMIC-REPEAT-001
**Area:** Analysis
**Algorithm:** Repeat Detection — Longest Repeated Substring (LRS) and all repeated substrings via suffix tree
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | CMU 15-451 Lecture #10: Suffix Trees and Arrays, §2.1 | 1 | https://www.cs.cmu.edu/~15451-f17/lectures/lec10-sufftree.pdf | 2026-06-13 |
| 2 | Wikipedia — Longest repeated substring problem | 4 | https://en.wikipedia.org/wiki/Longest_repeated_substring_problem | 2026-06-13 |
| 3 | GeeksforGeeks — Suffix Tree Application 3: Longest Repeated Substring | 3 | https://www.geeksforgeeks.org/dsa/suffix-tree-application-3-longest-repeated-substring/ | 2026-06-13 |
| 4 | Langmead (JHU) — Suffix Trees lecture notes (cites Gusfield 5.4) | 1 | https://www.cs.jhu.edu/~langmea/resources/lecture_notes/08_suffix_trees_v2.pdf | 2026-06-13 |

### 1.2 Key Evidence Points

1. Longest repeated substring = longest string `r` that occurs at least twice in `T`; found as the deepest internal node with ≥ 2 leaves under it (path-label depth = substring length). — Source 1 §2.1, corroborated by Sources 2, 3.
2. A substring occurring ≥ 2 times corresponds to an internal node (≥ 2 children/leaves); a substring occurring once is a leaf. — Source 1, Source 3.
3. Overlapping occurrences count: `AAAAAAAAAA` → `AAAAAAAAA` (pos 0,1 overlap); `ABABABA` → `ABABA` (pos 0,2 overlap). — Source 3.
4. When no substring repeats, there is no LRS (`ABCDEFG` → none). — Source 3.
5. Suffix-tree construction + query is linear, Θ(n). — Sources 1, 2, 3.
6. The suffix-tree repeat-finding application family is attributed to Gusfield ch. 5–7 ("Gusfield 5.4"). — Source 4.

### 1.3 Documented Corner Cases

- No repeat present → no LRS (empty result). — Source 3 (`ABCDEFG`).
- Overlapping occurrences are valid and counted. — Source 3 (`AAAAAAAAAA`, `ABABABA`).
- Ties: only "a longest repeated substring" is required; any equal-length winner is correct. — Source 2.

### 1.4 Known Failure Modes / Pitfalls

1. Requiring non-overlapping occurrences would wrongly reject `AAAAAAAAA` for `AAAAAAAAAA`. — Source 3.
2. Returning a leaf-path (occurs once) as a repeat is incorrect; only internal nodes (≥ 2 occurrences) qualify. — Source 1.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindLongestRepeat(DnaSequence)` | GenomicAnalyzer | **Canonical** | LRS via suffix tree deepest internal node. |
| `FindRepeats(DnaSequence, int minLength)` | GenomicAnalyzer | **Canonical** | All substrings occurring ≥ 2 times with length ≥ minLength. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned repeat occurs at least twice (`Count ≥ 2`). | Yes | Source 1 §2.1 ("occurs at least twice") |
| INV-2 | `RepeatInfo.Length == RepeatInfo.Sequence.Length` and equals the path depth. | Yes | Source 2 (depth = substring length) |
| INV-3 | Every position in `Positions` is a true 0-based occurrence start of `Sequence`. | Yes | Source 1 (leaves under node = occurrences) |
| INV-4 | `FindLongestRepeat` returns a sequence whose length ≥ every sequence length returned by `FindRepeats(.., 1)` (it is *a* longest). | Yes | Source 2 ("longest") |
| INV-5 | `FindRepeats` returns only substrings with `Length ≥ minLength`. | Yes | Method contract (minLength filter) |
| INV-6 | `Positions` is sorted ascending. | Yes | **ASSUMPTION** (output-shape convention) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | LRS Wikipedia example | `FindLongestRepeat("ATCGATCGA")` | Sequence `ATCGA`, Length 5, Count 2, Positions {0,4} | Source 2 (ATCGATCGA$→ATCGA) |
| M2 | LRS overlapping run | `FindLongestRepeat("AAAAAAAAAA")` | Sequence `AAAAAAAAA`, Length 9, Count 2, Positions {0,1} | Source 3 (AAAAAAAAAA→AAAAAAAAA) |
| M3 | LRS overlapping period-2 | `FindLongestRepeat("ATATATA")` | Sequence `ATATA`, Length 5, Count 2, Positions {0,2} | Source 3 (ABABABA→ABABA analog) |
| M4 | LRS no repeat | `FindLongestRepeat("ACGT")` | `RepeatInfo.None` / `IsEmpty` true | Source 3 (ABCDEFG→none analog) |
| M5 | LRS empty input | `FindLongestRepeat("")` | `RepeatInfo.None` / `IsEmpty` true | Definition (no substring twice in ε) |
| M6 | FindRepeats full enumeration | `FindRepeats("ACGTACGTTTTTACGT", 3)` | Exactly {ACGT@{0,4,12}, CGT@{1,5,13}, TACGT@{3,11}, TTT@{7,8,9}, TTTT@{7,8}} | Definition; each substring occurs ≥2× (verified against suffix tree) |
| M7 | FindRepeats all occur ≥ 2 & meet minLength | Every result of M6 has Count ≥ 2 and Length ≥ 3 | INV-1, INV-5 hold for all | Source 1 §2.1; method contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | minLength above all repeats | `FindRepeats("ACGTACGT", 5)` | Empty (no repeat ≥ length 5) | minLength filter upper boundary |
| S2 | minLength boundary inclusive | `FindRepeats("ACGTACGT", 4)` contains `ACGT` (len 4 = minLength) | `ACGT`@{0,4} present | INV-5 boundary (≥, not >) |
| S3 | Invariant property over results | All `RepeatInfo` from M6 satisfy INV-1..INV-3, INV-6 | All hold | Property test for O(n²) enumeration |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | minLength ≤ 0 degenerate | `FindRepeats("ACGTACGT", 0)` | Only substrings occurring ≥ 2× (no zero-length); `ACGT`@{0,4} present, none with Length 0 | Degenerate parameter robustness |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzerTests.cs` contains three weak repeat tests: `FindLongestRepeat_SimpleRepeat_FindsIt`, `FindLongestRepeat_NoRepeat_ReturnsEmpty`, `FindRepeats_MultipleRepeats_FindsAll`.
- No `GenomicAnalyzer_FindRepeats_Tests.cs` (the canonical `{Class}_{Method}_Tests.cs` form) existed before this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| GenomicAnalyzerTests.FindLongestRepeat_SimpleRepeat_FindsIt | ⚠ Weak | No positions/length checks, no assertion messages; subsumed by M1/M2/M3. Remove. |
| GenomicAnalyzerTests.FindLongestRepeat_NoRepeat_ReturnsEmpty | ⚠ Weak | No message; subsumed by M4. Remove. |
| GenomicAnalyzerTests.FindRepeats_MultipleRepeats_FindsAll | ⚠ Weak | Permissive `Any(...)`, no full set / positions; subsumed by M6/M7. Remove. |
| M1 LRS Wikipedia | ❌ Missing | Implement in canonical file. |
| M2 LRS overlapping run | ❌ Missing | Implement. |
| M3 LRS overlapping period-2 | ❌ Missing | Implement. |
| M4 LRS no repeat | ❌ Missing | Implement. |
| M5 LRS empty | ❌ Missing | Implement. |
| M6 FindRepeats full enumeration | ❌ Missing | Implement. |
| M7 FindRepeats count/minLength invariant | ❌ Missing | Implement. |
| S1 minLength above all | ❌ Missing | Implement. |
| S2 minLength inclusive boundary | ❌ Missing | Implement. |
| S3 invariant property | ❌ Missing | Implement. |
| C1 minLength ≤ 0 | ❌ Missing | Implement. |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindRepeats_Tests.cs` — all MUST/SHOULD/COULD cases for both `FindLongestRepeat` and `FindRepeats`.
- **Remove:** the three weak repeat tests in `GenomicAnalyzerTests.cs` (subsumed; keeping them would duplicate weakly). Motif/other tests in that file are untouched (different unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `GenomicAnalyzer_FindRepeats_Tests.cs` | Canonical for GENOMIC-REPEAT-001 | 11 |
| `GenomicAnalyzerTests.cs` | Other GenomicAnalyzer units (motif, etc.); 3 repeat tests removed | (repeat tests = 0) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact LRS + positions | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented overlapping run | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented period-2 overlap | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented no-repeat | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented empty input | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented full set + positions | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented count≥2 & minLength | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented minLength-above | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented inclusive boundary | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented invariant property | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented minLength≤0 | ✅ Done |
| 12 | weak×3 | ⚠ Weak | Removed from GenomicAnalyzerTests.cs | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `FindLongestRepeat_WikipediaExample_ReturnsAtcga` |
| M2 | ✅ Covered | `FindLongestRepeat_OverlappingRun_ReturnsLengthNine` |
| M3 | ✅ Covered | `FindLongestRepeat_OverlappingPeriodTwo_ReturnsAtata` |
| M4 | ✅ Covered | `FindLongestRepeat_NoRepeat_ReturnsNone` |
| M5 | ✅ Covered | `FindLongestRepeat_EmptySequence_ReturnsNone` |
| M6 | ✅ Covered | `FindRepeats_MinLengthThree_ReturnsExactSetWithPositions` |
| M7 | ✅ Covered | `FindRepeats_AllResults_OccurAtLeastTwiceAndMeetMinLength` |
| S1 | ✅ Covered | `FindRepeats_MinLengthAboveAllRepeats_ReturnsEmpty` |
| S2 | ✅ Covered | `FindRepeats_MinLengthEqualsRepeatLength_IncludesRepeat` |
| S3 | ✅ Covered | `FindRepeats_AllResults_SatisfyInvariants` |
| C1 | ✅ Covered | `FindRepeats_MinLengthZero_ReturnsNoZeroLengthRepeats` |
| weak×3 removed | ✅ Covered | Removed from `GenomicAnalyzerTests.cs` |

Total in-scope cases: 11 (+1 cleanup). ✅ = 11. No ❌/⚠ remain.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Tie-breaking among equal-length longest repeats is unspecified by sources (only "a" longest required). | INV-4 framing; no MUST asserts a unique winner among ties (all MUST inputs have a single longest). |
| 2 | `Positions` listed ascending (output-shape convention; the *set* is fixed). | INV-6, position assertions in M1/M2/M3/M6. |

---

## 7. Open Questions / Decisions

1. **Search reuse:** `FindLongestRepeat`/`FindRepeats` already build on the repository `SuffixTree` (`LongestRepeatedSubstring`, `FindAllOccurrences`, `GetAllSuffixes`). Decision: keep the suffix-tree implementation — it is exactly the algorithmically appropriate structure per Sources 1–3 (deepest internal node, many occurrence queries against one text). No naive scan introduced. Recorded in the algorithm doc §5.2.
