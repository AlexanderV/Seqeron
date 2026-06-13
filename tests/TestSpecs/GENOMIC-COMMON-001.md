# Test Specification: GENOMIC-COMMON-001

**Test Unit ID:** GENOMIC-COMMON-001
**Area:** Analysis
**Algorithm:** Longest Common Substring / Common Region Detection (generalized suffix tree)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Gusfield (1997), Algorithms on Strings, Trees and Sequences (ISBN 0-521-58519-8) | 1 | citation retrieved from https://en.wikipedia.org/wiki/Longest_common_substring | 2026-06-13 |
| 2 | Wikipedia — Longest common substring | 4 | https://en.wikipedia.org/wiki/Longest_common_substring | 2026-06-13 |
| 3 | GeeksforGeeks — Suffix Tree Application 5 (LCS) | 3 | https://www.geeksforgeeks.org/dsa/suffix-tree-application-5-longest-common-substring-2/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. LCS = "a longest string which is substring of both S and T" — Wikipedia (source 2).
2. The common region is **contiguous**, unlike a subsequence — Wikipedia (source 2).
3. Found via a generalized suffix tree: deepest node whose subtree has leaves from both strings; Θ(n+m) — Gusfield 1997 / Wikipedia (sources 1, 2); GeeksforGeeks O(M+N) (source 3).
4. Ties exist: "BADANAT"/"CANADAS" share two maximal substrings "ADA" and "ANA" — Wikipedia (source 2).
5. Worked length example: "xabxac"/"abcabxabcd" → "abxa", length 4 — GeeksforGeeks (source 3).

### 1.3 Documented Corner Cases

- Multiple maximal substrings (ties): a deterministic representative must be chosen — Wikipedia (source 2).
- No common character → only the empty string qualifies → length-0 LCS — derived from the Wikipedia definition (source 2).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing longest common **substring** (contiguous) with longest common **subsequence** (gapped) — Wikipedia (source 2). The repository XML doc previously mislabeled the method "subsequence"; the implementation computes substring.
2. Nondeterministic tie-break — must return a documented, deterministic representative — Wikipedia (source 2) + SuffixTree XML doc.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindLongestCommonRegion(seq1, seq2)` | GenomicAnalyzer | Canonical | longest contiguous common substring with positions |
| `FindCommonRegions(seq1, seq2, minLength)` | GenomicAnalyzer | Canonical | all distinct common substrings of length ≥ minLength |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | The returned substring is contiguous and occurs in **both** sequences at the reported 0-based positions | Yes | Wikipedia definition (source 2) |
| INV-2 | No common contiguous substring strictly longer than the returned one exists | Yes | Wikipedia definition (source 2) |
| INV-3 | Empty / no-shared-character input → `CommonRegion.None` (empty, length 0, positions −1) | Yes | Wikipedia definition (source 2) |
| INV-4 | Every region from `FindCommonRegions(minLength)` has length ≥ minLength and is contained in sequence1 | Yes | definition extended to all common substrings (source 2) |
| INV-5 | On ties, the representative is "first found in sequence2" (deterministic) | Yes | **ASSUMPTION** (SuffixTree XML doc; see §6) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Unique LCS | `ACGTACGT` vs `TTACGTGG` | `TACGT`, length 5, pos1=3, pos2=1 | Wikipedia def (source 2); brute-force-verified |
| M2 | Tie | `CACAGAG` vs `TACATAGAT` | `ACA` (first-in-other of {ACA,AGA}), length 3, pos1=1, pos2=1 | Wikipedia tie property (source 2); SuffixTree tie-break |
| M3 | Identical | `ACGT` vs `ACGT` | `ACGT`, length 4, pos1=0, pos2=0 | Wikipedia def (source 2) |
| M4 | No common substring | `AAAA` vs `GGGG` | `CommonRegion.None` (empty, length 0, IsEmpty, pos −1/−1) | Wikipedia def (source 2) |
| M5 | All regions, minLength=4 | `FindCommonRegions(ACGTACGT, TTACGTGG, 4)` | `{(TACGT,3,1),(ACGT,0,2)}` | def of all common substrings (source 2); brute-force-verified |
| M6 | All regions, minLength=3 | `FindCommonRegions(ACGTACGT, TTACGTGG, 3)` | `{(TACGT,3,1),(ACGT,0,2),(CGT,1,3)}` | def (source 2); brute-force-verified |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty sequence1 | `FindLongestCommonRegion("", ACGT)` | `CommonRegion.None` | empty handling (INV-3) |
| S2 | Empty sequence2 | `FindLongestCommonRegion(ACGT, "")` | `CommonRegion.None` | empty handling (INV-3) |
| S3 | `FindCommonRegions` no match | `FindCommonRegions(AAAA, GGGG, 2)` | empty enumeration | INV-4 |
| S4 | Single-char overlap | `FindLongestCommonRegion("A", "TTTAT")` | `A`, length 1, pos1=0, pos2=3 | minimal non-empty LCS |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV-1 property | returned substring occurs at reported positions in both | substring matches both slices | property test |
| C2 | INV-4 property | all FindCommonRegions results length ≥ minLength and ⊆ seq1 | all satisfy | O(n·m) property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `FindLongestCommonRegion` / `FindCommonRegions`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit, no prior test |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindCommonRegion_Tests.cs` — all cases for this unit.
- **Remove:** none (no pre-existing tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `GenomicAnalyzer_FindCommonRegion_Tests.cs` | canonical | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | S1 | ❌ Missing | implemented | ✅ Done |
| 8 | S2 | ❌ Missing | implemented | ✅ Done |
| 9 | S3 | ❌ Missing | implemented | ✅ Done |
| 10 | S4 | ❌ Missing | implemented | ✅ Done |
| 11 | C1 | ❌ Missing | implemented | ✅ Done |
| 12 | C2 | ❌ Missing | implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `FindLongestCommonRegion_UniqueLcs_ReturnsTacgt` |
| M2 | ✅ Covered | `FindLongestCommonRegion_Tie_ReturnsFirstInOther` |
| M3 | ✅ Covered | `FindLongestCommonRegion_Identical_ReturnsWholeSequence` |
| M4 | ✅ Covered | `FindLongestCommonRegion_NoCommon_ReturnsNone` |
| M5 | ✅ Covered | `FindCommonRegions_MinLengthFour_ReturnsTwoRegions` |
| M6 | ✅ Covered | `FindCommonRegions_MinLengthThree_ReturnsThreeRegions` |
| S1 | ✅ Covered | `FindLongestCommonRegion_EmptyFirst_ReturnsNone` |
| S2 | ✅ Covered | `FindLongestCommonRegion_EmptySecond_ReturnsNone` |
| S3 | ✅ Covered | `FindCommonRegions_NoMatch_ReturnsEmpty` |
| S4 | ✅ Covered | `FindLongestCommonRegion_SingleCharOverlap_ReturnsA` |
| C1 | ✅ Covered | `FindLongestCommonRegion_ReturnedSubstring_OccursInBothAtPositions` |
| C2 | ✅ Covered | `FindCommonRegions_AllRegions_SatisfyMinLengthAndContainment` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Tie-break = "first found in sequence2" (deterministic, documented in SuffixTree XML doc; does not change maximal length, only the representative) | INV-5, M2 |

---

## 7. Open Questions / Decisions

1. The repository method `FindLongestCommonRegion` computes the longest common **substring** (contiguous) but its XML doc previously said "subsequence". Corrected in Phase 5; behavior unchanged (the SuffixTree call already computed substring).
