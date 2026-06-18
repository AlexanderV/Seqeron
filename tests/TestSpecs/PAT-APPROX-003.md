# Test Specification: PAT-APPROX-003

**Test Unit ID:** PAT-APPROX-003
**Area:** Matching
**Algorithm:** Best Match and Frequency Analysis (Approximate Pattern Matching / Frequent Words with Mismatches)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | ROSALIND BA1I — Most Frequent Words with Mismatches (Compeau & Pevzner, ch. 1) | 1 | https://rosalind.info/problems/ba1i/ | 2026-06-13 |
| 2 | ROSALIND BA1H — Approximate Occurrences of a Pattern | 1 | https://rosalind.info/problems/ba1h/ | 2026-06-13 |
| 3 | ROSALIND BA1N — d-Neighborhood of a String | 1 | https://rosalind.info/problems/ba1n/ | 2026-06-13 |
| 4 | charlesreid1 go-rosalind `rosalind_ba1.go` | 3 | https://raw.githubusercontent.com/charlesreid1/go-rosalind/master/rosalind/rosalind_ba1.go | 2026-06-13 |
| 5 | zonghui0228 Rosalind-Solutions `rosalind_ba1h.py` | 3 | https://github.com/zonghui0228/Rosalind-Solutions/blob/master/code/rosalind_ba1h.py | 2026-06-13 |

### 1.2 Key Evidence Points

1. Count_d(Text, Pattern) = number of occurrences of Pattern in Text with at most d mismatches (Hamming) — Source 1.
2. An approximate occurrence at position i exists iff HammingDistance(Pattern, Text[i..i+k]) ≤ d — Source 2.
3. The most-frequent mismatch k-mer need not appear exactly in Text; counting is over the d-neighborhood — Source 1.
4. When several k-mers tie for the maximum mismatch-count, ALL are returned — Source 1 (sample returns three).
5. Neighbors(Pattern, d) = all k-mers within Hamming distance ≤ d, and always includes Pattern itself — Sources 3, 4.
6. Approximate matching is O(n·m) (Hamming distance per window) — Sources 4, 5.

### 1.3 Documented Corner Cases

- d = 0 degenerates Count_d to exact occurrence count and FrequentWords to exact frequent k-mers (Neighbors(P,0) = {P}) — Sources 1, 3.
- Pattern absent as exact substring is valid; matching/counting is over the Hamming ball — Sources 1, 2.
- Ties in BA1I must all be reported — Source 1.

### 1.4 Known Failure Modes / Pitfalls

1. Returning only one k-mer when multiple tie for the max count (BA1I) — Source 1.
2. Counting only exact-substring neighbors instead of the full d-neighborhood — Source 1.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindFrequentKmersWithMismatches(sequence, k, d)` | ApproximateMatcher | **Canonical** | BA1I; tally over d-neighborhoods, return all maxima |
| `CountApproximateOccurrences(sequence, pattern, maxMismatches)` | ApproximateMatcher | **Canonical** | Count_d (BA1H); delegates to FindWithMismatches |
| `FindBestMatch(sequence, pattern)` | ApproximateMatcher | **Canonical** | leftmost minimum-Hamming-distance equal-length window |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Count_d ≥ exact occurrence count; Count_0 = exact count | Yes | Source 1 (Count_d definition); Neighbors(P,0)={P} (Source 3) |
| INV-2 | A position i is an approximate occurrence iff HammingDistance(Pattern, Text[i..i+\|Pattern\|]) ≤ d | Yes | Source 2 |
| INV-3 | FindFrequentKmersWithMismatches returns every k-mer achieving the maximum Count_d (all ties) | Yes | Source 1 (sample has 3) |
| INV-4 | FindBestMatch returns distance = min over all equal-length windows of HammingDistance; IsExact iff that min = 0 | Yes | Source 2 (Hamming definition) |
| INV-5 | FindBestMatch tie-break: the leftmost window achieving the minimum distance is returned | Yes | **ASSUMPTION** (API convention; does not change the returned distance) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | FrequentKmers BA1I sample | Text=ACGTTGCATGTCGCATGATGCATGAGAGCT, k=4, d=1 | set {GATG, ATGC, ATGT}, count 5 each | Source 1 sample |
| M2 | FrequentKmers all ties returned | same as M1 — exactly 3 k-mers returned | Count == 3 | Source 1 (INV-3) |
| M3 | Count_d BA1H sample | Text (99 nt), Pattern=ATTCTGGA, d=3 | Count == 5 | Source 2 sample |
| M4 | Approximate positions BA1H sample | FindWithMismatches positions for M3 | {6,7,26,27,78} | Source 2 sample |
| M5 | Count_1 worked example | Text=AACAAGCTGATAAACATTTAAAGAG, Pattern=AAAAA, d=1 | Count == 4 | Source 1 worked example |
| M6 | Count_0 = exact | Text=ACGTACGT, Pattern=ACGT, d=0 | Count == 2 | Source 1 (INV-1) |
| M7 | FindBestMatch exact | seq=ACGTACGT, pat=ACGT | Distance 0, IsExact true, Position 0 | Source 2 (INV-4) |
| M8 | FindBestMatch no exact | seq=TTTTTTTT, pat=ACGT | Distance 3, leftmost Position 0, MatchedSequence TTTT | Source 2 (INV-4, INV-5) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | FrequentKmers d=0 exact | seq=AAAAAA, k=4, d=0 | single (AAAA, 3) | degenerate exact (INV-1) |
| S2 | FindBestMatch leftmost tie | seq=ACGTACGA, pat=TTTT | distance 3 at leftmost minimal window | INV-5 convention |
| S3 | Count case-insensitive | lowercase input matches uppercase | same count as uppercase | implementation upper-cases input |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | FindBestMatch empty/too-short | seq="", or pattern longer than seq | null | contract |
| C2 | FrequentKmers invalid k/d | k=0 or d<0 | ArgumentOutOfRangeException | contract |
| C3 | FrequentKmers empty seq | seq="" | empty result | contract |
| C4 | Count empty/short inputs | empty seq or pattern longer than seq | 0 | contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/ApproximateMatcherTests.cs` contains pre-template tests for exactly these three methods (FindBestMatch, CountApproximateOccurrences, FindFrequentKmersWithMismatches). No other file covers them.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (BA1I sample) | ❌ Missing | old file only tests trivial AAAAAA/ACGT; no published sample |
| M2 (all ties) | ❌ Missing | not tested |
| M3 (Count_d BA1H) | ❌ Missing | old file uses ad-hoc ACGTACGT only |
| M4 (positions BA1H) | ❌ Missing | not tested |
| M5 (Count_1 example) | ❌ Missing | not tested |
| M6 (Count_0 exact) | ⚠ Weak | `CountApproximateOccurrences_ExactMatches` exists but no assertion message / not evidence-cited |
| M7 (best exact) | ⚠ Weak | `FindBestMatch_ExactMatch` lacks messages, no Position check |
| M8 (best no exact) | ⚠ Weak | `FindBestMatch_NoExactMatch` lacks messages / Position / MatchedSequence |
| S1 (d=0 exact freq) | ⚠ Weak | `FindFrequentKmersWithMismatches_RepeatSequence` lacks messages |
| S2 (leftmost tie) | ❌ Missing | not tested |
| S3 (case-insensitive) | ❌ Missing | not tested |
| C1 (best empty/short) | ⚠ Weak | exists, no messages |
| C2 (invalid k/d) | ⚠ Weak | exists, no messages |
| C3 (freq empty) | ❌ Missing | not tested |
| C4 (count empty/short) | ❌ Missing | not tested |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ApproximateMatcher_FindBestMatch_Tests.cs` — all PAT-APPROX-003 cases (M/S/C), `#region` per method, evidence-cited.
- **Remove:** `tests/Seqeron/Seqeron.Genomics.Tests/ApproximateMatcherTests.cs` — its entire contents cover exactly these three methods with weak, non-evidence-based assertions; superseded by the canonical file (avoids duplication).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| ApproximateMatcher_FindBestMatch_Tests.cs | canonical PAT-APPROX-003 | 18 |
| ApproximateMatcherTests.cs | removed | 0 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ | implemented (BA1I sample) | ✅ Done |
| 2 | M2 | ❌ | implemented (3 ties) | ✅ Done |
| 3 | M3 | ❌ | implemented (Count_d BA1H) | ✅ Done |
| 4 | M4 | ❌ | implemented (positions) | ✅ Done |
| 5 | M5 | ❌ | implemented (Count_1) | ✅ Done |
| 6 | M6 | ⚠ | rewritten with evidence | ✅ Done |
| 7 | M7 | ⚠ | rewritten with Position+message | ✅ Done |
| 8 | M8 | ⚠ | rewritten with Position+MatchedSequence | ✅ Done |
| 9 | S1 | ⚠ | rewritten with messages | ✅ Done |
| 10 | S2 | ❌ | implemented | ✅ Done |
| 11 | S3 | ❌ | implemented | ✅ Done |
| 12 | C1 | ⚠ | rewritten | ✅ Done |
| 13 | C2 | ⚠ | rewritten | ✅ Done |
| 14 | C3 | ❌ | implemented | ✅ Done |
| 15 | C4 | ❌ | implemented | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | FrequentKmers_Ba1iSample_ReturnsThreeMostFrequentKmers |
| M2 | ✅ | FrequentKmers_Ba1iSample_ReturnsAllTiedKmers |
| M3 | ✅ | CountApproximateOccurrences_Ba1hSample_ReturnsFive |
| M4 | ✅ | FindWithMismatches_Ba1hSample_ReturnsExpectedPositions |
| M5 | ✅ | CountApproximateOccurrences_Count1WorkedExample_ReturnsFour |
| M6 | ✅ | CountApproximateOccurrences_ExactZeroMismatch_EqualsExactCount |
| M7 | ✅ | FindBestMatch_ExactMatch_ReturnsZeroDistanceAtPositionZero |
| M8 | ✅ | FindBestMatch_NoExactMatch_ReturnsLeftmostMinimumDistanceWindow |
| S1 | ✅ | FrequentKmers_ZeroMismatch_ReducesToExactFrequentKmer |
| S2 | ✅ | FindBestMatch_TiedMinima_ReturnsLeftmostWindow |
| S3 | ✅ | CountApproximateOccurrences_LowercaseInput_MatchesUppercase |
| C1 | ✅ | FindBestMatch_EmptyOrTooShort_ReturnsNull |
| C2 | ✅ | FrequentKmers_InvalidKOrD_Throws |
| C3 | ✅ | FrequentKmers_EmptySequence_ReturnsEmpty |
| C4 | ✅ | CountApproximateOccurrences_EmptyOrTooLongPattern_ReturnsZero |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | FindBestMatch returns the leftmost minimum-distance window (tie-break). Does NOT change the returned distance value, only which equal-distance window is reported; an API convention, not a correctness-affecting algorithm parameter. | INV-5, M8, S2 |

---

## 7. Open Questions / Decisions

1. Suffix-tree reuse: evaluated and rejected for approximate matching (SuffixTree supports exact match only); decision recorded in the algorithm doc §5.2. No open questions remain.
