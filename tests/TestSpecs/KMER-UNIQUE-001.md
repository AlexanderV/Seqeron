# Test Specification: KMER-UNIQUE-001

**Test Unit ID:** KMER-UNIQUE-001
**Area:** K-mer
**Algorithm:** Unique K-mers / K-mers with Minimum Count (frequency filtering)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wikipedia — K-mer | 4 | https://en.wikipedia.org/wiki/K-mer | 2026-06-14 |
| 2 | BioInfoLogics — k-mer counting, part I | 3 | https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/ | 2026-06-14 |
| 3 | Compeau & Pevzner — Bioinformatics Algorithms (2nd ed.) | 1 | https://www.amazon.com/BIOINFORMATICS-ALGORITHMS-Phillip-Compeau/dp/0990374637 | 2026-06-14 |

### 1.2 Key Evidence Points

1. A k-mer is a substring of length k; a sequence of length L has **L − k + 1** total (overlapping) k-mers — Source 1.
2. **Unique k-mers** = k-mers that appear exactly once (frequency = 1); **distinct** = each different k-mer counted once — Source 2.
3. Worked table (ATCGATCAC, k=3): 7 total, 6 distinct, 5 unique = {TCG, CGA, GAT, TCA, CAC}; ATC has Count 2 and is excluded — Source 2.
4. Count(Text, Pattern) = number of overlapping occurrences; min-count filtering selects k-mers with Count ≥ t (recurrent k-mers) — Source 3.

### 1.3 Documented Corner Cases

- k > L ⇒ L − k + 1 ≤ 0 ⇒ zero k-mers (Source 1).
- A k-mer with Count ≥ 2 is distinct but NOT unique (Source 2, ATC example).

### 1.4 Known Failure Modes / Pitfalls

1. Counting non-overlapping occurrences instead of overlapping ones — would undercount repeats — Source 1/3.
2. Treating "distinct" as "unique" — would wrongly include repeated k-mers in the unique set — Source 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindUniqueKmers(string sequence, int k)` | KmerAnalyzer | Canonical | Returns k-mers with Count = 1 |
| `FindKmersWithMinCount(string sequence, int k, int minCount)` | KmerAnalyzer | Canonical | Returns (k-mer, Count) with Count ≥ minCount, ordered by Count descending |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every k-mer returned by `FindUniqueKmers` has Count exactly 1 in the sequence | Yes | Source 2 (unique = appears once) |
| INV-2 | `FindUniqueKmers` ⊆ distinct k-mers; |unique| ≤ |distinct| ≤ L − k + 1 | Yes | Source 1, 2 |
| INV-3 | Every pair returned by `FindKmersWithMinCount` has Count ≥ minCount and the reported Count equals the true overlapping occurrence count | Yes | Source 3 (Count ≥ t) |
| INV-4 | `FindKmersWithMinCount` output is ordered by Count in non-increasing order | Yes | Implementation contract (recurrent-first); consistent with Source 3 most-frequent ranking |
| INV-5 | `FindKmersWithMinCount(seq, k, 1)` set of keys equals the set of distinct k-mers | Yes | Source 2 (distinct), Source 3 (Count ≥ 1) |
| INV-6 | k > L or empty sequence ⇒ both methods return empty | Yes | Source 1 (L − k + 1 ≤ 0) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | FindUniqueKmers ATCGATCAC k=3 | Published worked table | {TCG, CGA, GAT, TCA, CAC}; ATC excluded | Source 2 table |
| M2 | FindUniqueKmers AGAT k=2 | All 2-mers distinct | {AG, GA, AT} | Source 1 AGAT example |
| M3 | FindUniqueKmers AAAAA k=3 | Homopolymer | empty (AAA Count = 3) | Source 2 (Count>1 ⇒ not unique) |
| M4 | FindKmersWithMinCount ACGTACGT k=4 min=2 | Recurrent filter | {(ACGT, 2)} | Source 3 + enumerated occurrences |
| M5 | FindKmersWithMinCount ACGTACGT k=4 min=1 | All distinct, ordered | 4 pairs; ACGT(2) first, then 3 Count-1 pairs | Source 3 ordering / Source 2 distinct |
| M6 | FindKmersWithMinCount Counts correct | Count value matches occurrences | ACGT → 2 | Source 1/3 overlapping count |
| M7 | FindUniqueKmers empty / k>L | Edge | empty | Source 1 (L−k+1≤0) |
| M8 | FindKmersWithMinCount empty / k>L | Edge | empty | Source 1 |
| M9 | FindUniqueKmers k≤0 | Invalid | ArgumentOutOfRangeException | Definition: k must be positive |
| M10 | FindKmersWithMinCount k≤0 | Invalid | ArgumentOutOfRangeException | Definition: k must be positive |
| M11 | INV-1 property: every returned unique k-mer has Count 1 | Cross-check vs CountKmers | all Count == 1 | Source 2 |
| M12 | INV-5: min=1 keys == distinct keys | Cross-check vs CountKmers | equal sets | Source 2/3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | FindUniqueKmers lower-case input | Case normalisation | same as upper-case | Evidence Assumption 2 |
| S2 | FindKmersWithMinCount order non-increasing | INV-4 | Counts sorted desc | Implementation contract |
| S3 | FindKmersWithMinCount min above max count | Threshold beyond data | empty | Source 3 (Count ≥ t) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | FindUniqueKmers k=1 monomers | Single-char k-mers unique iff appear once | per definition | Source 1 monomer example |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `KmerAnalyzerTests.cs` — two weak `FindKmersWithMinCount` tests (no messages, `GreaterThanOrEqualTo`, no edge coverage).
- `KmerAnalyzer_Find_Tests.cs` (KMER-FIND-001) — owns `FindUniqueKmers` deep tests under a different unit. Left untouched; this unit adds its own evidence-keyed `FindUniqueKmers` cases (ATCGATCAC / AGAT) per the registry assignment (see §7).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| FindKmersWithMinCount filters (legacy) | ⚠ Weak | No assertion messages; not evidence-keyed → rewrite as M4 |
| FindKmersWithMinCount ordering (legacy) | ⚠ Weak | `GreaterThanOrEqualTo` on first vs last only → rewrite as S2 |
| M1–M3, M11 (FindUniqueKmers) | ❌ Missing | New |
| M4–M6, M12, S2, S3 (FindKmersWithMinCount) | ❌ Missing | New |
| M7–M10, S1, C1 (edges/case/monomers) | ❌ Missing | New |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_FindUniqueAndMinCount_Tests.cs` — all KMER-UNIQUE-001 cases.
- **Remove:** the two weak `FindKmersWithMinCount` tests from `KmerAnalyzerTests.cs` (consolidated here).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| KmerAnalyzer_FindUniqueAndMinCount_Tests.cs | Canonical (this unit) | 15 |
| KmerAnalyzerTests.cs | Auxiliary (min-count tests removed) | unchanged otherwise |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Added ATCGATCAC k=3 → 5 unique, ATC excluded | ✅ Done |
| 2 | M2 | ❌ Missing | Added AGAT k=2 → {AG,GA,AT} | ✅ Done |
| 3 | M3 | ❌ Missing | Added AAAAA k=3 → empty | ✅ Done |
| 4 | M4 | ⚠ Weak→rewrite | ACGTACGT k=4 min=2 → {(ACGT,2)} with messages | ✅ Done |
| 5 | M5/M6 | ❌ Missing | min=1 all distinct, ordered, counts checked | ✅ Done |
| 6 | M7 | ❌ Missing | empty / k>L → empty | ✅ Done |
| 7 | M8 | ❌ Missing | empty / k>L → empty | ✅ Done |
| 8 | M9 | ❌ Missing | k≤0 throws | ✅ Done |
| 9 | M10 | ❌ Missing | k≤0 throws | ✅ Done |
| 10 | M11 | ❌ Missing | every unique k-mer has count 1 (INV-1) | ✅ Done |
| 11 | M12 | ❌ Missing | min=1 keys == distinct keys (INV-5) | ✅ Done |
| 12 | S1 | ❌ Missing | lower-case input matches | ✅ Done |
| 13 | S2 | ⚠ Weak→rewrite | counts non-increasing across full result (INV-4) | ✅ Done |
| 14 | S3 | ❌ Missing | threshold above max count → empty | ✅ Done |
| 15 | C1 | ❌ Missing | monomers k=1 → letters appearing once | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | FindUniqueKmers_AtcgatcacK3_ReturnsFiveUniqueExcludingRepeated |
| M2 | ✅ | FindUniqueKmers_AgatK2_ReturnsAllThreeTwoMers |
| M3 | ✅ | FindUniqueKmers_HomopolymerK3_ReturnsEmpty |
| M4 | ✅ | FindKmersWithMinCount_Acgtacgt_K4_Min2_ReturnsOnlyRepeated |
| M5/M6 | ✅ | FindKmersWithMinCount_Acgtacgt_K4_Min1_ReturnsAllDistinctOrderedDesc |
| M7 | ✅ | FindUniqueKmers_EmptyOrKExceedsLength_ReturnsEmpty |
| M8 | ✅ | FindKmersWithMinCount_EmptyOrKExceedsLength_ReturnsEmpty |
| M9 | ✅ | FindUniqueKmers_NonPositiveK_Throws |
| M10 | ✅ | FindKmersWithMinCount_NonPositiveK_Throws |
| M11 | ✅ | FindUniqueKmers_EveryReturnedKmer_HasCountOne |
| M12 | ✅ | FindKmersWithMinCount_Min1Keys_EqualDistinctKmerSet |
| S1 | ✅ | FindUniqueKmers_LowerCaseInput_MatchesUpperCase |
| S2 | ✅ | FindKmersWithMinCount_OutputCounts_AreNonIncreasing |
| S3 | ✅ | FindKmersWithMinCount_ThresholdAboveMaxCount_ReturnsEmpty |
| C1 | ✅ | FindUniqueKmers_MonomersK1_ReturnsLettersAppearingOnce |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | minCount ≤ 1 returns all distinct k-mers (Count ≥ t consistent extension) | M5, M12, INV-5 |
| 2 | Case normalisation (upper-cased) so case variants are the same k-mer | S1 |

---

## 7. Open Questions / Decisions

1. `FindUniqueKmers` is also referenced by the pre-existing KMER-FIND-001 TestSpec/test file. Per the Processing Registry method index, both `FindUniqueKmers` and `FindKmersWithMinCount` are assigned to KMER-UNIQUE-001; this unit's canonical file owns the deep `FindKmersWithMinCount` coverage and adds evidence-based `FindUniqueKmers` cases keyed to this unit's evidence (ATCGATCAC, AGAT). The weak legacy `FindKmersWithMinCount` tests in `KmerAnalyzerTests.cs` are consolidated into this unit's canonical file.
