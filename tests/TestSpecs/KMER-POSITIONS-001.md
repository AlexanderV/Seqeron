# Test Specification: KMER-POSITIONS-001

**Test Unit ID:** KMER-POSITIONS-001
**Area:** K-mer
**Algorithm:** K-mer Positions (find all start positions of a k-mer in a sequence)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Rosalind BA1D — Find All Occurrences of a Pattern in a String | 4 | https://rosalind.info/problems/ba1d/ | 2026-06-14 |
| 2 | Wikipedia — k-mer | 4 | https://en.wikipedia.org/wiki/K-mer | 2026-06-14 |
| 3 | Compeau & Pevzner, Bioinformatics Algorithms (Pattern Matching) | 1 | https://gerdos.web.elte.hu/edu/bioinformatics_algorithms/week1.pdf | 2026-06-14 |

### 1.2 Key Evidence Points

1. Output is all starting positions where the pattern appears as a substring, **0-based** — Rosalind BA1D.
2. Overlapping occurrences are all reported (self-overlap counted) — Rosalind BA1D.
3. A k-mer is a substring of length k; there are **L − k + 1** start positions in a length-L sequence — Wikipedia k-mer.
4. Worked example: `ATAT` in `GATATATGCATATACTT` → `1 3 9` — Rosalind BA1D.

### 1.3 Documented Corner Cases

- Overlapping self-occurrences (e.g. `AA` in `AAAA` → 0,1,2) — Rosalind BA1D.
- Pattern longer than text → no valid start positions (L − k + 1 ≤ 0) → empty — Wikipedia count formula.
- Pattern absent → empty; pattern equals whole sequence → one position (0).

### 1.4 Known Failure Modes / Pitfalls

1. Reporting only non-overlapping occurrences (incorrect) — must report every overlapping start — Rosalind BA1D.
2. Off-by-one / 1-based indexing — output must be 0-based — Rosalind BA1D.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindKmerPositions(string sequence, string kmer)` | KmerAnalyzer | Canonical | Returns 0-based start positions, ascending, overlapping. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned position p satisfies sequence[p..p+|kmer|] == kmer (case-folded). | Yes | Rosalind BA1D |
| INV-2 | Positions are 0-based and strictly ascending. | Yes | Rosalind BA1D |
| INV-3 | Result count equals the overlapping occurrence count of kmer in sequence (self-overlaps included). | Yes | Rosalind BA1D |
| INV-4 | All positions lie in [0, |sequence| − |kmer|]; empty when |kmer| > |sequence|. | Yes | Wikipedia k-mer (L−k+1) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Rosalind BA1D sample | `ATAT` in `GATATATGCATATACTT` | `[1, 3, 9]` | Rosalind BA1D |
| M2 | Overlapping self-occurrence | `AA` in `AAAA` | `[0, 1, 2]` | BA1D overlapping rule |
| M3 | "ana" in "banana" | classic overlapping pattern | `[1, 3]` | SuffixTree doc / BA1D |
| M4 | AGAT 2-mers | `AG`,`GA`,`AT` each occur once | AG→`[0]`, GA→`[1]`, AT→`[2]` | Wikipedia k-mer |
| M5 | Pattern absent | `GG` in `ATATAT` | `[]` | BA1D (only matching starts) |
| M6 | Ascending order | overlapping matches returned in increasing index order | strictly ascending | BA1D / INV-2 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Pattern longer than text | `ACGT` in `AC` | `[]` | L−k+1 ≤ 0 |
| S2 | Pattern equals whole sequence | `ACGT` in `ACGT` | `[0]` | single occurrence |
| S3 | Case-insensitive match | `atat` in `GATATATGCATATACTT` | `[1, 3, 9]` | case-folding assumption |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | null/empty sequence | null or "" sequence | `[]` | null/empty assumption |
| C2 | null/empty kmer | null or "" kmer | `[]` | null/empty assumption |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior canonical test file for `FindKmerPositions`. Existing `KmerAnalyzer_*_Tests.cs` cover sibling methods only; `KmerAnalyzer_Find_Tests.cs` covers FindMostFrequentKmers/FindUniqueKmers/FindClumps, not `FindKmerPositions`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New unit |
| M2 | ❌ Missing | New unit |
| M3 | ❌ Missing | New unit |
| M4 | ❌ Missing | New unit |
| M5 | ❌ Missing | New unit |
| M6 | ❌ Missing | New unit |
| S1 | ❌ Missing | New unit |
| S2 | ❌ Missing | New unit |
| S3 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |
| C2 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_FindKmerPositions_Tests.cs` — all cases for this unit.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `KmerAnalyzer_FindKmerPositions_Tests.cs` | Canonical | 11 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented | ✅ Done |
| 10 | C1 | ❌ Missing | Implemented | ✅ Done |
| 11 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | FindKmerPositions_RosalindBA1D_ReturnsExpectedPositions |
| M2 | ✅ Covered | FindKmerPositions_OverlappingSelfOccurrence_ReturnsAllStarts |
| M3 | ✅ Covered | FindKmerPositions_AnaInBanana_ReturnsOverlappingStarts |
| M4 | ✅ Covered | FindKmerPositions_AgatTwoMers_ReturnsEachStart |
| M5 | ✅ Covered | FindKmerPositions_PatternAbsent_ReturnsEmpty |
| M6 | ✅ Covered | FindKmerPositions_OverlappingMatches_ReturnedInAscendingOrder |
| S1 | ✅ Covered | FindKmerPositions_PatternLongerThanText_ReturnsEmpty |
| S2 | ✅ Covered | FindKmerPositions_PatternEqualsSequence_ReturnsZeroOnly |
| S3 | ✅ Covered | FindKmerPositions_LowercaseKmer_MatchesCaseInsensitively |
| C1 | ✅ Covered | FindKmerPositions_NullOrEmptySequence_ReturnsEmpty |
| C2 | ✅ Covered | FindKmerPositions_NullOrEmptyKmer_ReturnsEmpty |

**✅ count: 11 = total in-scope cases.**

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | 0-based indexing | INV-2, M1–M6 |
| 2 | Case-insensitive matching | S3 |
| 3 | Null/empty input → empty result | C1, C2 |

---

## 7. Open Questions / Decisions

1. **Search reuse decision:** SuffixTree `FindAllOccurrences` was evaluated. It counts overlapping occurrences correctly but returns positions in unordered leaf-collection order and requires O(n) construction per text; for a single k-mer query against one text the naive O(n·m) scan is simpler and yields ascending order directly. Naive scan retained. Recorded in the algorithm doc §4.3 / §5.2.
