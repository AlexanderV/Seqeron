# Test Specification: KMER-STATS-001

**Test Unit ID:** KMER-STATS-001
**Area:** K-mer
**Algorithm:** K-mer Statistics (`KmerAnalyzer.AnalyzeKmers`)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wikipedia — K-mer (L−k+1; example tables AGAT, GTAGAGCTGT) | 4 | https://en.wikipedia.org/wiki/K-mer | 2026-06-14 |
| 2 | BioInfoLogics — k-mer counting part I (ATCGATCAC counts; distinct vs unique) | 4 | https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/ | 2026-06-14 |
| 3 | Manca et al. (2021) Spectral concepts in genome informational analysis (k-entropy E_k = −Σ p log₂ p, p=mult/(L−k+1)) | 1 | https://arxiv.org/abs/2106.15351 | 2026-06-14 |
| 4 | arXiv:2511.05300 Entropy–Rank Ratio (single-sequence k-mer Shannon entropy form) | 1 | https://arxiv.org/html/2511.05300 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Total number of k-mers in a length-L sequence is L−k+1 (overlapping windows) — Wikipedia K-mer; BioInfoLogics.
2. "Distinct" k-mers = each different k-mer counted once; the `UniqueKmers` field of `AnalyzeKmers` reports this distinct count — Wikipedia example tables (GTAGAGCTGT k=2 → 7 distinct); BioInfoLogics (ATCGATCAC k=3 → 6 distinct).
3. K-mer Shannon entropy E_k = −Σ p(α) log₂ p(α) with p(α) = mult(α)/(L−k+1) — Manca et al. (2021); corroborated by arXiv:2511.05300 as H_k(s) = −Σ p_i log₂ p_i.
4. Average k-mer multiplicity = total/distinct = (L−k+1)/distinct — derived from totals in the Wikipedia table.
5. Max/Min count are the extremes of the multiplicity distribution — Wikipedia example table (GTAGAGCTGT k=1 → max 4, min 1).

### 1.3 Documented Corner Cases

- k > L → L−k+1 ≤ 0 → no k-mers → all-zero statistics (Wikipedia L−k+1 formula).
- Single distinct k-mer (homopolymer) → entropy = −1·log₂1 = 0; max==min==total (k-entropy definition).
- All windows distinct → entropy = log₂(L−k+1); max==min==1 (k-entropy definition).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing "distinct" with "unique (count==1)": the `UniqueKmers` field is the distinct count, NOT the count-1 count — BioInfoLogics distinct/unique distinction.
2. Using log base e or 10 instead of base 2 changes entropy units — Manca et al.; arXiv:2511.05300 (bits, log₂).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `AnalyzeKmers(string sequence, int k)` | KmerAnalyzer | **Canonical** | Returns `KmerStatistics(TotalKmers, UniqueKmers, MaxCount, MinCount, AverageCount, Entropy)` |

<!-- Type values: -->
<!-- **Canonical** — deep evidence-based testing -->
<!-- **Delegate** — smoke verification only (1–2 tests proving delegation) -->
<!-- **Internal** — tested indirectly via canonical methods -->

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | TotalKmers = L − k + 1 for L ≥ k (number of overlapping windows) | Yes | Wikipedia K-mer; BioInfoLogics |
| INV-2 | TotalKmers = sum over all distinct k-mers of their counts | Yes | Definition of multiplicity (Manca p(α)=mult/(L−k+1)) |
| INV-3 | UniqueKmers = number of distinct k-mers (CountKmers key count) | Yes | Wikipedia example tables; BioInfoLogics |
| INV-4 | MinCount ≤ AverageCount ≤ MaxCount; AverageCount = TotalKmers/UniqueKmers | Yes | Arithmetic / derived from totals |
| INV-5 | 0 ≤ Entropy ≤ log₂(UniqueKmers); Entropy = 0 iff one distinct k-mer; = log₂(distinct) iff all counts equal | Yes | Manca et al. k-entropy; Shannon bounds |
| INV-6 | k > L or empty sequence ⇒ all fields = 0 | Yes | L−k+1 ≤ 0 (Wikipedia formula) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | GTAGAGCTGT k=1 | Full statistics for monomer counts | Total=10, Unique=4, Max=4, Min=1, Avg=2.5, Entropy=1.846439344671 | Wikipedia table (G4 T3 A2 C1) |
| M2 | GTAGAGCTGT k=2 | 2-mer statistics with two doubled k-mers | Total=9, Unique=7, Max=2, Min=1, Avg=1.29, Entropy=2.725480556998 | Wikipedia table (GT,AG ×2) |
| M3 | GTAGAGCTGT k=3 | All 8 windows distinct | Total=8, Unique=8, Max=1, Min=1, Avg=1.0, Entropy=3.0 (=log₂8) | Wikipedia table (8 distinct) |
| M4 | ATCGATCAC k=3 | Distinct=6 with one doubled k-mer | Total=7, Unique=6, Max=2, Min=1, Avg=1.17, Entropy=2.521640636343 | BioInfoLogics table (ATC=2) |
| M5 | AGAT k=2 | All distinct, uniform | Total=3, Unique=3, Max=1, Min=1, Avg=1.0, Entropy=log₂3=1.584962500721 | Wikipedia AGAT example |
| M6 | INV-1/INV-2 cross-check | TotalKmers == L−k+1 AND == sum of CountKmers values | both hold for a longer seq | INV-1, INV-2 |
| M7 | INV-3 cross-check | UniqueKmers == CountKmers key count | holds | INV-3 |
| M8 | Empty sequence | k=3 over "" | all fields 0 | INV-6 (L−k+1≤0) |
| M9 | k > length | "ACG" k=5 | all fields 0 | INV-6 |
| M10 | k ≤ 0 | "ACGT" k=0 | throws ArgumentOutOfRangeException | KmerAnalyzer k>0 contract; CountKmers validation |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | AAAA k=2 homopolymer | single distinct k-mer | Total=3, Unique=1, Max=3, Min=3, Avg=3.0, Entropy=0 | INV-5 lower bound |
| S2 | INV-4 ordering | Min ≤ Avg ≤ Max and Avg=Total/Unique | holds for GTAGAGCTGT k=2 | derived |
| S3 | INV-5 entropy upper bound | Entropy ≤ log₂(Unique) for several inputs | holds | Shannon bound |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Case-insensitivity | "gtagagctgt" k=1 equals upper-case stats | identical KmerStatistics | implementation upper-cases internally |
| C2 | null sequence | null k=3 | all fields 0 | CountKmers treats null as empty |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `AnalyzeKmers` / `KmerStatistics`. No existing canonical test file for KMER-STATS-001. Sibling K-mer fixtures exist (`KmerAnalyzer_FindUniqueAndMinCount_Tests.cs`, `KmerAnalyzer_CountKmersBothStrands_Tests.cs`, etc.) but none cover `AnalyzeKmers`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 GTAGAGCTGT k=1 | ❌ Missing | new |
| M2 GTAGAGCTGT k=2 | ❌ Missing | new |
| M3 GTAGAGCTGT k=3 | ❌ Missing | new |
| M4 ATCGATCAC k=3 | ❌ Missing | new |
| M5 AGAT k=2 | ❌ Missing | new |
| M6 INV-1/INV-2 | ❌ Missing | new |
| M7 INV-3 | ❌ Missing | new |
| M8 empty | ❌ Missing | new |
| M9 k>L | ❌ Missing | new |
| M10 k≤0 throws | ❌ Missing | new |
| S1 homopolymer | ❌ Missing | new |
| S2 INV-4 ordering | ❌ Missing | new |
| S3 entropy upper bound | ❌ Missing | new |
| C1 case-insensitivity | ❌ Missing | new |
| C2 null | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_AnalyzeKmers_Tests.cs` — all KMER-STATS-001 tests.
- **Remove:** none (no prior tests for this unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `KmerAnalyzer_AnalyzeKmers_Tests.cs` | Canonical fixture for KMER-STATS-001 | 15 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | S1 | ❌ Missing | Implemented | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented | ✅ Done |
| 13 | S3 | ❌ Missing | Implemented | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented | ✅ Done |
| 15 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | AnalyzeKmers_GtagagctgtK1_... |
| M2 | ✅ Covered | AnalyzeKmers_GtagagctgtK2_... |
| M3 | ✅ Covered | AnalyzeKmers_GtagagctgtK3_AllDistinct_... |
| M4 | ✅ Covered | AnalyzeKmers_AtcgatcacK3_... |
| M5 | ✅ Covered | AnalyzeKmers_AgatK2_... |
| M6 | ✅ Covered | AnalyzeKmers_TotalKmers_EqualsLMinusKPlus1AndSumOfCounts |
| M7 | ✅ Covered | AnalyzeKmers_UniqueKmers_EqualsDistinctCount |
| M8 | ✅ Covered | AnalyzeKmers_EmptySequence_ReturnsAllZero |
| M9 | ✅ Covered | AnalyzeKmers_KExceedsLength_ReturnsAllZero |
| M10 | ✅ Covered | AnalyzeKmers_NonPositiveK_Throws |
| S1 | ✅ Covered | AnalyzeKmers_HomopolymerK2_EntropyZero |
| S2 | ✅ Covered | AnalyzeKmers_MinAvgMax_Ordering |
| S3 | ✅ Covered | AnalyzeKmers_Entropy_WithinLog2DistinctBound |
| C1 | ✅ Covered | AnalyzeKmers_LowerCaseInput_MatchesUpperCase |
| C2 | ✅ Covered | AnalyzeKmers_NullSequence_ReturnsAllZero |

**Total in-scope cases:** 15 | **✅:** 15

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | AverageCount rounded to 2 decimals (presentation only; exact ratio also verified) | M2, M4 expected values |
| 2 | Entropy reported unrounded in bits (log base 2) | M1–M5, S1, S3 |

---

## 7. Open Questions / Decisions

1. The `UniqueKmers` field name denotes the **distinct** k-mer count (not the count-1 "unique" set of KMER-UNIQUE-001). This naming is retained for API stability; documented in the algorithm doc to avoid the distinct/unique confusion (failure mode 1.4.1).
