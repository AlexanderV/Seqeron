# Test Specification: MOTIF-DISCOVER-001

**Test Unit ID:** MOTIF-DISCOVER-001
**Area:** Matching
**Algorithm:** Motif Discovery via Overrepresented k-mers (observed/expected enrichment)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Compeau & Pevzner, *Bioinformatics Algorithms*, expected k-mer occurrences | 1 | https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string | 2026-06-14 |
| 2 | monaLisa `getKmerFreq` (O/E ratio, log2 enrichment) | 3 | https://fmicompbio.github.io/monaLisa/reference/getKmerFreq.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. Expected occurrences of a specific k-mer in a length-N string under the i.i.d. uniform (each base p=1/4) background is `E = (N − k + 1) / 4^k` — Source 1.
2. Overrepresentation is the observed/expected (O/E) ratio; value > 1 means overrepresented — Sources 1, 2.
3. `N − k + 1` is the number of length-k windows; `4^k` is the number of distinct DNA k-mers — Source 1.

### 1.3 Documented Corner Cases

- `k > N`: zero length-k windows, no k-mer can be counted (no motifs returned) — Source 1.
- The self-overlap caveat affects only the *probability* statistic, not the deterministic observed count or the O/E denominator — Source 1.

### 1.4 Known Failure Modes / Pitfalls

1. Using a floor/clamp on the expected count (e.g. `max(E, 0.1)`) is not part of the published statistic and distorts the O/E ratio — Source 1 (formula has no clamp).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DiscoverMotifs(DnaSequence sequence, int k = 6, int minCount = 2)` | MotifFinder | Canonical | Returns `DiscoveredMotif` records (Sequence, Count, Positions, Enrichment) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Count equals the number of occurrences of the k-mer; Positions are the 0-based window starts of those occurrences | Yes | Source 1 (window enumeration) |
| INV-2 | Enrichment = Count / E where E = (N − k + 1) / 4^k | Yes | Source 1 |
| INV-3 | Every returned motif has Count ≥ minCount | Yes | Method contract |
| INV-4 | Enrichment > 0 for every returned motif (E > 0 when any k-mer exists) | Yes | Source 1 (N−k+1 ≥ 1) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Repeated k-mer found | "ATGCATGCATGC", k=4, minCount=2 | result contains "ATGC" with Count=3 | Source 1 (window count) |
| M2 | Positions reported | "ATGCATGCATGC", k=4, minCount=2 | "ATGC" Positions = {0,4,8} exactly | INV-1 |
| M3 | Exact O/E enrichment (tandem) | "ATGCATGCATGC", k=4 | "ATGC" Enrichment = 768/9 = 85.3333… (3 / (9/256)) | Source 1 formula |
| M4 | Exact O/E enrichment (homopolymer) | "AAAAAAAAAA", k=3 | "AAA" Count=8, Enrichment = 64.0 (8 / (8/64)) | Source 1 formula |
| M5 | minCount filter | "ATGCAAAA", k=4, minCount=2 | every returned motif has Count ≥ 2 | INV-3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null sequence | sequence = null | ArgumentNullException | Validation contract |
| S2 | k < 1 | k = 0 | ArgumentOutOfRangeException | Validation contract |
| S3 | k > N | "AAA", k=5 | empty result | Corner case (no windows) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | No floor on expected count | very long sequence so E > 0.1, plus a long-k case where the old `max(E,0.1)` clamp would have changed the value | Enrichment = Count / ((N−k+1)/4^k) exactly | Guards against re-introducing the clamp defect |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinderTests.cs` — `#region Motif Discovery Tests` contains `DiscoverMotifs_FindsRepeatedKmer`, `DiscoverMotifs_ReturnsPositions`, `DiscoverMotifs_CalculatesEnrichment`, `DiscoverMotifs_FiltersByMinCount`, plus `DiscoverMotifs_NullSequence_ThrowsException`, `DiscoverMotifs_ZeroK_ThrowsException`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (repeated k-mer found) | ⚠ Weak | Existing `DiscoverMotifs_FindsRepeatedKmer` uses `Any(...)`, no Count assertion |
| M2 (positions) | ⚠ Weak | `DiscoverMotifs_ReturnsPositions` uses `Does.Contain` (no exact set / boundary) |
| M3 (exact O/E tandem) | ❌ Missing | No exact enrichment value asserted anywhere |
| M4 (exact O/E homopolymer) | ⚠ Weak | `DiscoverMotifs_CalculatesEnrichment` only asserts `> 1` (permissive) |
| M5 (minCount filter) | ✅ Covered | `DiscoverMotifs_FiltersByMinCount` asserts Count ≥ 2 (kept; mirrored in canonical file) |
| S1 (null) | ✅ Covered | `DiscoverMotifs_NullSequence_ThrowsException` (mirrored) |
| S2 (k < 1) | ✅ Covered | `DiscoverMotifs_ZeroK_ThrowsException` (mirrored) |
| S3 (k > N) | ❌ Missing | Not tested |
| C1 (no floor) | ❌ Missing | Not tested |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_DiscoverMotifs_Tests.cs` — all MUST/SHOULD/COULD cases with exact evidence values.
- **Remove:** the six weak/duplicate `DiscoverMotifs_*` tests from `MotifFinderTests.cs` (their coverage is superseded by the canonical file with exact values).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MotifFinder_DiscoverMotifs_Tests.cs` | Canonical | 9 |
| `MotifFinderTests.cs` (`Motif Discovery` region) | removed | 0 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | Rewrote with exact Count=3 | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote with exact Positions {0,4,8} | ✅ Done |
| 3 | M3 | ❌ Missing | Added exact Enrichment 768/9 | ✅ Done |
| 4 | M4 | ⚠ Weak | Rewrote with exact Enrichment 64.0 + Count 8 | ✅ Done |
| 5 | M5 | ✅ Covered | Mirrored into canonical file | ✅ Done |
| 6 | S1 | ✅ Covered | Mirrored | ✅ Done |
| 7 | S2 | ✅ Covered | Mirrored | ✅ Done |
| 8 | S3 | ❌ Missing | Added k>N empty-result test | ✅ Done |
| 9 | C1 | ❌ Missing | Added no-floor exact-ratio test | ✅ Done |

**Total items:** 9
**✅ Done:** 9 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact Count=3 in canonical file |
| M2 | ✅ | Exact Positions {0,4,8} |
| M3 | ✅ | Enrichment = 768/9 within 1e-10 |
| M4 | ✅ | Enrichment = 64.0, Count=8 |
| M5 | ✅ | minCount filter asserted |
| S1 | ✅ | ArgumentNullException |
| S2 | ✅ | ArgumentOutOfRangeException |
| S3 | ✅ | k>N empty |
| C1 | ✅ | No-floor exact O/E |

**Total in-scope cases:** 9 | **✅:** 9

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | `minCount` is a presentation filter, not part of the published statistic (not correctness-affecting per record) | M5 |

---

## 7. Open Questions / Decisions

1. The checklist signature reads `DiscoverMotifs(sequences, k)`, but the registered Type is "Overrepresented k-mers" and the implemented canonical method operates on a single `DnaSequence` (per-sequence k-mer overrepresentation). The single-sequence overrepresentation method is the canonical one; the cross-sequence variant is the separate unit MOTIF-SHARED-001 (`FindSharedMotifs`). External evidence (Compeau & Pevzner expected-count formula) defines the single-sequence statistic, so testing targets `DiscoverMotifs(DnaSequence, k, minCount)`.
