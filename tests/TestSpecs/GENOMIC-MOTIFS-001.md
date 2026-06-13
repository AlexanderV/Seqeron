# Test Specification: GENOMIC-MOTIFS-001

**Test Unit ID:** GENOMIC-MOTIFS-001
**Area:** Analysis
**Algorithm:** Known Motif Search (multi-pattern exact substring matching)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Gusfield (1997) exact-matching definition via Tufts COMP 150GEN | 1 | https://www.cs.tufts.edu/comp/150GEN/classpages/exact.html | 2026-06-13 |
| 2 | Biopython `Bio.Seq` (`search`, `count_overlap`) reference impl | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py | 2026-06-13 |
| 3 | Wikipedia "Restriction site" (EcoRI = GAATTC) | 4 | https://en.wikipedia.org/wiki/Restriction_site | 2026-06-13 |

### 1.2 Key Evidence Points

1. Exact matching = "find all occurrences of a pattern string P ... in a text string T"; the answer is the set of all start positions. — Source 1.
2. Occurrences include overlaps: P=aaa, T=aaaaa → 3 occurrences (positions 0,1,2). — Source 1.
3. Multi-motif search: Biopython `Seq.search` accepts a set of motifs and yields per-motif `(index, substring)` hits. — Source 2.
4. Overlapping count is correct: `Seq("AAAA").count_overlap("AA")` → 3 (vs non-overlap 2). — Source 2.
5. DNA is processed upper-cased. — Source 2.
6. EcoRI motif is `GAATTC` (length 6), a real biological motif. — Source 3.

### 1.3 Documented Corner Cases

- Overlapping occurrences must all be reported (Source 1).
- Absent pattern → empty occurrence set (Source 1).
- Multiple query motifs each map to their own position set (Source 2).
- Case normalization to uppercase (Source 2).

### 1.4 Known Failure Modes / Pitfalls

1. Non-overlapping (greedy skip-ahead) search under-reports occurrences of self-overlapping motifs — defect vs Source 1.
2. Returning positions in nondeterministic (DFS) order from the suffix tree breaks a stable public contract — positions must be sorted ascending.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindKnownMotifs(DnaSequence, IEnumerable<string>)` | GenomicAnalyzer | Canonical | Multi-motif exact search; returns dict motif → sorted 0-based positions |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported position p for motif m satisfies `T[p..p+|m|-1] == m` (0-based) | Yes | Source 1 (exact-matching definition) |
| INV-2 | All occurrences are reported, including overlapping ones | Yes | Source 1 (P=aaa/T=aaaaa→3); Source 2 (count_overlap) |
| INV-3 | Positions for each motif are sorted strictly ascending and distinct | Yes | Source 1 (set of positions) + stable-contract requirement |
| INV-4 | A motif with zero occurrences is absent from the result dictionary | Yes | Source 1 (empty set); Source 2 (per-motif hits only when found) |
| INV-5 | Result keys are upper-cased motif strings | Yes | Source 2 (uppercase processing) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Overlapping homopolymer | T=`AAAAA`, motifs={`AAA`} | `{ "AAA": [0,1,2] }` | Source 1 worked example (3 overlapping occurrences) |
| M2 | EcoRI biological motif | T=`GAATTCAAAGAATTC`, motifs={`GAATTC`} | `{ "GAATTC": [0,9] }` | Source 3 (EcoRI=GAATTC) + Source 1 (positions) |
| M3 | Multi-motif set, some absent | T=`ACGTACGTAA`, motifs={`ACGT`,`AA`,`TTT`} | `{ "ACGT":[0,4], "AA":[8] }` (TTT omitted) | Source 2 (per-motif hits, absent omitted) |
| M4 | Positions sorted ascending | T=`ACGTACGTAA`, motif `ACGT` | list is `[0,4]` strictly ascending | Source 1 (set) + INV-3 |
| M5 | Absent motif omitted | T=`ACGT`, motifs={`TTTT`} | empty dictionary | Source 1 (empty set); INV-4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty motif set | motifs={} | empty dictionary | degenerate input |
| S2 | Lower/mixed-case motif normalized | T=`GAATTC`, motifs={`gaattc`} | `{ "GAATTC": [0] }` | Source 2 uppercase; key upper-cased (INV-5) |
| S3 | Empty/whitespace motif skipped | motifs={`""`} | empty dictionary | ASSUMPTION: empty-motif policy (Evidence Assumptions §1) |
| S4 | Single-char motif counts all positions | T=`AAAAA`, motifs={`A`} | `{ "A":[0,1,2,3,4] }` | overlap rule boundary (Source 1) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Duplicate motif in input | motifs={`AA`,`aa`} normalize to same key | one `AA` entry, positions correct | dedup via dict key |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `FindKnownMotifs`. No existing test file/cases for this method. Canonical file `GenomicAnalyzer_FindKnownMotifs_Tests.cs` will be created.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindKnownMotifs_Tests.cs` — all cases above.
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| GenomicAnalyzer_FindKnownMotifs_Tests.cs | Canonical | 10 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented overlapping homopolymer test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented EcoRI test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented multi-motif set test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented sorted-ascending test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented absent-motif test | ✅ Done |
| 6 | S1 | ❌ Missing | Implemented empty-set test | ✅ Done |
| 7 | S2 | ❌ Missing | Implemented case-normalization test | ✅ Done |
| 8 | S3 | ❌ Missing | Implemented empty-motif-skipped test | ✅ Done |
| 9 | S4 | ❌ Missing | Implemented single-char overlap test | ✅ Done |
| 10 | C1 | ❌ Missing | Implemented duplicate-motif dedup test | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | FindKnownMotifs_OverlappingMotif_ReportsAllStarts |
| M2 | ✅ Covered | FindKnownMotifs_EcoRIMotif_FindsAllSites |
| M3 | ✅ Covered | FindKnownMotifs_MotifSet_ReturnsPerMotifPositions |
| M4 | ✅ Covered | FindKnownMotifs_Positions_AreSortedAscending |
| M5 | ✅ Covered | FindKnownMotifs_AbsentMotif_OmittedFromResult |
| S1 | ✅ Covered | FindKnownMotifs_EmptyMotifSet_ReturnsEmpty |
| S2 | ✅ Covered | FindKnownMotifs_LowercaseMotif_NormalizedAndMatched |
| S3 | ✅ Covered | FindKnownMotifs_EmptyMotif_Skipped |
| S4 | ✅ Covered | FindKnownMotifs_SingleCharMotif_ReportsEveryPosition |
| C1 | ✅ Covered | FindKnownMotifs_DuplicateMotifKeys_Deduplicated |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Empty/whitespace motif contributes no result entry (repository policy, mirrors `FindMotif`) | S3 |
| 2 | Result key is the upper-cased motif (Biopython case handling) | S2, C1, INV-5 |

---

## 7. Open Questions / Decisions

1. Existing `FindKnownMotifs` returned suffix-tree DFS-order (unsorted) positions and did not skip empty motifs — corrected to sorted ascending + empty-motif skip to satisfy INV-3 and the empty-motif policy. Decision recorded in algorithm doc §5.2.
