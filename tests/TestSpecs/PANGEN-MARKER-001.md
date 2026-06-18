# Test Specification: PANGEN-MARKER-001

**Test Unit ID:** PANGEN-MARKER-001
**Area:** PanGenome
**Algorithm:** Phylogenetic Marker Selection (single-copy core genes ranked by parsimony-informative sites)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Ding, Baumdicker & Neher (2018), panX, *Nucleic Acids Research* 46(1):e5 | 1 | https://doi.org/10.1093/nar/gkx977 (https://pmc.ncbi.nlm.nih.gov/articles/PMC5758898/) | 2026-06-13 |
| 2 | Page et al. (2015), Roary, *Bioinformatics* 31(22):3691 | 1 | https://doi.org/10.1093/bioinformatics/btv421 (https://pmc.ncbi.nlm.nih.gov/articles/PMC4817141/) | 2026-06-13 |
| 3 | Roary documentation (Sanger Pathogens) | 3 | https://sanger-pathogens.github.io/Roary/ | 2026-06-13 |
| 4 | Zvelebil & Baum (2008), *Understanding Bioinformatics* (via Wikipedia "Informative site") | 4/1 | https://en.wikipedia.org/wiki/Informative_site | 2026-06-13 |

### 1.2 Key Evidence Points

1. Phylogenetic markers are **single-copy core gene clusters** — "gene clusters in which all strains are represented exactly once" — panX [1].
2. Core = present in **≥ 99%** of samples (default 0.99) — Roary [2].
3. Paralog-containing clusters are **filtered out** of the core gene alignment (markers must be single-copy orthologs) — Roary docs [3], Roary [2].
4. panX **extracts all variable positions** from the alignments of single-copy core genes (invariant columns carry no signal) — panX [1].
5. **Parsimony-informative site:** a column with "at least two different character states and each of those states occurs in at least two of the sequences" — Zvelebil 2008 [4].
6. Monomorphic and singleton columns are **not** parsimony-informative — Zvelebil 2008 [4].

### 1.3 Documented Corner Cases

- Not single-copy (genome with 0 or ≥ 2 genes in the cluster) → excluded [1][2][3].
- Below core threshold (absent from ≥ 1 genome) → excluded [2].
- Fully conserved single-copy core cluster (0 parsimony-informative sites) → excluded (no variable positions) [1].
- Monomorphic column → 0 PIS; singleton column → 0 PIS; two-state-each-≥2 column → 1 PIS [4].
- Member sequences of unequal length → no common alignment → PIS 0 (Assumption 1).

### 1.4 Known Failure Modes / Pitfalls

1. Counting singleton variants as informative (they are not) — Zvelebil 2008 [4].
2. Treating multi-copy / paralogous clusters as markers — panX/Roary forbid this [1][2][3].
3. Selecting conserved clusters with no variable positions — panX excludes them [1].

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CountParsimonyInformativeSites(alignedSequences)` | PanGenomeAnalyzer | Canonical | Column-wise PIS count per Zvelebil 2008 definition |
| `SelectPhylogeneticMarkers(genomes, coreClusters, totalGenomes, maxMarkers)` | PanGenomeAnalyzer | Canonical | Single-copy core selection + PIS ranking (panX/Roary) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | PIS count ≥ 0 and ≤ alignment length | Yes | Zvelebil 2008 [4] (column count bound) |
| INV-2 | A monomorphic or singleton column contributes 0 to the PIS count | Yes | Zvelebil 2008 [4] |
| INV-3 | A column with ≥ 2 states, each in ≥ 2 sequences, contributes exactly 1 | Yes | Zvelebil 2008 [4] |
| INV-4 | Every selected marker is single-copy core: present in all `totalGenomes` with exactly one gene per genome | Yes | panX [1]; Roary [2][3] |
| INV-5 | Every selected marker has ≥ 1 parsimony-informative site | Yes | panX "variable positions" [1] |
| INV-6 | Selected markers are ordered by descending PIS; result size ≤ `maxMarkers` | Yes | panX (most informative positions) [1]; informativeness ↔ PIS [4] |
| INV-7 | PIS is invariant to row (sequence) order and to a bijective relabeling of states | Yes | Column-content property of the definition [4] |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | PIS worked alignment | s1=`AAAAA`, s2=`AAACA`, s3=`AACCG`, s4=`ACCTG` | PIS = 2 (cols 3 & 5) | [4] + Evidence dataset |
| M2 | Monomorphic column | all rows `AAAA` (col) | 0 | [4] |
| M3 | Singleton column | col A,A,A,C | 0 | [4] |
| M4 | Minimal informative column | col A,A,C,C | 1 | [4] |
| M5 | Four-singleton column | col A,C,G,T | 0 (no state ≥ 2) | [4] |
| M6 | Excludes paralog cluster | core cluster where a genome has 2 genes | not selected | [1][2][3] |
| M7 | Excludes non-core cluster | cluster absent from one genome (below threshold) | not selected | [1][2] |
| M8 | Excludes conserved single-copy core (0 PIS) | identical members across all genomes | not selected | [1] |
| M9 | Selects single-copy core with PIS ≥ 1 | one informative marker among others | selected | [1][4] |
| M10 | Ranking by descending PIS | two markers, PIS 2 vs 1 | higher-PIS first | [1][4] |
| M11 | maxMarkers cap | 2 eligible markers, maxMarkers=1 | exactly 1 (the most informative) | [1] |
| M12 | Null/empty inputs | null clusters / null genomes / empty | empty result, no exception | corner cases |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Unequal-length members | cluster members of differing length | PIS 0 → not selected | Assumption 1 |
| S2 | Single sequence / fewer than 2 rows | one aligned sequence | PIS 0 | no state can occur in ≥ 2 rows |
| S3 | Empty alignment list / empty strings | no rows or zero-length rows | PIS 0 | INV-1 lower bound |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Row-order invariance | shuffle rows of M1 alignment | PIS unchanged (=2) | INV-7 |
| C2 | State-relabel invariance | swap A↔C in M1 alignment | PIS unchanged (=2) | INV-7 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzerTests.cs` — `#region Phylogenetic Marker Tests`: two tests (`SelectPhylogeneticMarkers_FiltersAndLimits`, `SelectPhylogeneticMarkers_PrefersLongerSequences`) exercising the old, unsourced identity-band + consensus-length heuristic.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1..M5 (PIS) | ❌ Missing | no PIS method existed |
| M6 (exclude paralog) | ❌ Missing | old code had no single-copy check |
| M7 (exclude non-core) | ❌ Missing | old code ignored core membership |
| M8 (exclude conserved) | ❌ Missing | old code had no PIS filter |
| M9 (select informative) | ❌ Missing | new contract |
| M10 (rank by PIS) | ⚠ Weak | old `SelectPhylogeneticMarkers_PrefersLongerSequences` ranked by length, not PIS — unsourced |
| M11 (maxMarkers cap) | ⚠ Weak | old `SelectPhylogeneticMarkers_FiltersAndLimits` used permissive `LessThanOrEqualTo` + unsourced band |
| M12, S1–S3, C1–C2 | ❌ Missing | not covered |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_SelectPhylogeneticMarkers_Tests.cs` — all PIS + marker-selection cases.
- **Remove:** `#region Phylogenetic Marker Tests` (both tests) from `PanGenomeAnalyzerTests.cs` — they assert an unsourced heuristic the implementation no longer uses.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PanGenomeAnalyzer_SelectPhylogeneticMarkers_Tests.cs` | Canonical (this unit) | 17 |
| `PanGenomeAnalyzerTests.cs` | Phylogenetic Marker region removed | 0 (for this unit) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ | Implemented worked-alignment PIS test | ✅ Done |
| 2 | M2 | ❌ | Implemented monomorphic-column test | ✅ Done |
| 3 | M3 | ❌ | Implemented singleton-column test | ✅ Done |
| 4 | M4 | ❌ | Implemented minimal-informative-column test | ✅ Done |
| 5 | M5 | ❌ | Implemented four-singleton-column test | ✅ Done |
| 6 | M6 | ❌ | Implemented paralog-exclusion test | ✅ Done |
| 7 | M7 | ❌ | Implemented non-core-exclusion test | ✅ Done |
| 8 | M8 | ❌ | Implemented conserved-exclusion test | ✅ Done |
| 9 | M9 | ❌ | Implemented select-informative test | ✅ Done |
| 10 | M10 | ⚠ | Rewrote ranking to descending PIS (exact) | ✅ Done |
| 11 | M11 | ⚠ | Rewrote cap test to exact count + identity | ✅ Done |
| 12 | M12 | ❌ | Implemented null/empty tests | ✅ Done |
| 13 | S1 | ❌ | Implemented unequal-length test | ✅ Done |
| 14 | S2 | ❌ | Implemented single-row test | ✅ Done |
| 15 | S3 | ❌ | Implemented empty-alignment test | ✅ Done |
| 16 | C1 | ❌ | Implemented row-order-invariance test | ✅ Done |
| 17 | C2 | ❌ | Implemented state-relabel-invariance test | ✅ Done |
| 18 | (cleanup) | 🔁 | Removed old heuristic region from PanGenomeAnalyzerTests.cs | ✅ Done |

**Total items:** 18
**✅ Done:** 18 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | `CountParsimonyInformativeSites_WorkedAlignment_ReturnsTwo` |
| M2 | ✅ | `CountParsimonyInformativeSites_MonomorphicColumn_ReturnsZero` |
| M3 | ✅ | `CountParsimonyInformativeSites_SingletonColumn_ReturnsZero` |
| M4 | ✅ | `CountParsimonyInformativeSites_MinimalInformativeColumn_ReturnsOne` |
| M5 | ✅ | `CountParsimonyInformativeSites_FourSingletonsColumn_ReturnsZero` |
| M6 | ✅ | `SelectPhylogeneticMarkers_ParalogCluster_Excluded` |
| M7 | ✅ | `SelectPhylogeneticMarkers_NonCoreCluster_Excluded` |
| M8 | ✅ | `SelectPhylogeneticMarkers_ConservedSingleCopyCore_Excluded` |
| M9 | ✅ | `SelectPhylogeneticMarkers_InformativeSingleCopyCore_Selected` |
| M10 | ✅ | `SelectPhylogeneticMarkers_OrdersByDescendingPis` |
| M11 | ✅ | `SelectPhylogeneticMarkers_MaxMarkers_CapsToMostInformative` |
| M12 | ✅ | `SelectPhylogeneticMarkers_NullOrEmpty_ReturnsEmpty` |
| S1 | ✅ | `SelectPhylogeneticMarkers_UnequalLengthMembers_NotSelected` |
| S2 | ✅ | `CountParsimonyInformativeSites_SingleSequence_ReturnsZero` |
| S3 | ✅ | `CountParsimonyInformativeSites_EmptyOrNull_ReturnsZero` |
| C1 | ✅ | `CountParsimonyInformativeSites_RowOrderPermuted_Unchanged` |
| C2 | ✅ | `CountParsimonyInformativeSites_StateRelabeled_Unchanged` |

✅ count = 17 in-scope planned cases (M1–M12, S1–S3, C1–C2); all ✅. No ❌/⚠ remain.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Per-cluster member alignment = equal-length member sequences (ungapped, position-wise); unequal lengths → PIS 0. Affects only how the alignment is obtained, not the PIS criterion. | `CountParsimonyInformativeSites`, `SelectPhylogeneticMarkers` (S1) |

---

## 7. Open Questions / Decisions

1. **Signature change (decided):** the canonical method now takes `(genomes, coreClusters, totalGenomes, maxMarkers)` so it can recover per-cluster member sequences for PIS, mirroring `CreateCoreGenomeAlignment`/`CreatePresenceAbsenceMatrix`. The old `(coreClusters, maxMarkers, minIdentity, maxIdentity)` identity-band signature used unsourced thresholds and is removed; the MCP wrapper is updated in lock-step (in-scope: direct consumer of the changed method).
