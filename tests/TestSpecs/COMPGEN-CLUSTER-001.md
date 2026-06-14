# Test Specification: COMPGEN-CLUSTER-001

**Test Unit ID:** COMPGEN-CLUSTER-001
**Area:** Comparative
**Algorithm:** Conserved Gene Clusters (common intervals of permutations)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Bui-Xuan, Habib, Paul (2013). MinMax-Profiles. arXiv:1304.5140 | 1 | https://arxiv.org/abs/1304.5140 | 2026-06-14 |
| 2 | Uno, Yagiura (2000). Fast Algorithms to Enumerate All Common Intervals of Two Permutations. Algorithmica 26(2):290–309 | 1 | https://doi.org/10.1007/s004539910014 | 2026-06-14 |
| 3 | Didier, Schmidt, Stoye, Tsur (2013). Extending Common Intervals Searching from Permutations to Sequences. arXiv:1310.4290 | 1 | https://arxiv.org/abs/1310.4290 | 2026-06-14 |
| 4 | Heber, Stoye (2001). Finding All Common Intervals of k Permutations. CPM 2001, LNCS 2089:207–218 | 1 | https://doi.org/10.1007/3-540-48194-X_19 | 2026-06-14 |

### 1.2 Key Evidence Points

1. A **common interval** of a family of permutations is "a set of integers that is an interval of each Pk" — i.e. a set of ortholog-group labels contiguous in *every* genome. — Source 1, Definition 1 (citing Source 2).
2. An interval `[i,j]` is "defined only for 1 ≤ i < j ≤ n", so common intervals have size ≥ 2; the whole set is always a (trivial) common interval. — Source 1, §2.
3. Worked Example 1: P1=Id7, P2=(7 2 1 3 6 4 5) → common intervals {1,2}, {1,2,3}, {3,4,5,6}, {4,5}, {4,5,6}, {1,…,6}, {1,…,7}. — Source 1, Example 1 (independently recomputed by brute force).
4. For sequences with repeated labels (paralogs), a set is a common interval iff *some* contiguous window in each genome has exactly that label set. — Source 3, Definition 1 + Example 1.
5. The model generalises to k ≥ 2 permutations (all conserved in *all* genomes). — Source 4.

### 1.3 Documented Corner Cases

- **Trivial intervals:** whole set always common; singletons excluded by `i < j` (Source 1, §2).
- **K < 2 genomes:** common interval is a family notion; undefined / vacuous for a single genome (Source 1, §2 family definition).
- **Repeated labels:** any matching contiguous window counts as a location (Source 3, Example 1).

### 1.4 Known Failure Modes / Pitfalls

1. Counting a set as a cluster when it is contiguous in only *one* genome — the set must be contiguous in **all** genomes. — Source 1, Definition 1.
2. Allowing foreign genes inside the window — a common interval is the set of **all** elements in the window, so no non-member group may sit inside. — Source 1, §2 (interval = Set of window).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindConservedClusters(genomes, orthologGroups, minClusterSize, maxGap)` | ComparativeGenomics | **Canonical** | Common-interval cluster detection over K genomes |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Every returned cluster is a set of ortholog-group labels that is contiguous (forms a window) in **every** input genome | Yes | Source 1, Definition 1 |
| INV-02 | Every returned cluster has size ≥ `minClusterSize` (and ≥ 2) | Yes | Source 1, §2 (i < j; size threshold) |
| INV-03 | A set contiguous in some but not all genomes is NOT returned | Yes | Source 1, Definition 1 |
| INV-04 | With fewer than 2 genomes the result is empty | Yes | Source 1, §2 (family definition) |
| INV-05 | Result is deterministic and independent of cluster discovery order (set semantics) | Yes | Source 2 (well-defined enumeration) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | GoldenVector_AllCommonIntervals | P1=Id7, P2=(7 2 1 3 6 4 5), minClusterSize=2 | Returned label-sets exactly = {1,2},{1,2,3},{3,4,5,6},{4,5},{4,5,6},{1..6},{1..7} | Source 1, Example 1 |
| M2 | SetContiguousInOneGenomeOnly_NotReturned | {2,3}: contiguous in P1, split in P2 | {2,3} is NOT among returned clusters | Source 1, Example 1 (positions of 2,3 in P2 = 2,4) |
| M3 | RepeatedLabels_WindowMatch | Genomes with a duplicated group label where a window has the same set in both | The set is reported despite duplicates | Source 3, Example 1 |
| M4 | MinClusterSize_FiltersSmall | Golden vector with minClusterSize=4 | Only sets of size ≥ 4 returned ({3,4,5,6},{1..6},{1..7}) | Source 1, §2 + Example 1 |
| M5 | AllMembersContiguousInAllGenomes | Cluster window contains no foreign group between members in all genomes | Returned cluster's set = exactly the window's group set | Source 1, §2 (interval = Set of window) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | SingleGenome_ReturnsEmpty | One genome only | Empty | INV-04 |
| S2 | IdenticalOrderThreeGenomes | Three genomes, identical group order 1..5, minClusterSize=3 | All windows of size ≥ 3 conserved: {1,2,3},{2,3,4},{3,4,5},{1,2,3,4},{2,3,4,5},{1,2,3,4,5} | Identity vs identity, all intervals common (Source 1) |
| S3 | ThreeGenomes_OnlyAllConservedReported | Set conserved in genomes 1,2 but split in genome 3 | Not returned | INV-01/INV-03, Source 4 (k-permutation) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism_StableAcrossRuns | Same inputs twice | Identical cluster set both runs | INV-05 |
| C2 | NoConservedCluster_ReturnsEmpty | Genomes with no shared contiguous set ≥ minClusterSize | Empty | Boundary |
| C3 | NullArguments_Throw | null genomes / null orthologGroups | ArgumentNullException | Defensive contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomicsTests.cs` — `#region FindConservedClusters Tests` (lines 103–147): two tests, `FindConservedClusters_ConservedCluster_ReturnsCluster` and `FindConservedClusters_SingleGenome_ReturnsEmpty`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| existing: ConservedCluster_ReturnsCluster | ⚠ Weak | Permissive `GreaterThanOrEqualTo(1)` / `>= 3`; no exact set checked; not evidence-derived. Rewrite. |
| existing: SingleGenome_ReturnsEmpty | 🔁 Duplicate | Same as S1; will be the canonical S1 in new file. Remove from old file. |
| M1 GoldenVector | ❌ Missing | |
| M2 ContiguousInOneGenomeOnly | ❌ Missing | |
| M3 RepeatedLabels | ❌ Missing | |
| M4 MinClusterSize | ❌ Missing | |
| M5 AllMembersContiguous | ❌ Missing | |
| S1 SingleGenome | ❌ Missing (in canonical file) | |
| S2 IdenticalOrderThreeGenomes | ❌ Missing | |
| S3 ThreeGenomes_OnlyAllConserved | ❌ Missing | |
| C1 Determinism | ❌ Missing | |
| C2 NoConservedCluster | ❌ Missing | |
| C3 NullArguments | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindConservedClusters_Tests.cs` — all COMPGEN-CLUSTER-001 cases.
- **Remove:** the `#region FindConservedClusters Tests` block (both tests) from `ComparativeGenomicsTests.cs` — weak/duplicate, superseded by the canonical file.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ComparativeGenomics_FindConservedClusters_Tests.cs` | Canonical COMPGEN-CLUSTER-001 | 12 |
| `ComparativeGenomicsTests.cs` | Other ComparativeGenomics methods (FindConservedClusters region removed) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented golden-vector exact-set test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented one-genome-only negative test ({2,3}) | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented repeated-label window test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented minClusterSize=4 filter test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented no-foreign-group window test | ✅ Done |
| 6 | S1 | ❌ Missing | Implemented single-genome empty test | ✅ Done |
| 7 | S2 | ❌ Missing | Implemented identical-order three-genome test | ✅ Done |
| 8 | S3 | ❌ Missing | Implemented split-in-third-genome test | ✅ Done |
| 9 | C1 | ❌ Missing | Implemented determinism test | ✅ Done |
| 10 | C2 | ❌ Missing | Implemented no-cluster empty test | ✅ Done |
| 11 | C3 | ❌ Missing | Implemented null-arg throw tests | ✅ Done |
| 12 | old weak/dup | ⚠/🔁 | Removed FindConservedClusters region from ComparativeGenomicsTests.cs | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact set match against Source 1 Example 1 |
| M2 | ✅ | {2,3} excluded |
| M3 | ✅ | Repeated label window reported |
| M4 | ✅ | Only size ≥ 4 sets returned |
| M5 | ✅ | Foreign group breaks the window |
| S1 | ✅ | Single genome empty |
| S2 | ✅ | Identical order all windows conserved |
| S3 | ✅ | Split-in-third excluded |
| C1 | ✅ | Deterministic |
| C2 | ✅ | No cluster empty |
| C3 | ✅ | Null throws |

All in-scope cases ✅ (11 cases + cleanup). Count of ✅ = total in-scope cases.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | `maxGap` parameter retained for API/MCP compatibility; validated behaviour is the strict gap-free common-interval model (gene-teams gapped extension not implemented — source not retrievable). API-shape only. | Method signature |

---

## 7. Open Questions / Decisions

1. **Decision:** Implement the strict common-interval model (Uno & Yagiura 2000; Heber & Stoye 2001) as the canonical, validated behaviour. `maxGap` is kept in the signature for backward/MCP compatibility but is not used to relax the strict contiguity contract (Assumption 1). The gene-teams δ-gap generalisation (Bergeron, Corteel & Raffinot 2002) is out of scope because its source could not be retrieved in this session.
2. **Decision:** Cluster identity is the set of ortholog-group labels (deduplicated). Output order of clusters is sorted canonically for determinism (INV-05).
