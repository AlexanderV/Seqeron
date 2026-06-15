# Test Specification: PHYLO-BOOT-001

**Test Unit ID:** PHYLO-BOOT-001
**Area:** Phylogenetic
**Algorithm:** Phylogenetic Bootstrap Analysis (Felsenstein's Bootstrap Proportions)
**Status:** ☑ Validated (PHYLO-BOOT-001 session, 2026-06-15)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Felsenstein J. (1985). Confidence Limits on Phylogenies. Evolution 39(4):783–791 | 1 | https://doi.org/10.1111/j.1558-5646.1985.tb00420.x (text via https://www.osti.gov/biblio/6044842) | 2026-06-13 |
| 2 | Lemoine et al. (2018). Renewing Felsenstein's Phylogenetic Bootstrap. Nature 556:452–456 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC6030568/ | 2026-06-13 |
| 3 | Biopython `Bio.Phylo.Consensus` (`bootstrap`, `get_support`) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/Phylo/Consensus.py | 2026-06-13 |

### 1.2 Key Evidence Points

1. Each bootstrap replicate resamples the alignment **columns (sites) with replacement** while keeping all taxa, producing a pseudo-alignment of the **same length** as the original. — Felsenstein (1985) abstract; Lemoine (2018); Biopython `bootstrap_trees`.
2. Support for a clade = **proportion of bootstrap trees that contain a clade with the identical terminal (leaf) set**. — Lemoine (2018); Biopython `get_support` (`(t+1)*100/size`, here as a fraction).
3. The entities scored are the **non-trivial clades of the reference tree** built from the original data; clades appearing only in replicates are not reported. — Lemoine (2018); Biopython `find_clades(terminal=False)`.
4. A group appearing in 100% of replicates has support 1.0; thresholds like 95%/70% are interpretive, not part of the computation. — Felsenstein (1985).
5. Clades are compared by **leaf-name set only** (branch lengths and internal labels irrelevant). — Biopython `_clade_to_bitstr`.

### 1.3 Documented Corner Cases

- Reference tree fixes the clade set; replicate-only clades are ignored (Lemoine).
- Binary per-replicate scoring: exact match counted, otherwise not (Lemoine).
- Pseudo-alignment length must equal original length (Felsenstein; Biopython).
- All-identical / zero-distance alignments reproduce one topology every replicate → support 1.0 for every reported clade.

### 1.4 Known Failure Modes / Pitfalls

1. Non-deterministic results without a fixed RNG seed (resampling is randomized). — Felsenstein (1985).
2. Comparing clades by anything other than the leaf-name set produces wrong support. — Biopython `get_support`.
3. Fewer than 2 sequences cannot form a tree; `BuildTree` requires ≥2. — repository contract / source-aligned.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `Bootstrap(sequences, replicates, distanceMethod, treeMethod, seed)` | PhylogeneticAnalyzer | Canonical | Deep evidence-based testing; randomized → fixed documented seed |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned support value is in [0,1] | Yes | Support = count/replicates, count∈[0,replicates] — Biopython `get_support`; Lemoine "proportion" |
| INV-2 | Each support value equals k/replicates for some integer k (0 ≤ k ≤ replicates) | Yes | Count is an integer number of matching replicate trees — Biopython |
| INV-3 | Result keys equal exactly the non-trivial clades of the reference (original-data) tree | Yes | Biopython scores `find_clades(terminal=False)`; Lemoine "branches of the reference tree" |
| INV-4 | Deterministic: identical (sequences, replicates, methods, seed) → identical result | Yes | Randomized algorithm requires fixed seed for reproducibility — Felsenstein (1985) |
| INV-5 | A clade present in every replicate has support exactly 1.0 | Yes | Felsenstein (1985): 100% occurrence → P=1 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Two-group support = 1.0 | A=B=`AAAAAAAAAA`, C=D=`GGGGGGGGGG`, UPGMA, JC, 100 reps, seed 42. Distances invariant under column resampling → same topology every replicate | support({A,B})=1.0 and support({C,D})=1.0 (Within 1e-10) | Felsenstein (1985); derivation in Evidence Dataset 1 |
| M2 | All values in [0,1] | Mixed alignment, fixed seed | every value ≥ 0.0 and ≤ 1.0 | INV-1; Biopython `get_support` |
| M3 | Quantized to k/replicates | replicates=20, fixed seed | every value × 20 is (within 1e-9) a non-negative integer ≤ 20 | INV-2; Biopython count/size |
| M4 | Keys = reference clades | Build reference tree from same data; collect non-trivial clades | result keys set == reference non-trivial clade set | INV-3; Biopython `find_clades(terminal=False)` |
| M5 | Determinism (same seed) | Run twice with seed 42, identical inputs | both dictionaries equal key-for-key and value-for-value | INV-4; Felsenstein randomized method |
| M6 | All-identical sequences | A=B=C=`ACGTACGT`, fixed seed | every reported clade has support 1.0 | Evidence Dataset 2; INV-5 |
| M7 | NeighborJoining branch | Two-group dataset, `treeMethod=NeighborJoining`, JC, 50 reps, seed 42 | support({A,B})=1.0 and support({C,D})=1.0 | Felsenstein (1985): support procedure is tree-method-agnostic; distances invariant under column resampling (exercises the NJ branch of `treeMethod`) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null sequences | `Bootstrap(null, ...)` | throws `ArgumentException` (or `ArgumentNullException`) | BuildTree requires non-null ≥2 |
| S2 | Fewer than 2 sequences | single-sequence dictionary | throws `ArgumentException` | Cannot build a tree from 1 taxon |
| S3 | Replicates < 1 | replicates = 0 | throws `ArgumentException` | Denominator must be ≥ 1 |
| S3b | Negative replicates | replicates = −5 | throws `ArgumentException` | Boundary below the ≥1 contract |
| S4 | Different seeds may differ but stay valid | seed 1 vs seed 7 on a partly-informative alignment | both results satisfy INV-1/INV-2 | Determinism is per-seed |
| S5 | Unequal-length sequences | A=`ACGT`, B=`ACG` | throws `ArgumentException` | Bootstrap resamples alignment columns; unequal lengths are not an alignment (Felsenstein 1985; surfaces from `BuildTree`) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Replicate alignment length preserved | Verify each replicate resamples exactly `alignmentLength` columns (indirect via support quantization at known length) | covered by M3 | Same-length resampling (Felsenstein) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_TreeComparison_Tests.cs` — `#region Bootstrap Tests (Different Scope - Kept for Coverage)` (2 tests: `Bootstrap_ReturnsSupportsInValidRange`, `Bootstrap_DistinctGroups_HighSupport`). Both are weak (permissive `InRange`, `Any(v => v >= 0.5)`, no assertion messages, no exact evidence values, non-deterministic w.r.t. seed exposure).
- No canonical `PhylogeneticAnalyzer_Bootstrap_Tests.cs` existed before this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | No exact-value support test exists |
| M2 | ⚠ Weak | `Bootstrap_ReturnsSupportsInValidRange` uses `InRange` w/o messages; rewrite from scratch in canonical file |
| M3 | ❌ Missing | No quantization check |
| M4 | ❌ Missing | No reference-clade-key check |
| M5 | ❌ Missing | No determinism check; seed was not exposed |
| M6 | ❌ Missing | No all-identical case |
| S1 | ❌ Missing | No null validation test (method also lacked validation) |
| S2 | ❌ Missing | No <2 sequence test |
| S3 | ❌ Missing | No replicates<1 test |
| S4 | ❌ Missing | No multi-seed test (seed not exposed) |
| C1 | ❌ Missing | Covered indirectly by M3 |
| `Bootstrap_DistinctGroups_HighSupport` | ⚠ Weak | Permissive `Any(v >= 0.5)`; superseded by exact M1 — remove |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_Bootstrap_Tests.cs` — all M/S/C cases for `Bootstrap`.
- **Remove:** the entire `#region Bootstrap Tests` from `PhylogeneticAnalyzer_TreeComparison_Tests.cs` (2 weak tests) — bootstrap is PHYLO-BOOT-001's scope, not tree comparison's.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PhylogeneticAnalyzer_Bootstrap_Tests.cs` | Canonical PHYLO-BOOT-001 fixture | 13 |
| `PhylogeneticAnalyzer_TreeComparison_Tests.cs` | PHYLO-COMPARE — bootstrap region removed | (unchanged otherwise) |

> **Validation update (2026-06-15, PHYLO-BOOT-001 session):** added M7 (NeighborJoining
> branch — `treeMethod` had no coverage), S3b (negative replicates), and S5 (unequal-length
> sequences — a documented Stage-A edge case that was previously untested). Fixture 11→13 tests.
> Note: the §5.4 "11" above was the count claimed at authoring; the file actually held 10 before
> this session (M1–M6, M4/M5, S1–S4 → 10), now 13.

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact 1.0 support test (two groups, seed 42) | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote from scratch in canonical file with messages | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented quantization check (k/replicates) | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented reference-clade-key equality test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented determinism (same seed) test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented all-identical → 1.0 test | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented null validation test | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented <2 sequence test | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented replicates<1 test | ✅ Done |
| 10 | S4 | ❌ Missing | Implemented multi-seed validity test | ✅ Done |
| 11 | C1 | ❌ Missing | Folded into M3 (same-length resampling implies integer k) | ✅ Done |
| 12 | `Bootstrap_DistinctGroups_HighSupport` | ⚠ Weak | Removed (superseded by exact M1) | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact-value support test |
| M2 | ✅ Covered | Rewritten range test with messages |
| M3 | ✅ Covered | Quantization test |
| M4 | ✅ Covered | Reference-clade-key test |
| M5 | ✅ Covered | Determinism test |
| M6 | ✅ Covered | All-identical test |
| S1 | ✅ Covered | Null validation |
| S2 | ✅ Covered | <2 sequences |
| S3 | ✅ Covered | replicates<1 |
| S4 | ✅ Covered | Multi-seed validity |
| C1 | ✅ Covered | Via M3 |
| Old weak region | ✅ Covered | Removed from TreeComparison file |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Rooted-clade scoring (subtree leaf-set) rather than unrooted bipartitions; matches Biopython `get_support` | INV-3, M4 |
| 2 | Support returned as proportion in [0,1] rather than percentage (×100 gives the published value) | INV-1, M2 |

Both are units/representation choices, not correctness-affecting on which clades are reported or their ranking; documented in Evidence Assumptions and algorithm doc §5.4.

---

## 7. Open Questions / Decisions

1. **Decision:** Added an explicit `seed` parameter to `Bootstrap` so tests are deterministic with a documented seed (default 42 preserves prior behavior for existing callers). The original signature in the checklist (`Bootstrap(sequences, nReplicates, treeMethod)`) omits `distanceMethod`; the implemented signature keeps `distanceMethod` (present in code and required to compute distances) and appends `seed`. Conflict noted; checklist entry is workflow-only.
2. **Decision:** No suffix-tree reuse — this is a resampling + distance-tree + clade-matching unit, not a substring-search unit. N/A.
