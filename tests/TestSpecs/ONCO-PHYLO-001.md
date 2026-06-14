# Test Specification: ONCO-PHYLO-001

**Test Unit ID:** ONCO-PHYLO-001
**Area:** Oncology
**Algorithm:** Tumor Phylogeny Reconstruction — clonal tree from CCF clusters (sum rule + lineage precedence)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Popic et al. (2015), LICHeE, *Genome Biology* 16:91 | 1 | https://doi.org/10.1186/s13059-015-0647-8 (PDF https://arxiv.org/pdf/1412.8574) | 2026-06-15 |
| 2 | Zheng et al. (2022), PICTograph, *Bioinformatics* 38(15):3677–3683 | 1 | https://doi.org/10.1093/bioinformatics/btac367 (https://pmc.ncbi.nlm.nih.gov/articles/PMC9344857/) | 2026-06-15 |

### 1.2 Key Evidence Points

1. Edge `(u→v)` valid iff for every sample i: `u.CCF[i] ≥ v.CCF[i] − ϵ` and `u.CCF[i]=0 ⇒ v.CCF[i]=0` (ancestor ≥ descendant; presence pattern) — Popic 2015 Eq. 2; Zheng 2022 lineage precedence.
2. Sum rule: for every node u and every sample i, `Σ_{children v} v.CCF[i] ≤ u.CCF[i] + ϵ` — Popic 2015 Eq. 5; Zheng 2022 sum condition.
3. Constraint (1): a cluster present in more samples cannot descend from one present in fewer — Popic 2015.
4. Trunk = clusters on the path from the root that are present (CCF>0) across all samples; branches = the rest — Popic 2015 (common predecessor present in all samples).

### 1.3 Documented Corner Cases

- Sibling groups whose CCFs sum above the parent CCF cannot all be its children (Popic 2015): forces re-placement/chain.
- Private (single-sample) clusters are under-constrained → deterministic tie-break required (Popic 2015).
- Inequality (not equality) tolerates unobserved branches / noise (Popic 2015).

### 1.4 Known Failure Modes / Pitfalls

1. Treating a child whose CCF exceeds the parent as valid — violates Eq. 2 (Popic 2015).
2. Allowing children CCF sum to exceed the parent — violates Eq. 5 (Popic 2015; Zheng 2022).
3. Placing a multi-sample cluster under a single-sample cluster — violates constraint (1) (Popic 2015).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ReconstructPhylogeny(IReadOnlyList<CcfCluster>, tolerance)` | OncologyAnalyzer | Canonical | Builds the clonal tree from CCF clusters (Eq. 2 + Eq. 5). |
| `IdentifyTrunkMutations(ClonalPhylogeny)` | OncologyAnalyzer | Canonical | Clusters on the trunk (root→first branch point) present in all samples. |
| `IdentifyBranchMutations(ClonalPhylogeny)` | OncologyAnalyzer | Canonical | Non-trunk (subclonal) clusters. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | On every edge (u→v): `u.CCF[i] ≥ v.CCF[i] − ϵ` for all samples i. | Yes | Popic 2015 Eq. 2 |
| INV-2 | For every node u, sample i: `Σ_children v.CCF[i] ≤ u.CCF[i] + ϵ`. | Yes | Popic 2015 Eq. 5; Zheng 2022 |
| INV-3 | Result is a single rooted tree: every cluster has exactly one parent; no cycles; the (synthetic) root has CCF=1 in all samples. | Yes | Popic 2015 (spanning tree) |
| INV-4 | Trunk ∩ Branch = ∅ and Trunk ∪ Branch = all input clusters. | Yes | Popic 2015 (partition) |
| INV-5 | Deterministic: same input ⇒ same tree (deepest-valid-ancestor, id tie-break). | Yes | ASSUMPTION (tie-break); constraints leave a valid set, choice documented |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Linear chain | Single sample, CCFs A=1.0,B=0.6,C=0.3 | Edges Normal→A, A→B, B→C | Popic 2015 Eq.2 |
| M2 | Branching | 2 samples, A=[1,1] trunk, B=[0.6,0], C=[0,0.7] | Edges Normal→A, A→B, A→C (B,C siblings) | Popic 2015 (1)(2)(3) |
| M3 | Sum-rule forces chain | Single sample, A=1.0,B=0.6,C=0.6 | C cannot be A's 2nd child (0.6+0.6>1.0) → Normal→A→B→C | Popic 2015 Eq.5 |
| M4 | Trunk identification | M2 tree | Trunk = {A} | Popic 2015 |
| M5 | Branch identification | M2 tree | Branches = {B, C} | Popic 2015 |
| M6 | INV-1 holds | M1 & M2 trees | every edge ancestor CCF ≥ descendant CCF per sample | Popic 2015 Eq.2 |
| M7 | INV-2 holds (property) | random valid CCF inputs (fixed seed) | per-node children CCF sum ≤ parent CCF | Popic 2015 Eq.5 |
| M8 | Single cluster | one cluster A=1.0 | Normal→A; trunk={A}; branches={} | Popic 2015 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty input | no clusters | tree = root only; trunk={}, branches={} | boundary |
| S2 | Tolerance admits near-violation | A=1.0, B=1.05 with ε=0.1 | B is child of A (1.0 ≥ 1.05−0.1) | Popic 2015 ϵ |
| S3 | Tolerance ε=0 rejects | A=1.0, B=1.05 | B cannot be A's child via Eq.2; attaches to root | strict |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism | M2 input run twice | identical edge set | INV-5 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No `ReconstructPhylogeny` / `TumorEvolutionAnalyzer` / phylogeny code or tests exist. New unit. Searched `src/` and `tests/` for `ReconstructPhylogeny`, `Phylo`, `ClonalTree` — none found. By-area definition (checklist §ONCO-PHYLO-001) references `TumorEvolutionAnalyzer`; per task instruction the methods are added to the existing `OncologyAnalyzer` class (sibling units ONCO-CLONAL-001 etc. live there).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M8 | ❌ Missing | new unit |
| S1–S3 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ReconstructPhylogeny_Tests.cs` — all cases.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_ReconstructPhylogeny_Tests.cs | canonical | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented (property, seed 42) | ✅ Done |
| 8 | M8 | ❌ Missing | implemented | ✅ Done |
| 9 | S1 | ❌ Missing | implemented | ✅ Done |
| 10 | S2 | ❌ Missing | implemented | ✅ Done |
| 11 | S3 | ❌ Missing | implemented | ✅ Done |
| 12 | C1 | ❌ Missing | implemented | ✅ Done |
| 13 | Validation (null/NaN/range/ragged) | ❌ Missing | implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | linear chain edges asserted |
| M2 | ✅ | branching edges asserted |
| M3 | ✅ | sum-rule chain asserted |
| M4 | ✅ | trunk={A} |
| M5 | ✅ | branches={B,C} |
| M6 | ✅ | INV-1 per edge |
| M7 | ✅ | INV-2 property, seed 42 |
| M8 | ✅ | single cluster |
| S1 | ✅ | empty input |
| S2 | ✅ | tolerance admits |
| S3 | ✅ | strict rejects |
| C1 | ✅ | determinism |
| Validation | ✅ | null/NaN/range/ragged exceptions |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Deepest-valid-ancestor + id tie-break for under-constrained placement | INV-5, M1–M3, determinism |
| 2 | Default noise margin ε = 0 (configurable) | S2, S3 |

---

## 7. Open Questions / Decisions

1. By-area definition names `TumorEvolutionAnalyzer`; methods placed in `OncologyAnalyzer` per the session instruction and to match sibling Oncology units. Recorded in §5.1.
