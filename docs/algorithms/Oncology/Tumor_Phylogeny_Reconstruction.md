# Tumor Phylogeny Reconstruction

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Clonal evolution |
| Test Unit ID | ONCO-PHYLO-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Complete (LICHeE port) |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Reconstructs a rooted clonal (tumor) phylogeny — the ancestor/descendant ordering of subclones — from cancer cell fraction (CCF) clusters. Each cluster is a candidate clone with a CCF in every sequenced sample; the clusters themselves are produced upstream (ONCO-CCF-001). The method builds the tree by applying two deterministic lineage constraints from the multi-sample perfect-phylogeny model — lineage precedence (an ancestor's CCF is at least its descendant's) and the sum rule (children CCFs may not exceed the parent's) [1][2]. The implementation is a port of LICHeE's tree search (Popic et al. 2015; reference code github.com/viq854/lichee, cell-prevalence mode `-cp` with pre-computed clusters): build the constraint network, enumerate every spanning tree that satisfies the sum rule (Gabow–Myers), rank them by the LICHeE error score and return the top tree. It is not a probabilistic joint clustering method (PyClone/PhyloWGS/CITUP).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A tumor evolves clonally: subclones arise by acquiring additional somatic mutations on top of an ancestral clone. Bulk/multi-region sequencing yields, for each mutation cluster, a cancer cell fraction (the fraction of tumor cells carrying it) in each sample. Under the infinite-sites / perfect-phylogeny assumption (a mutation arises once and is never lost), these CCFs constrain which clones can be ancestors of which others [1].

### 2.2 Core Model

Let cluster `v` have CCF `v.CCF[i]` in sample `i`, and let `ε` be a noise margin. The clonal tree `T` (rooted at a synthetic normal node with CCF = 1 in all samples) must satisfy:

- **Lineage precedence (ancestor ≥ descendant), Eq. 2 [1]:** an edge `u→v` is admissible only if for every sample `i`, `u.CCF[i] ≥ v.CCF[i] − ε`, and (presence pattern) `u.CCF[i] = 0 ⇒ v.CCF[i] = 0`. Equivalently, "the CCF of any mutation cannot exceed the CCF of its ancestor" [2].
- **Sum rule, Eq. 5 [1]:** for every node `u` and every sample `i`, `Σ_{v : (u→v)∈T} v.CCF[i] ≤ u.CCF[i] + ε`. Equivalently, "the CCF of an ancestral clone must be greater than or equal to the sum of CCFs of its descendants" [2]. This is the cell-fraction analogue of the pigeonhole principle: disjoint sibling subclones cannot together occupy more cells than their common parent [1].

These constraints define a *set* of valid spanning trees; they do not uniquely determine one (private/single-sample clusters are under-constrained) [1]. LICHeE enumerates the set over a constraint network and ranks it by the error score `√(Σ_u Σ_i max(0, Σ_children v.CCF[i] − u.CCF[i])²)` (`PHYTree.computeErrorScore`) [3]; the set may be empty, in which case no tree is reported.

- **Trunk:** "alterations that are in the trunk of the tree must be present in all cells of the tumour" [4] — a cluster is truncal iff it lies on the root path with CCF = 1 (≥ 1 − ε) in every sample. A single-child descendant with CCF < 1 is subclonal (its parent keeps a residual population without it).

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Infinite sites / perfect phylogeny: each mutation arises once and is never lost [1]. | Convergent or back-mutation breaks the ancestor ≥ descendant ordering; reconstructed edges may be wrong. |
| ASM-02 | CCFs are correctly estimated and clustered upstream (ONCO-CCF-001). | Mis-clustered CCFs propagate to wrong edges; this unit does not re-estimate them. |
| ASM-03 | Among equal-score trees the first one enumerated by LICHeE's Gabow–Myers search is returned; presence profiles are visited in first-appearance order (LICHeE: Java `HashMap` order). | A different (still valid, equal-score) tree could be returned for another profile order; topology among ambiguous private clusters is not unique. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every edge `u→v` satisfies `u.CCF[i] ≥ v.CCF[i] − ε` for all samples `i`. | Network edges only exist when LICHeE `checkAndAddEdge` admits them (Eq. 2). |
| INV-02 | For every node (root included), per sample, children CCF sum ≤ parent CCF + ε. | Every tree edge passes `PHYTree.checkConstraint`; no tree ⇒ `InvalidOperationException` / `TryReconstructPhylogeny` false. |
| INV-03 | The result is a single rooted tree: every cluster has exactly one parent, no cycles. | Gabow–Myers grows spanning trees rooted at the normal node [3]. |
| INV-04 | Trunk and Branch partition the clusters (disjoint, union = all). | `IdentifyBranchMutations` = clusters ∉ trunk set. |
| INV-05 | Deterministic: identical input ⇒ identical edge set. | Fixed node order, LIFO enumeration and stable ranking. |

### 2.5 Comparison with Related Methods

| Aspect | This method | PyClone / PhyloWGS / CITUP |
|--------|-------------|----------------------------|
| Inference | LICHeE exhaustive search over the constraint network, ranked by error score | Probabilistic (MCMC / ILP) joint clustering + tree |
| Clustering | Consumes clusters (ONCO-CCF-001) | Performs clustering jointly |
| Output | Top-ranked tree + number of valid trees + error score | Posterior over trees / optimal ILP tree |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `clusters` | `IReadOnlyList<CcfCluster>` | required | CCF clusters to place | each `CcfPerSample` same non-zero length; values in [0,1]; unique ids |
| `tolerance` | `double` | `0.0` | Noise margin ε for Eq. 2 and Eq. 5 | ≥ 0, not NaN |
| `clusters` (summary overloads, F42) | `IReadOnlyList<CcfClusterSummary>` | required | Centroid CCF + per-sample SD (population, divisor n) + member count n, e.g. `CcfClusterSummary.FromMembers(id, memberCcfs)` (LICHeE `Cluster.recomputeCentroidAndStdDev`) | centroids as above; SD finite ≥ 0, same length; n ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `RootId` | `int` | Synthetic normal root id (distinct from all cluster ids) |
| `Clusters` | `IReadOnlyList<CcfCluster>` | Input clusters in input order |
| `Edges` | `IReadOnlyList<ClonalEdge>` | parent→child edges, one per non-root cluster |
| `SampleCount` | `int` | Samples per cluster |
| `Tolerance` | `double` | ε used (also the trunk threshold 1 − ε) |
| `ErrorScore` | `double` | LICHeE error score of the returned tree |
| `ValidTreeCount` | `int` | Valid trees enumerated (capped at 100 000) |
| `UsedCompleteNetwork` | `bool` | The complete-network fallback was needed |

`IdentifyTrunkMutations` returns trunk cluster ids (root downward); `IdentifyBranchMutations` returns the remaining (subclonal) ids in input order.

### 3.3 Preconditions and Validation

No valid tree → `InvalidOperationException` (`TryReconstructPhylogeny` returns false). Null `clusters` or null cluster CCF list → `ArgumentNullException`. Empty/ragged CCF lists, NaN/out-of-[0,1] CCF, or duplicate ids → `ArgumentException`. Negative/NaN `tolerance` → `ArgumentOutOfRangeException`. Empty cluster list → a root-only tree (no edges, no trunk/branch). CCF indexing is 0-based per sample; ids are caller-assigned integers.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; pick a synthetic root id below the minimum cluster id; the root has CCF = 1 in every sample (LICHeE `-cp`: `VAF_MAX = 1`).
2. Constraint network (`PHYNetwork` constructor): level = number of samples with CCF > 0; root level = k+1. `checkAndAddEdge(u,v)` over samples: stop at the first sample with `u = 0 ≠ v`; count samples with `u ≥ v − ε`; if both directions pass, keep u→v iff `Σ max(0, v−u) < Σ max(0, u−v)` (otherwise v→u). Edges: all pairs within a presence profile; each level to the next non-empty lower level; then every node without a parent is tried against levels `level+2 … k+1`, else linked to the root.
3. Enumerate spanning trees from the root (Gabow–Myers `grow`), rejecting a partial tree as soon as a parent violates Eq. 5 with margin ε; stop at 100 000 trees (`MAX_NUM_TREES`) or 10⁸ grow calls.
4. Rank by the error score (stable); return the first minimum.
5. If no tree was found, rebuild with the complete network (`ALL_EDGES`: every higher level to every lower level ≥ 1) and search again; if still none, report failure.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Lineage precedence `u.CCF[i] ≥ v.CCF[i] − ε` and presence `u.CCF[i]=0 ⇒ v.CCF[i]=0` [1] Eq. 2.
- Sum rule `Σ_children v.CCF[i] > u.CCF[i] + ε ⇒ reject` [1] Eq. 5 (`PHYTree.checkConstraint`).
- Per-cluster edge margin (summary overloads, F42; LICHeE `PHYNetwork.getAAFErrorMargin`): the Eq. 2 test of u→v in sample i uses `max(ε, se_u,i + se_v,i)`, `se = 1.96·sd_i/√n` for a cluster (sd = 0 where the centroid is 0, `PHYNode.getStdDev`) and `se = ε` for the root. It applies to every `checkAndAddEdge` call (within-profile, inter-level, orphan attachment, complete network); the sum rule (Eq. 5, `PHYTree.checkConstraint`) keeps the static ε. With all SD = 0 this is the static network.
- Default `ε = 0` (strict); source defaults are ϵ (LICHeE) / ε₁=0.1, ε₂=0.2 (PICTograph) [1][2], exposed via `tolerance`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ReconstructPhylogeny` | O(n² · k) network + O(T · n² · k) search | O(n² + n · k) | T = number of spanning trees visited (exponential in the worst case, capped at 10⁵ trees / 10⁸ grow calls as in LICHeE). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.PhylogenyHeterogeneity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.PhylogenyHeterogeneity.cs) (private `LicheeNetwork`).

- `OncologyAnalyzer.ReconstructPhylogeny(IReadOnlyList<CcfCluster>, double)`: top-ranked LICHeE tree; throws `InvalidOperationException` when none is valid.
- `OncologyAnalyzer.TryReconstructPhylogeny(IReadOnlyList<CcfCluster>, out ClonalPhylogeny, double)`: non-throwing variant.
- `OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(IReadOnlyList<CcfClusterSummary>, double)` / `TryReconstructPhylogenyFromClusterSummaries(...)` (F42): LICHeE with per-cluster `1.96·sd/√n` edge margins; `CcfClusterSummary.FromMembers` builds centroid/SD/n from member CCFs.
- `OncologyAnalyzer.IdentifyTrunkMutations(ClonalPhylogeny)`: truncal clusters (root path, CCF ≥ 1 − ε in every sample).
- `OncologyAnalyzer.IdentifyBranchMutations(ClonalPhylogeny)`: subclonal/branch clusters (the rest).

### 5.2 Current Behavior

Cross-checked against the original `lichee.jar` (PHYNetwork + getLineageTrees + evaluateLineageTrees driven through a harness, clusters marked robust, profiles in first-appearance order): 4 503 / 4 503 random inputs (1–4 samples, 1–11 clusters, ε ∈ {0, 0.02, 0.05, 0.1, 0.2}) give identical feasibility, network mode, tree count, top tree and bit-identical error score (incl. 83 capped at 10⁵ trees, 81 complete-network fallbacks, 1 623 infeasible). The former greedy "deepest valid ancestor" build returned a sum-rule-violating tree (root fallback / FP budget debit) on 939 / 2 280 random inputs. **Per-cluster margins (F42):** the same harness, with SNVGroups built from member CCF rows and `Cluster.recomputeCentroidAndStdDev`, on 9 000 random member-level inputs (1–3 samples, 2–6 clusters of 2–5 members, ε ∈ {0, 0.02, 0.05, 0.1}): 9 000 / 9 000 identical (feasibility, network mode, tree count, top tree, error score) and 57 594 / 57 594 per-sample centroids and SDs bit-identical; the margins changed the result vs the static ε on 16 of them (locked: p1, t02150, t01061, t00298). Because a margin edge u→v with `v − u > ε` in some sample can never satisfy the static-ε sum rule, the margins act through edge orientation (`checkAndAddEdge` keeps the smaller one-sided excess once both directions pass) and orphan attachment, often ending in the complete network or no tree. **Search reuse:** the repository suffix tree was evaluated and is **not** applicable — this unit performs no substring/pattern search; it is a numeric constraint-satisfaction tree build over CCF vectors.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Lineage precedence edge rule `u.CCF[i] ≥ v.CCF[i] − ε` and presence pattern `u=0 ⇒ v=0` [1] Eq. 2.
- Sum rule `Σ_children v.CCF[i] ≤ u.CCF[i] + ε` per node, per sample [1] Eq. 5; [2] sum condition.
- Rooted spanning tree over all clusters [1]; LICHeE constraint network, Gabow–Myers enumeration, error-score ranking, complete-network fallback [3].

**Deviations from LICHeE:**

- `ReconstructPhylogeny(CcfCluster…)` carries point CCFs only, so it uses the static margin ε (= LICHeE with zero-variance clusters); the per-cluster `1.96·sd/√n` margins are available through `ReconstructPhylogenyFromClusterSummaries` (F42).
- LICHeE's `fixNetwork` step (dropping non-robust clusters when no tree exists) is not applied — every cluster is kept (clusters treated as robust); the caller receives a failure instead of a silently reduced tree.
- Hitting the 10⁸ grow-call cap stops the search (LICHeE continues from a partially unwound state).
- Default `ε = 0` (LICHeE default 0.1); pass `tolerance` to reproduce it.

**Not implemented:**

- Probabilistic clone clustering / posterior over trees (PyClone, PhyloWGS, CITUP); **users should rely on:** ONCO-CCF-001 for CCF clustering and dedicated external tools for posterior tree inference.
- CNA-aware multiplicity corrections to CCF inside this unit (it consumes already-corrected CCF clusters); **users should rely on:** `DeriveMultiplicity` (ONCO-ASCAT-001) + `EstimateCcf` with the local copy number (ONCO-CCF-001) upstream.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Profile order = first appearance (LICHeE: HashMap order) | Assumption | Tie-break among equal-score trees | accepted | ASM-03 |
| 2 | Default ε = 0 | Assumption | Stricter admissibility than LICHeE's 0.1 | accepted | configurable via `tolerance` |
| 3 | No `fixNetwork` cluster removal | Deviation | Infeasible inputs fail instead of dropping clusters | accepted | needs SNV-level robustness data |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty cluster list | root-only tree; trunk={}, branches={} | nothing to place |
| Single cluster | root→cluster; trunk={cluster} iff CCF = 1 in every sample | trunk definition [4] |
| Two equal-CCF clusters under one parent | only one nests as child; the other chains below it (LICHeE orients later→earlier) | sum rule (Eq. 5) [1] |
| No tree satisfies Eq. 5 (e.g. ε=0, A=[0.5,0.5], B=[0.55,0]) | `InvalidOperationException`; `TryReconstructPhylogeny` false | LICHeE reports no valid tree |
| Cluster present in more samples than candidate parent | not a descendant of that parent (presence pattern) | constraint (1) [1] |
| Null/ragged/NaN/out-of-range CCF, duplicate id | exception | validation |

### 6.2 Limitations

Assumes correct upstream CCF clustering and the infinite-sites model; does not model copy-number-driven CCF distortions, mutation loss, or convergent evolution; for ambiguous private clusters the topology is one of several valid trees, not a unique biological truth.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var clusters = new[]
{
    new OncologyAnalyzer.CcfCluster(1, new[] { 1.0, 1.0 }), // trunk
    new OncologyAnalyzer.CcfCluster(2, new[] { 0.6, 0.0 }), // private to sample 1
    new OncologyAnalyzer.CcfCluster(3, new[] { 0.0, 0.7 }), // private to sample 2
};
OncologyAnalyzer.ClonalPhylogeny tree = OncologyAnalyzer.ReconstructPhylogeny(clusters);
// Edges: root→1, 1→2, 1→3. Trunk = {1}; Branches = {2, 3}.
```

**Numerical walk-through:** levels: 1 (present in 2 samples) = level 2; 2 and 3 = level 1 (different profiles). Network: root→1; 1→2, 1→3 (1 ≥ child, presence OK); 2↔3 inadmissible (presence). Only spanning tree: root→1→{2,3}; sum rule under 1: s1 0.6 ≤ 1, s2 0.7 ≤ 1. LICHeE: 1 valid tree, error 0. Trunk {1} (CCF 1 in both samples).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_ReconstructPhylogeny_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_ReconstructPhylogeny_Tests.cs) — covers `INV-01`, `INV-02`, `INV-05`
- Evidence: [ONCO-PHYLO-001-Evidence.md](../../../docs/Evidence/ONCO-PHYLO-001-Evidence.md)
- Related algorithms: [Clonal_Subclonal_Classification](./Clonal_Subclonal_Classification.md)

## 8. References

1. Popic V, Salari R, Hajirasouliha I, Kashef-Haghighi D, West RB, Batzoglou S. 2015. Fast and scalable inference of multi-sample cancer lineages. *Genome Biology* 16:91. https://doi.org/10.1186/s13059-015-0647-8
2. Zheng L, Dang H, Niknafs N, et al. 2022. Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors (PICTograph). *Bioinformatics* 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac367
3. LICHeE reference implementation (MIT), github.com/viq854/lichee: `LICHeE/src/lineage/PHYNetwork.java`, `PHYTree.java`, `PHYNode.java`, `Parameters.java`, `LineageEngine.java`, `SNVDataStore.java`; binary `LICHeE/release/lichee.jar`. Gabow HN, Myers EW. 1978. Finding all spanning trees of directed and undirected graphs. *SIAM J Comput* 7(3):280–287.
4. Werner B et al. 2017. Detecting truly clonal alterations from multi-region profiling of tumours. *Sci Rep* 7:44991. https://doi.org/10.1038/srep44991
