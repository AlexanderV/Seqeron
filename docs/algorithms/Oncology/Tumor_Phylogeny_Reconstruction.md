# Tumor Phylogeny Reconstruction

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Clonal evolution |
| Test Unit ID | ONCO-PHYLO-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Complete (LICHeE port) |
| Last Reviewed | 2026-10-10 |

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
| ASM-03 | Among equal-score trees the first one enumerated by LICHeE's Gabow–Myers search is returned; presence-profile groups are visited in LICHeE's own order — the iteration order of its `HashMap<String, …>` of profile tags, emulated exactly (F45). | Equal-score trees are equally supported by the data; the choice is LICHeE's, not a biological preference — topology among ambiguous private clusters is not unique. |

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
| `tolerance` | `double` | `0.1` (`DefaultPhylogenyTolerance`, LICHeE `-e`; F44) | Noise margin ε for Eq. 2 and Eq. 5 | ≥ 0, not NaN |
| `removeNonRobustClusters` (summary overloads, F43) | `bool` | `true` | Run LICHeE `fixNetwork` (drop the smallest non-robust cluster while no tree exists) | robust ⇔ `RobustMemberCount ?? MemberCount` ≥ 2 |
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
| `RemovedClusterIds` | `IReadOnlyList<int>` | Clusters dropped by `fixNetwork` (removal order; absent from `Clusters`/`Edges`); empty for `ReconstructPhylogeny` (F43) |

`IdentifyTrunkMutations` returns trunk cluster ids (root downward); `IdentifyBranchMutations` returns the remaining (subclonal) ids in input order.

### 3.3 Preconditions and Validation

No valid tree → `InvalidOperationException` (`TryReconstructPhylogeny` returns false). Null `clusters` or null cluster CCF list → `ArgumentNullException`. Empty/ragged CCF lists, NaN/out-of-[0,1] CCF, or duplicate ids → `ArgumentException`. Negative/NaN `tolerance` → `ArgumentOutOfRangeException`. Empty cluster list → a root-only tree (no edges, no trunk/branch). CCF indexing is 0-based per sample; ids are caller-assigned integers.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; pick a synthetic root id below the minimum cluster id; the root has CCF = 1 in every sample (LICHeE `-cp`: `VAF_MAX = 1`).
2. Node order (F45): one group per presence profile (tag of '1'/'0' per sample); groups in the iteration order of LICHeE's `HashMap<String, …> tag2SNVs` (clusters-file loader: `containsKey` + `put` in first-appearance order, no removal; `LineageEngine` step 2 iterates `keySet()`), reproduced by `JavaStringHashMapOrder` (OpenJDK `String.hashCode`, spread `h ^ h>>>16`, bucket `(n−1) & hash`, resize at 0.75·n, bins of ≥ 9 treeified at n ≥ 64 with hash-then-`compareTo` red-black insertion, `moveRootToFront`, tree split/untreeify on resize); within a group, input order. Node ids follow this order and fix the edge-insertion order, the `checkAndAddEdge` argument order, the enumeration order, the `computeErrorScore` summation order and the `fixNetwork` scan.
   Constraint network (`PHYNetwork` constructor): level = number of samples with CCF > 0; root level = k+1. `checkAndAddEdge(u,v)` over samples: stop at the first sample with `u = 0 ≠ v`; count samples with `u ≥ v − ε`; if both directions pass, keep u→v iff `Σ max(0, v−u) < Σ max(0, u−v)` (otherwise v→u). Edges: all pairs within a presence profile; each level to the next non-empty lower level; then every node without a parent is tried against levels `level+2 … k+1`, else linked to the root.
3. Enumerate spanning trees from the root (Gabow–Myers `grow`), rejecting a partial tree as soon as a parent violates Eq. 5 with margin ε; stop at 100 000 trees (`MAX_NUM_TREES`) or 10⁸ grow calls (`MAX_NUM_GROW_CALLS`). At the cap LICHeE's `grow` returns from the current frame while its callers keep popping their stack edges, but a recursion is only entered while the counter is below the cap, so no tree is added after it: the result is the best of the trees completed within the budget, which the port returns by stopping at once (F44, jar-locked with reduced caps).
4. Rank by the error score (stable); return the first minimum.
5. If no tree was found and cluster sizes/robustness are known (summary overloads, F43): `fixNetwork` — scan nodes in id order and drop the first non-robust cluster with the strictly smallest member count (robust = ≥ `MIN_ROBUST_CLUSTER_SUPPORT` = 2 robust members), rebuild the default network, search; repeat until a tree is found or no non-robust cluster is left.
6. If still no tree, rebuild with the complete network (`ALL_EDGES`: every higher level to every lower level ≥ 1) on the remaining clusters (LICHeE rebuilds from the groups mutated by `fixNetwork`) and search again; if still none, report failure.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Lineage precedence `u.CCF[i] ≥ v.CCF[i] − ε` and presence `u.CCF[i]=0 ⇒ v.CCF[i]=0` [1] Eq. 2.
- Sum rule `Σ_children v.CCF[i] > u.CCF[i] + ε ⇒ reject` [1] Eq. 5 (`PHYTree.checkConstraint`).
- Per-cluster edge margin (summary overloads, F42; LICHeE `PHYNetwork.getAAFErrorMargin`): the Eq. 2 test of u→v in sample i uses `max(ε, se_u,i + se_v,i)`, `se = 1.96·sd_i/√n` for a cluster (sd = 0 where the centroid is 0, `PHYNode.getStdDev`) and `se = ε` for the root. It applies to every `checkAndAddEdge` call (within-profile, inter-level, orphan attachment, complete network); the sum rule (Eq. 5, `PHYTree.checkConstraint`) keeps the static ε. With all SD = 0 this is the static network.
- Default `ε = 0.1` = LICHeE `-e` default (`Parameters.VAF_ERROR_MARGIN`), unchanged by `-cp` (which only sets `VAF_MAX = MAX_ALLOWED_VAF = 1`) [3] (F44; PICTograph uses ε₁=0.1, ε₂=0.2 [2]); `tolerance: 0` gives the strict inequalities.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ReconstructPhylogeny` | O(n² · k) network + O(T · n² · k) search | O(n² + n · k) | T = number of spanning trees visited (exponential in the worst case, capped at 10⁵ trees / 10⁸ grow calls as in LICHeE). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.PhylogenyHeterogeneity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.PhylogenyHeterogeneity.cs) (private `LicheeNetwork`).

- `OncologyAnalyzer.ReconstructPhylogeny(IReadOnlyList<CcfCluster>, double)`: top-ranked LICHeE tree; throws `InvalidOperationException` when none is valid.
- `OncologyAnalyzer.TryReconstructPhylogeny(IReadOnlyList<CcfCluster>, out ClonalPhylogeny, double)`: non-throwing variant.
- `OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(IReadOnlyList<CcfClusterSummary>, double, bool removeNonRobustClusters = true)` / `TryReconstructPhylogenyFromClusterSummaries(...)` (F42/F43): LICHeE with per-cluster `1.96·sd/√n` edge margins and `fixNetwork`; `CcfClusterSummary.FromMembers` builds centroid/SD/n from member CCFs; `RobustMemberCount` (null = all members robust; 0 = LICHeE `--clustersFile` input, never robust).
- `OncologyAnalyzer.IdentifyTrunkMutations(ClonalPhylogeny)`: truncal clusters (root path, CCF ≥ 1 − ε in every sample).
- `OncologyAnalyzer.IdentifyBranchMutations(ClonalPhylogeny)`: subclonal/branch clusters (the rest).

### 5.2 Current Behavior

Cross-checked against the original `lichee.jar` (PHYNetwork + getLineageTrees + evaluateLineageTrees driven through a harness, clusters marked robust, profiles in first-appearance order): 4 503 / 4 503 random inputs (1–4 samples, 1–11 clusters, ε ∈ {0, 0.02, 0.05, 0.1, 0.2}) give identical feasibility, network mode, tree count, top tree and bit-identical error score (incl. 83 capped at 10⁵ trees, 81 complete-network fallbacks, 1 623 infeasible). The former greedy "deepest valid ancestor" build returned a sum-rule-violating tree (root fallback / FP budget debit) on 939 / 2 280 random inputs. **Per-cluster margins (F42):** the same harness, with SNVGroups built from member CCF rows and `Cluster.recomputeCentroidAndStdDev`, on 9 000 random member-level inputs (1–3 samples, 2–6 clusters of 2–5 members, ε ∈ {0, 0.02, 0.05, 0.1}): 9 000 / 9 000 identical (feasibility, network mode, tree count, top tree, error score) and 57 594 / 57 594 per-sample centroids and SDs bit-identical; the margins changed the result vs the static ε on 16 of them (locked: p1, t02150, t01061, t00298). Because a margin edge u→v with `v − u > ε` in some sample can never satisfy the static-ε sum rule, the margins act through edge orientation (`checkAndAddEdge` keeps the smaller one-sided excess once both directions pass) and orphan attachment, often ending in the complete network or no tree. **fixNetwork (F43):** 4 000 random member-level inputs with 1–3-member clusters and ~10 % non-robust members (seed 3; 1 958 with at least one cluster dropped): 4 000 / 4 000 identical (removed set, network mode, tree count, top tree, error score) against lichee.jar run with `-XX:hashCode=2` (constant identity hash ⇒ `fixNetwork`'s `HashSet<SNVGroup>` iterates in insertion order = node order); with the default JVM identity hash 3 844 / 4 000 — the 156 others differ between the two JVM settings of LICHeE itself (the rebuilt network's group order follows identity hash codes, so orphan attachment / enumeration order there is not reproducible). Locked fixtures f1–f4, f6, f7 give the same jar output under both settings. **Group order (F45):** the harness now fills a `java.util.HashMap` with `containsKey`/`put` exactly as `SNVDataStore` does (the earlier harnesses used insertion order); against it the port is identical on 21 000 / 21 000 random inputs (the F18/F42/F43 sets r1 3 000, r2 6 000, m1 4 000 incl. the per-step removal order, plus 8 000 new inputs with 1–5 samples and 2–8 clusters), and `JavaStringHashMapOrder` reproduces `HashMap.keySet()` on 8 000 / 8 000 random key sets (1 404 with treeified bins). Insertion-order grouping would have differed on 5 / 3 000 static and 85 / 3 000 member-level inputs (locked: t02042, t02235, t00750, t01561 equal-score ties; t01396 error-score summation order). **Search reuse:** the repository suffix tree was evaluated and is **not** applicable — this unit performs no substring/pattern search; it is a numeric constraint-satisfaction tree build over CCF vectors.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Lineage precedence edge rule `u.CCF[i] ≥ v.CCF[i] − ε` and presence pattern `u=0 ⇒ v=0` [1] Eq. 2.
- Sum rule `Σ_children v.CCF[i] ≤ u.CCF[i] + ε` per node, per sample [1] Eq. 5; [2] sum condition.
- Rooted spanning tree over all clusters [1]; LICHeE constraint network, Gabow–Myers enumeration, error-score ranking, complete-network fallback [3].

**Deviations from LICHeE:**

- `ReconstructPhylogeny(CcfCluster…)` carries point CCFs only, so it uses the static margin ε (= LICHeE with zero-variance clusters); the per-cluster `1.96·sd/√n` margins are available through `ReconstructPhylogenyFromClusterSummaries` (F42).
- `ReconstructPhylogeny(CcfCluster…)` has no member counts, so every cluster is treated as robust — for robust clusters LICHeE's `fixNetwork` removes nothing, so this is LICHeE's own behaviour; the summary overloads run `fixNetwork` by default (F43) and report dropped clusters in `RemovedClusterIds` (`removeNonRobustClusters: false` keeps them).
- After a `fixNetwork` rebuild LICHeE orders groups by a `HashSet<SNVGroup>` (identity hash codes, `SNVGroup` overrides `equals` but not `hashCode`); the port keeps the pre-removal node order (= LICHeE with a constant identity hash). This is the only identity-hash dependence; it is not reproducible run to run in LICHeE itself (F43). LICHeE's SNV-level front end (SNV filtering, grouping SNVs by presence profile and clustering each group into sub-populations — `SNVGroup.setSubPopulations` size filtering / centroid collapsing) is outside ONCO-PHYLO-001, whose input is clusters (ALGORITHMS_CHECKLIST_V2 ONCO-PHYLO-001: "CCF clustering itself is ONCO-CCF-001"); since F42/F43 per-cluster summaries (member count, mean, sd) can be built from member CCFs with `CcfClusterSummary.FromMembers` and passed to `ReconstructPhylogenyFromClusterSummaries`.
- The group order emulates LICHeE's clusters-file path (profiles `put` in first-appearance order, never removed). On LICHeE's SNV-file path `tag2SNVs` also loses small groups and gains re-assigned ambiguous ones before `keySet()` is read, so its table capacity/order depends on that SNV-level history, which this cluster-level API does not have.

**Not implemented:**

- Probabilistic clone clustering — scope of the proposed unit ONCO-PYCLONE-001 (Dirichlet-process subclonal clustering; ALGORITHMS_CHECKLIST_V2.md); deterministic CCF clustering is ONCO-CCF-001. Posterior over trees (PhyloWGS, CITUP) is outside ONCO-PHYLO-001, a deterministic sum-rule method; **users should rely on:** dedicated external tools for posterior tree inference.
- CNA-aware multiplicity corrections to CCF inside this unit (it consumes already-corrected CCF clusters); **users should rely on:** `DeriveMultiplicity` (ONCO-ASCAT-001) + `EstimateCcf` with the local copy number (ONCO-CCF-001) upstream.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | ~~Profile order = first appearance (LICHeE: HashMap order)~~ | ~~Assumption~~ | Resolved by F45: Java `HashMap<String>` iteration emulated; 16 000 / 16 000 random inputs + 8 000 / 8 000 key sets identical to lichee.jar / java.util.HashMap | resolved | ASM-03 |
| 2 | ~~Default ε = 0~~ | ~~Assumption~~ | Resolved by F44: default = LICHeE's 0.1 | resolved | — |
| 5 | ~~Post-10⁸-cap state~~ | ~~Deviation~~ | Resolved by F44: no tree can be added after the cap in LICHeE either (equivalence + reduced-cap jar locks) | resolved | — |
| 3 | ~~No `fixNetwork` cluster removal~~ | ~~Deviation~~ | Resolved by F43 (summary overloads; `ReconstructPhylogeny` = all clusters robust, where `fixNetwork` is a no-op) | resolved | — |
| 4 | Group order after `fixNetwork` = pre-removal node order (LICHeE: identity-hash `HashSet`) | Assumption | equal to LICHeE with a constant identity hash; differs from a default JVM run on ~4 % of random inputs with removals | accepted | F43 |

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
