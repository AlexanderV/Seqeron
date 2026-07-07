# Tumor Phylogeny Reconstruction

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Clonal evolution |
| Test Unit ID | ONCO-PHYLO-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-15 |

## 1. Overview

Reconstructs a rooted clonal (tumor) phylogeny — the ancestor/descendant ordering of subclones — from cancer cell fraction (CCF) clusters. Each cluster is a candidate clone with a CCF in every sequenced sample; the clusters themselves are produced upstream (ONCO-CCF-001). The method builds the tree by applying two deterministic lineage constraints from the multi-sample perfect-phylogeny model — lineage precedence (an ancestor's CCF is at least its descendant's) and the sum rule (children CCFs may not exceed the parent's) [1][2]. It is a constraint-satisfaction / greedy construction, not a probabilistic Bayesian clustering method (PyClone/PhyloWGS/CITUP); among the trees that satisfy the cited constraints it returns one deterministic tree via a documented tie-break.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A tumor evolves clonally: subclones arise by acquiring additional somatic mutations on top of an ancestral clone. Bulk/multi-region sequencing yields, for each mutation cluster, a cancer cell fraction (the fraction of tumor cells carrying it) in each sample. Under the infinite-sites / perfect-phylogeny assumption (a mutation arises once and is never lost), these CCFs constrain which clones can be ancestors of which others [1].

### 2.2 Core Model

Let cluster `v` have CCF `v.CCF[i]` in sample `i`, and let `ε` be a noise margin. The clonal tree `T` (rooted at a synthetic normal node with CCF = 1 in all samples) must satisfy:

- **Lineage precedence (ancestor ≥ descendant), Eq. 2 [1]:** an edge `u→v` is admissible only if for every sample `i`, `u.CCF[i] ≥ v.CCF[i] − ε`, and (presence pattern) `u.CCF[i] = 0 ⇒ v.CCF[i] = 0`. Equivalently, "the CCF of any mutation cannot exceed the CCF of its ancestor" [2].
- **Sum rule, Eq. 5 [1]:** for every node `u` and every sample `i`, `Σ_{v : (u→v)∈T} v.CCF[i] ≤ u.CCF[i] + ε`. Equivalently, "the CCF of an ancestral clone must be greater than or equal to the sum of CCFs of its descendants" [2]. This is the cell-fraction analogue of the pigeonhole principle: disjoint sibling subclones cannot together occupy more cells than their common parent [1].

These constraints define a *set* of valid spanning trees; they do not uniquely determine one (private/single-sample clusters are under-constrained) [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Infinite sites / perfect phylogeny: each mutation arises once and is never lost [1]. | Convergent or back-mutation breaks the ancestor ≥ descendant ordering; reconstructed edges may be wrong. |
| ASM-02 | CCFs are correctly estimated and clustered upstream (ONCO-CCF-001). | Mis-clustered CCFs propagate to wrong edges; this unit does not re-estimate them. |
| ASM-03 | Deepest-valid-ancestor + ascending-id tie-break selects among valid trees (Evidence Assumption 1). | A different (still valid) tree could be returned under another tie-break; topology among ambiguous private clusters is not unique. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every edge `u→v` satisfies `u.CCF[i] ≥ v.CCF[i] − ε` for all samples `i`. | Edges are only added when `SatisfiesLineagePrecedence` returns true [1] Eq. 2. |
| INV-02 | For every node, per sample, children CCF sum ≤ parent CCF + ε. | Each child debits the parent's budget; admission requires `FitsSumRule` [1] Eq. 5. |
| INV-03 | The result is a single rooted tree: every cluster has exactly one parent, no cycles. | Each cluster is attached exactly once to an already-placed candidate (the root or an earlier node) [1]. |
| INV-04 | Trunk and Branch partition the clusters (disjoint, union = all). | `IdentifyBranchMutations` = clusters ∉ trunk set. |
| INV-05 | Deterministic: identical input ⇒ identical edge set. | Processing order and parent choice are total orders (total CCF desc, then id; deepest valid then id). |

### 2.5 Comparison with Related Methods

| Aspect | This method | PyClone / PhyloWGS / CITUP |
|--------|-------------|----------------------------|
| Inference | Deterministic constraint satisfaction over given CCF clusters | Probabilistic (MCMC / ILP) joint clustering + tree |
| Clustering | Consumes clusters (ONCO-CCF-001) | Performs clustering jointly |
| Output | One tree (documented tie-break) | Posterior over trees / optimal ILP tree |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `clusters` | `IReadOnlyList<CcfCluster>` | required | CCF clusters to place | each `CcfPerSample` same non-zero length; values in [0,1]; unique ids |
| `tolerance` | `double` | `0.0` | Noise margin ε for Eq. 2 and Eq. 5 | ≥ 0, not NaN |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `RootId` | `int` | Synthetic normal root id (distinct from all cluster ids) |
| `Clusters` | `IReadOnlyList<CcfCluster>` | Input clusters in input order |
| `Edges` | `IReadOnlyList<ClonalEdge>` | parent→child edges, one per non-root cluster |
| `SampleCount` | `int` | Samples per cluster |

`IdentifyTrunkMutations` returns trunk cluster ids (root downward); `IdentifyBranchMutations` returns the remaining (subclonal) ids in input order.

### 3.3 Preconditions and Validation

Null `clusters` or null cluster CCF list → `ArgumentNullException`. Empty/ragged CCF lists, NaN/out-of-[0,1] CCF, or duplicate ids → `ArgumentException`. Negative/NaN `tolerance` → `ArgumentOutOfRangeException`. Empty cluster list → a root-only tree (no edges, no trunk/branch). CCF indexing is 0-based per sample; ids are caller-assigned integers.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; pick a synthetic root id below the minimum cluster id; the root has CCF = 1 in every sample.
2. Initialize each node's per-sample sum-rule budget to its own CCF (root budget = 1).
3. Process clusters in descending total CCF (ties: ascending id) so ancestors precede descendants.
4. For each cluster, among the root and all already-placed clusters, pick the *deepest valid ancestor*: smallest total CCF that satisfies lineage precedence (Eq. 2) and still has per-sample budget for the child (Eq. 5); ties by ascending id.
5. Add the edge and debit the chosen parent's per-sample budget by the child's CCF.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Lineage precedence `u.CCF[i] ≥ v.CCF[i] − ε` and presence `u.CCF[i]=0 ⇒ v.CCF[i]=0` [1] Eq. 2.
- Sum rule on remaining budget `budget_u[i] ≥ v.CCF[i] − ε` [1] Eq. 5.
- Default `ε = 0` (strict); source defaults are ϵ (LICHeE) / ε₁=0.1, ε₂=0.2 (PICTograph) [1][2], exposed via `tolerance`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ReconstructPhylogeny` | O(n² · k) | O(n · k) | n clusters, k samples; each of n clusters scans up to n candidates, each candidate check is O(k). Matches the checklist O(n²·k). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ReconstructPhylogeny(IReadOnlyList<CcfCluster>, double)`: builds the rooted clonal tree.
- `OncologyAnalyzer.IdentifyTrunkMutations(ClonalPhylogeny)`: clonal/trunk clusters (root→first branch point).
- `OncologyAnalyzer.IdentifyBranchMutations(ClonalPhylogeny)`: subclonal/branch clusters (the rest).

### 5.2 Current Behavior

Single-sample CCFs typically yield a caterpillar (chain) because any child with CCF ≤ parent fits the parent's full budget; genuine branching arises when (a) the per-sample presence pattern forbids nesting (a cluster absent in a sample cannot descend from one present there) or (b) the sum rule exhausts a parent's budget. The synthetic root guarantees every cluster has a parent (spanning tree). **Search reuse:** the repository suffix tree was evaluated and is **not** applicable — this unit performs no substring/pattern search; it is a numeric constraint-satisfaction tree build over CCF vectors.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Lineage precedence edge rule `u.CCF[i] ≥ v.CCF[i] − ε` and presence pattern `u=0 ⇒ v=0` [1] Eq. 2.
- Sum rule `Σ_children v.CCF[i] ≤ u.CCF[i] + ε` per node, per sample [1] Eq. 5; [2] sum condition.
- Rooted spanning tree over all clusters [1].

**Intentionally simplified:**

- Tree selection: among the valid set, one tree is returned via a deterministic deepest-valid-ancestor / id tie-break; **consequence:** for under-constrained private clusters the exact parent may differ from a probabilistic method's MAP tree, though it always satisfies INV-01/INV-02.
- Default `ε = 0`; **consequence:** stricter than the source defaults (ϵ; ε₁=0.1, ε₂=0.2) unless the caller passes `tolerance`.

**Not implemented:**

- Probabilistic clone clustering / posterior over trees (PyClone, PhyloWGS, CITUP); **users should rely on:** ONCO-CCF-001 for CCF clustering and dedicated external tools for posterior tree inference.
- CNA-aware multiplicity corrections to CCF; **users should rely on:** ONCO-CCF-001 / ONCO-CNA-001.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Deepest-valid-ancestor tie-break | Assumption | Determines the single returned tree among valid ones | accepted | ASM-03; Evidence Assumption 1 |
| 2 | Default ε = 0 | Assumption | Stricter admissibility than source defaults | accepted | configurable via `tolerance` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty cluster list | root-only tree; trunk={}, branches={} | nothing to place |
| Single cluster | root→cluster; trunk={cluster}; branches={} | only the trunk exists |
| Two equal-CCF clusters under one parent | only one nests as child; the other chains below it | sum rule (Eq. 5) [1] |
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

**Numerical walk-through:** Process order by total CCF: cluster 1 (2.0), 3 (0.7), 2 (0.6). Cluster 1 attaches to root (1≥1). Cluster 3 ([0,0.7]) is valid under both root and 1; deepest valid = 1 → edge 1→3. Cluster 2 ([0.6,0]) cannot descend from 3 (3 has 0 in sample 1, but 2 has 0.6 → presence violated) → attaches to 1. Sum rule under 1: s1 = 0.6+0 = 0.6 ≤ 1; s2 = 0+0.7 = 0.7 ≤ 1. Result: branching tree with trunk {1}.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_ReconstructPhylogeny_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_ReconstructPhylogeny_Tests.cs) — covers `INV-01`, `INV-02`, `INV-05`
- Evidence: [ONCO-PHYLO-001-Evidence.md](../../../docs/Evidence/ONCO-PHYLO-001-Evidence.md)
- Related algorithms: [Clonal_Subclonal_Classification](./Clonal_Subclonal_Classification.md)

## 8. References

1. Popic V, Salari R, Hajirasouliha I, Kashef-Haghighi D, West RB, Batzoglou S. 2015. Fast and scalable inference of multi-sample cancer lineages. *Genome Biology* 16:91. https://doi.org/10.1186/s13059-015-0647-8
2. Zheng L, Dang H, Niknafs N, et al. 2022. Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors (PICTograph). *Bioinformatics* 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac367
