---
type: concept
title: "Tumor phylogeny reconstruction (clonal tree from CCF clusters — sum + lineage-precedence rules)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-PHYLO-001-Evidence.md
  - docs/algorithms/Oncology/Tumor_Phylogeny_Reconstruction.md
source_commit: abca521a486fd8b7aec2566f4f42e5dc27a99769
created: 2026-07-10
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-phylo-001-evidence
      evidence: "Test Unit ID: ONCO-PHYLO-001 ... Algorithm: Tumor Phylogeny Reconstruction — clonal tree from CCF clusters (sum rule + lineage-precedence rule)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:cancer-cell-fraction-clonal-clustering
      source: onco-phylo-001-evidence
      evidence: "The algorithm builds the clonal tree FROM CCF clusters — every test dataset is a table of (cluster, per-sample CCF) rows, the deconvoluted clusters ONCO-CCF-001 ClusterCCFValues produces; the reconstruction step ASCAT/CCF-clustering stops short of."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:clonal-subclonal-classification-ccf-posterior
      source: onco-phylo-001-evidence
      evidence: "Trunk identification = clusters on the path from root with CCF ≈ 1 in all samples (the clonal populations); branches = all non-trunk clusters (subclones) — the same clonal/subclonal distinction ONCO-CLONAL-001 makes per mutation, here made structurally over the tree."
      confidence: medium
      status: current
---

# Tumor phylogeny reconstruction (clonal tree from CCF clusters)

The **clonal-tree reconstruction layer** of the Oncology family: given per-cluster **cancer cell
fractions (CCF)** across one or more tumor samples, assemble the clusters into a **clonal-evolution
tree** — a rooted phylogeny whose root is the normal cell and whose edges encode
parent→child subclonal descent. This is the **reconstruction step** that CCF estimation and
clustering ([[cancer-cell-fraction-clonal-clustering]]) stop short of: it takes the CCF clusters and
decides **who descends from whom**. Validated under test unit **ONCO-PHYLO-001**; the
literature-traced record is [[onco-phylo-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade correctness
reference — [[scientific-rigor|research-grade]], **not for clinical or diagnostic use**.

This is a **constraint-satisfaction / perfect-phylogeny** tree builder (LICHeE, Popic 2015;
PICTograph, Zheng 2022), **not** a distance-based phylogenetics method — it does not compute a
distance matrix or run Neighbor-Joining / UPGMA (those operate on sequence divergence, not
cell-fraction ordering). The clonal tree is determined entirely by two per-sample CCF-ordering
constraints plus a presence pattern.

## 1. The two ordering constraints

An edge `(parent u → child v)` may exist only if, **for every sample i**, both hold (LICHeE Eq. 2 /
PICTograph lineage precedence):

- **Lineage precedence (ancestor ≥ descendant):** `u.CCFᵢ ≥ v.CCFᵢ − ε` — a mutation cannot have a
  higher CCF than its ancestor, since every cell carrying the descendant also carries the ancestor.
  Plus the **presence pattern**: if `u.CCFᵢ = 0` then `v.CCFᵢ = 0` — a descendant is present only in
  (a subset of) the samples where its ancestor is present (LICHeE constraint 1).

The tree as a whole must additionally satisfy, at **every node** and **every sample** (LICHeE Eq. 5 /
PICTograph "sum condition"):

- **Sum rule (pigeonhole generalization):** `Σ_{v : (u→v)} v.CCFᵢ ≤ u.CCFᵢ + ε` — the CCFs of all
  **children** of a node may not sum to more than the parent's CCF, because the disjoint subclones
  descending from a common parent must together fit inside the parent's cell population. This is the
  cell-fraction analogue of the pigeonhole principle applied per node. It is stated as an **inequality,
  not equality**, because not every true lineage branch need have been observed.

The sum rule is what forces **branching vs nesting**: two sibling subclones whose CCFs would overflow
the parent budget cannot both attach to it; at least one must instead nest deeper (become a
descendant of the other), converting a would-be sibling pair into a chain.

## 2. Trunk vs branch identification

- **Trunk** = the clusters on the path from the root that are present with CCF ≈ 1 across all samples
  — the **clonal** populations (the common predecessor shared by every cell). In the linear oracle
  below, Trunk = {A}.
- **Branches** = all non-trunk clusters — the **subclonal** lineages. Branches = {B, C}.

This is the same clonal/subclonal distinction that
[[clonal-subclonal-classification-ccf-posterior]] makes **per mutation** (probabilistically) and that
[[cancer-cell-fraction-clonal-clustering]] makes by **highest-centroid = clonal** — here it is made
**structurally** as a property of position in the reconstructed tree.

## 3. Deterministic tie-break — deepest valid ancestor (assumption)

Popic et al. note that **private / single-sample** clusters are **under-constrained**: multiple tree
nodes can validly serve as their ancestor, so a valid tree is not unique. To make the output
deterministic (an explicit **ASSUMPTION**, source-consistent), the unit attaches each cluster to its
**deepest valid ancestor** — the candidate parent with the smallest total CCF whose per-sample
sum-rule budget still admits the child — with **ties broken by ascending cluster id**. This selects
the most-recent common ancestor consistent with all cited constraints; it changes only **which single
valid tree is returned**, never **which trees are valid**.

## 4. Noise margin ε (assumption)

The cited sources relax both inequalities by a configurable ε (LICHeE ϵ; PICTograph ε₁ = 0.1 lineage,
ε₂ = 0.2 sum). Because this unit consumes **already-clustered CCF point estimates** (clustering, with
its noise model, is [[cancer-cell-fraction-clonal-clustering|ONCO-CCF-001]]), the default comparison
uses **ε = 0** (strict inequalities), exposed as an optional tolerance parameter so callers can supply
the source defaults. Setting ε > 0 only **widens admissibility** — it never flips a strictly-satisfied
relationship.

## Worked datasets (deterministic, hand-derived)

- **Linear (chain) evolution — single sample:** Normal(1.0), A(1.0), B(0.6), C(0.3) → edges
  Normal→A→B→C (each child ≤ parent, budget never exceeded). Trunk = {A}, Branches = {B, C}.
- **Branching evolution — two samples:** Normal(1,1), A(1,1), B(0.6, 0.0), C(0.0, 0.7) → edges
  Normal→A, A→B, A→C. A is the trunk parent of two **sibling** branches: sum rule under A holds per
  sample (s1 = 0.6, s2 = 0.7, both ≤ 1.0), and B, C cannot be ancestor/descendant of each other
  because each is **absent** in the other's sample (presence-pattern / lineage-precedence violated
  both directions).
- **Sum rule forces a chain instead of siblings — single sample:** Normal(1.0), A(1.0), B(0.6),
  C(0.6). B and C cannot **both** be children of A (0.6 + 0.6 = 1.2 > 1.0 → sum-rule violation). The
  deterministic rule attaches B under A first, then C must **nest under B** (0.6 ≤ 0.6, budget
  0.6 ≥ 0.6) → chain Normal→A→B→C. The sum rule converts a would-be sibling pair into a chain.

## Reconstructed-tree invariants (test oracles)

Two invariants must hold on **every** reconstructed tree, and are the primary correctness oracles:

1. **Ancestor CCF ≥ descendant CCF** on every edge, per sample (LICHeE Eq. 2).
2. **Per-node sum rule** — children CCF sum ≤ parent CCF, per sample (LICHeE Eq. 5 / PICTograph).

Boundary cases: **empty input** → tree with only the root, no trunk/branch mutations; **single
cluster** → child of root, and it is the trunk.

## Implementation surface (ONCO-PHYLO-001 spec)

Entry points on `OncologyAnalyzer` (`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`):

- `ReconstructPhylogeny(IReadOnlyList<CcfCluster> clusters, double tolerance = 0.0)` → `ClonalPhylogeny`
  with `RootId` (synthetic normal root, id below the minimum cluster id, CCF = 1 in every sample),
  `Clusters` (input order), `Edges` (one `ClonalEdge` parent→child per non-root cluster), `SampleCount`.
- `IdentifyTrunkMutations(ClonalPhylogeny)` → trunk (clonal) cluster ids from the root downward.
- `IdentifyBranchMutations(ClonalPhylogeny)` → the remaining (subclonal) ids in input order; trunk and
  branch **partition** the clusters (disjoint, union = all).

**Construction order:** process clusters by descending total CCF (ties: ascending id) so ancestors
precede descendants; initialise each node's per-sample sum-rule budget to its own CCF (root = 1); for
each cluster pick the deepest valid ancestor and **debit** that parent's per-sample budget by the
child's CCF.

**Contract / validation:** each `CcfCluster.CcfPerSample` must be same non-zero length with values in
[0,1] and unique ids; `tolerance` ≥ 0 and not NaN. Null `clusters`/CCF list → `ArgumentNullException`;
empty/ragged CCF, NaN/out-of-[0,1], or duplicate id → `ArgumentException`; negative/NaN `tolerance` →
`ArgumentOutOfRangeException`; empty cluster list → a root-only tree.

**Complexity:** `ReconstructPhylogeny` is **O(n²·k)** time, **O(n·k)** space (n clusters, k samples —
each cluster scans up to n placed candidates, each check O(k)). **Search reuse:** the repository suffix
tree was evaluated and rejected as inapplicable — this is a numeric constraint-satisfaction build over
CCF vectors, not substring/pattern search.

## Relationship to the neighbouring Oncology units

- **Upstream — CCF clustering** [[cancer-cell-fraction-clonal-clustering]] (ONCO-CCF-001) produces the
  CCF clusters this unit consumes. That unit **deconvolutes** CCFs into clones/subclones; this unit
  **orders** those clusters into a tree. `depends_on`.
- **Sibling — probabilistic clonal classification**
  [[clonal-subclonal-classification-ccf-posterior]] (ONCO-CLONAL-001) labels each mutation clonal vs
  subclonal on its own read-count uncertainty; this unit answers the structural version (trunk vs
  branch) via the tree topology. `relates_to`.
- **Downstream — heterogeneity metrics** [[intratumor-heterogeneity-metrics]] (ONCO-HETERO-001)
  reduces the reconstructed clonal structure to scalar ITH numbers (subclone count, Shannon diversity)
  — this unit supplies that structure.

## Scope and limitations

Sources are mutually consistent: Popic et al. (2015, *Genome Biology*, LICHeE) and Zheng et al.
(2022, *Bioinformatics*, PICTograph) state the **same** two constraints in VAF and CCF form
respectively, with PICTograph naming the sum rule the pigeonhole generalization. **Two flagged
assumptions**, both source-consistent rather than source-mandated: (1) the **deepest-valid-ancestor**
deterministic tie-break for under-constrained private clusters (does not change the valid-tree set),
and (2) the default **ε = 0** strict inequalities (the sources' ε only widens admissibility).
Distance-based sequence phylogenetics (Neighbor-Joining / UPGMA, and its bootstrap-support layer
[[phylogenetic-bootstrap-support]]) is a **different problem** and out of scope — that family orders
taxa by sequence divergence and resamples alignment columns, whereas this unit orders clusters by
per-sample cell-fraction inequalities and never builds a distance matrix.
**Not for clinical or diagnostic use.**
