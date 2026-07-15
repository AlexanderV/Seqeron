---
type: concept
title: "Tree statistics (leaf count, total tree length, tree depth/height)"
tags: [phylogenetics, algorithm]
mcp_tools:
  - tree_depth
sources:
  - docs/Evidence/PHYLO-STATS-001-Evidence.md
source_commit: 956d8f52e81160361eaf4673e2b2dedcc906ea08
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: phylo-stats-001-evidence
      evidence: "Test Unit ID: PHYLO-STATS-001 ... Algorithm: Tree Statistics (leaves, total tree length, tree height/depth)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:distance-based-tree-construction
      source: phylo-stats-001-evidence
      evidence: "operate on the same rooted binary PhyloNode that the UPGMA/NJ builder produces — the statistics summarize the tree the construction step (PHYLO-TREE-001) emits"
      confidence: high
      status: current
---

# Tree statistics (leaf count, total tree length, tree depth/height)

A phylogenetics-family (`PHYLO-*`) unit, **PHYLO-STATS-001** — three **descriptive summaries read off a
single, already-built phylogenetic tree**: **`GetLeaves`** (its terminal nodes), **`CalculateTreeLength`**
(total tree length = sum of every branch length), and **`GetTreeDepth`** (tree height = edges on the
longest root-to-leaf path). These are whole-tree *summary statistics*, not inference and not comparison —
genuinely separate from [[tree-comparison-metrics]] (PHYLO-COMP-001), which computes *relationships between*
taxa (RF distance, MRCA, patristic distance). They operate on the same rooted binary `PhyloNode`
(`Left`/`Right` children) that the [[distance-based-tree-construction]] step (PHYLO-TREE-001, UPGMA/NJ)
produces from an [[evolutionary-distance-matrix]] (PHYLO-DIST-001) and that
[[phylogenetic-bootstrap-support]] resamples. Validated under test unit
**PHYLO-STATS-001**; the literature-traced record is [[phylo-stats-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the evidence-artifact pattern.
Research-grade correctness reference ([[scientific-rigor|research-grade]]), not for clinical use.

## 1. Leaves (`GetLeaves`)

The **terminal nodes** — nodes with no children (graph theory: a *leaf* is a vertex with no children /
a degree-1 external vertex). `GetLeaves` returns exactly those nodes; a leaf count is
`GetLeaves(...).Count()`, matching Biopython's `get_terminals` / `count_terminals`.

## 2. Total tree length (`CalculateTreeLength`)

The **sum of all branch lengths** across every edge in the tree (Biopython `total_branch_length`,
DendroPy `Tree.length()`). This is the quantity the **minimum-evolution** criterion minimizes (Rzhetsky &
Nei 1992). Convention: an **undefined or default (0) branch length counts as 0** in the sum (DendroPy:
"edges with no lengths defined will be considered to have a length of 0").

## 3. Tree depth / height (`GetTreeDepth`)

The **height** of the tree — the **number of edges** on the longest downward (root-to-leaf) path. Height
counts edges, not nodes: in a rooted tree the height of a vertex is "the length of the longest downward
path to a leaf," and the tree's height is the height of its root. This unit names the whole-tree quantity
*depth*, but it is the graph-theory **height** value (they coincide at the deepest leaf: the depth of the
deepest leaf equals the height of the root).

## Conventions, invariants, and test oracles

| Tree | Leaves | Total length | Depth/height |
|------|--------|--------------|--------------|
| Balanced 4-taxon `((A:1,B:1):1,(C:1,D:1):1)` | 4 | 6 | 2 |
| Caterpillar `(A:1,(B:1,(C:1,D:1)))` (internal 0.5) | 4 | 5 | 3 |
| Single leaf `A` | 1 | its branch length | 0 |
| Null root (empty tree) | ∅ | 0 | −1 |

- **Single-node tree** (a root that is also a leaf) → height and depth **0**.
- **Empty tree** → height **−1** by convention.
- **Height counts edges** — the deepest-leaf path is measured in edges, so a two-level balanced tree has
  height 2, not 3.

## Empty-tree ↔ null mapping (one assumption)

The only assumption in the source: a **`null` `PhyloNode`** is the repository's representation of "no tree",
so it maps onto the graph-theory **empty-tree** convention — `GetTreeDepth(null)` → **−1**, and
`GetLeaves(null)` / `CalculateTreeLength(null)` → empty / **0**. This is the natural mapping of the
empty-tree convention onto the nullable node model; no cited source addresses a C# null reference directly.
No deviations from the literature and no source contradictions.

## Relationship to the rest of the PHYLO family

These statistics *summarize* a tree; [[tree-comparison-metrics]] *compares or queries* two taxa or two
trees; [[phylogenetic-bootstrap-support]] attaches *confidence* to a tree by resampling; and
[[evolutionary-distance-matrix]] supplies the *pairwise-distance substrate* the UPGMA/NJ builder turns into
the `PhyloNode` tree. The `CalculateTreeLength` sum-over-all-branches is the whole-tree analogue of the
per-path branch-length sum that [[tree-comparison-metrics|patristic distance]] computes between one taxon
pair. Trees are serialized to/from **Newick** text by the family's I/O layer, PHYLO-NEWICK-001
([[phylo-newick-001-evidence]]).

## Reference tools

Definitions trace to the Wikipedia articles on **Tree (graph theory)** and **Tree (abstract data type)**
(leaf / depth / height definitions, single-node-tree height 0, empty-tree height −1), **Biopython**
`Bio.Phylo.BaseTree` (`get_terminals`, `count_terminals`, `is_terminal`, `total_branch_length`),
**DendroPy** `Tree.length()` (sum of edge lengths, missing length → 0), and **Rzhetsky & Nei (1992)** via
Wikipedia **Minimum evolution** (total branch length as the ME criterion). No source contradictions.
