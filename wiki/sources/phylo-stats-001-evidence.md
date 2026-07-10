---
type: source
title: "Evidence: PHYLO-STATS-001 (Tree Statistics — leaves, total tree length, tree depth/height)"
tags: [validation, phylogenetics]
doc_path: docs/Evidence/PHYLO-STATS-001-Evidence.md
sources:
  - docs/Evidence/PHYLO-STATS-001-Evidence.md
source_commit: 956d8f52e81160361eaf4673e2b2dedcc906ea08
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PHYLO-STATS-001

The validation-evidence artifact for test unit **PHYLO-STATS-001** — **Tree Statistics**: three
descriptive summaries read off a single (already-built) phylogenetic tree — **`GetLeaves`** (the
terminal/leaf nodes), **`CalculateTreeLength`** (the total tree length = sum of all branch lengths), and
**`GetTreeDepth`** (the tree height = number of edges on the longest root-to-leaf path). This is a
per-algorithm instance of the templated [[algorithm-validation-evidence|evidence artifact]] pattern; the
definitions, conventions, and edge cases are synthesized in the dedicated concept [[tree-statistics]].
See [[test-unit-registry]] for how units are tracked. Distinct from the comparison/query operations of
[[tree-comparison-metrics]] (PHYLO-COMP-001), which read *relationships between* taxa rather than
whole-tree descriptive summaries.

## What this file records

- **Online sources:** Wikipedia **Tree (graph theory)** and **Tree (abstract data type)** (leaf = vertex
  with no children / degree-1 external node; **depth** of a vertex = length of its root path, root depth 0;
  **height** of a vertex = length of the longest downward path to a leaf, height of tree = height of root;
  single-node tree has height and depth 0; **empty tree has height −1** by convention); **Biopython**
  `Bio.Phylo.BaseTree` (`get_terminals`, `count_terminals`, `is_terminal`, `total_branch_length` = sum of
  all branch lengths); **DendroPy** `Tree.length()` (sum of edge lengths, **edges with no length → 0**);
  and Wikipedia **Minimum evolution** (context: ME selects the branching pattern of smallest total branch
  length; Rzhetsky & Nei 1992) establishing "total tree length = sum of branch lengths" as the ME quantity.
- **Key definitions:** leaf = node with no children (`GetLeaves` returns exactly the terminal nodes);
  tree length = Σ branch lengths across all edges; tree depth = height = edge count on the longest
  root-to-deepest-leaf path. Height counts **edges**, not nodes.
- **Test datasets (worked by hand):** balanced 4-taxon `((A:1,B:1):1,(C:1,D:1):1)` → **4 leaves**,
  **length 6**, **height 2**; caterpillar/ladder `(A:1,(B:1,(C:1,D:1)))` (internal edges 0.5) → **4 leaves**,
  **length 5**, **height 3**; single-leaf `A` → 1 leaf, length = its branch length, **height 0**; null root
  → no leaves, length 0, **height −1**.
- **Corner cases (§ Documented Corner Cases):** single-node tree (root that is also a leaf) → height and
  depth 0; empty tree → height −1; undefined/default (0) branch length treated as 0 when summing length.

## Deviations and assumptions

No deviations from the literature. **One assumption (§Assumptions):** a `null` `PhyloNode` reference is the
repository's representation of "no tree", so it maps onto the graph-theory **empty-tree** convention —
`GetTreeDepth(null)` → −1, and `GetLeaves(null)` / `CalculateTreeLength(null)` → empty/zero. This is the
natural mapping of the empty-tree convention onto the nullable node model; no cited source addresses a C#
null reference directly. No source contradictions.
