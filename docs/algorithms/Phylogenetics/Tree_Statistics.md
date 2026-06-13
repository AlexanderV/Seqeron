# Tree Statistics

| Field | Value |
|-------|-------|
| Algorithm Group | Phylogenetics |
| Test Unit ID | PHYLO-STATS-001 |
| Related Projects | Seqeron.Genomics.Phylogenetics |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Tree Statistics provides three exact, deterministic descriptive measures over a (rooted, binary) phylogenetic tree: the set of leaves (terminal taxa), the total tree length (sum of all branch lengths), and the tree height/depth (number of edges on the longest root-to-leaf path). Each is computed by a single linear traversal of the tree and is used to summarize a tree produced by UPGMA, Neighbor-Joining, or Newick parsing. The measures are combinatorial/exact (no estimation), and follow the standard graph-theory and phylogenetics definitions used by reference libraries such as Biopython and DendroPy.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A phylogenetic tree is a rooted (here, binary) tree whose leaves are taxa/operational taxonomic units and whose edges carry branch lengths (e.g., expected substitutions per site or time). Three quantities recur in downstream analysis: the leaf set (which taxa the tree spans), the total branch length (the amount of evolution represented, minimized by the minimum-evolution criterion [5]), and the tree height (how many speciation/coalescent levels deep the tree is).

### 2.2 Core Model

- **Leaf:** a leaf is "a vertex with no children" [1]; equivalently a terminal/external node [2]. Biopython's `get_terminals` returns "all of this tree's terminal (leaf) nodes" [3].
- **Total tree length:** the sum of the branch lengths of every edge of the tree — DendroPy `Tree.length()` returns "the sum of edge lengths of self" [4]; Biopython `total_branch_length` "calculate[s] the sum of all the branch lengths in this tree" [3]. This is the quantity minimized by minimum evolution [5].
- **Tree height (depth):** "the length of the longest downward path to a leaf from [the root]" measured in edges; "the height of the root is the height of the tree" [2]. The root has depth 0 [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every node returned by `GetLeaves` has no children (`IsLeaf`). | Definition of a leaf [1][2]. |
| INV-02 | For an N-leaf fully-bifurcating tree, `GetLeaves` returns exactly N nodes. | Traversal visits every terminal node once [3]. |
| INV-03 | `CalculateTreeLength` = Σ over all nodes of `node.BranchLength`. | Sum-of-edge-lengths definition [3][4]. |
| INV-04 | `CalculateTreeLength ≥ 0` when all branch lengths are ≥ 0. | Sum of non-negative terms [4]. |
| INV-05 | `GetTreeDepth` = number of edges on the longest root→leaf path; a leaf-only tree = 0. | Height definition; single node has height 0 [1][2]. |
| INV-06 | `GetTreeDepth(null) = -1`; `GetLeaves(null)` empty; `CalculateTreeLength(null) = 0`. | Empty-tree convention (height −1) [1][2]; null = "no tree" (ASM-01). |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | This unit's `GetTreeDepth` | Biopython `depths()` |
|--------|---------------------------|----------------------|
| Quantity | Topological height in **edges** | Distance from root **by branch length** [3] |
| Use | Tree shape/balance | Patristic depth / ultrametricity |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| root | `PhylogeneticAnalyzer.PhyloNode` | required | Root of the (binary) tree | May be `null` (treated as empty tree) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `GetLeaves` | `IEnumerable<PhyloNode>` | Terminal nodes in left-to-right (pre-order) order |
| `CalculateTreeLength` | `double` | Sum of all branch lengths (same units as `BranchLength`) |
| `GetTreeDepth` | `int` | Height in edges to the deepest leaf; `-1` for a null tree |

### 3.3 Preconditions and Validation

No exceptions are thrown for tree shape. `null` root is accepted and maps to the empty-tree case: `GetLeaves` yields nothing, `CalculateTreeLength` returns `0`, `GetTreeDepth` returns `-1`. Branch lengths default to `0` and are summed as-is. The model is a binary tree (`Left`/`Right`); multifurcations are not representable (enforced at parse time elsewhere).

## 4. Algorithm

### 4.1 High-Level Steps

1. **GetLeaves:** recurse pre-order; yield a node iff it `IsLeaf`; on `null`, yield nothing.
2. **CalculateTreeLength:** return `0` for `null`; otherwise `root.BranchLength + length(Left) + length(Right)`.
3. **GetTreeDepth:** return `-1` for `null`, `0` for a leaf; otherwise `1 + max(depth(Left), depth(Right))`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| GetLeaves / CalculateTreeLength / GetTreeDepth | O(n) | O(h) | n = nodes, h = tree height (recursion stack) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PhylogeneticAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs)

- `PhylogeneticAnalyzer.GetLeaves(PhyloNode)`: yields terminal nodes via pre-order recursion.
- `PhylogeneticAnalyzer.CalculateTreeLength(PhyloNode)`: recursive sum of `BranchLength`.
- `PhylogeneticAnalyzer.GetTreeDepth(PhyloNode)`: recursive edge-count height; `EmptyTreeHeight = -1` constant for the null/empty case.

### 5.2 Current Behavior

`GetLeaves` is implemented with iterator (`yield`) pre-order recursion. `GetTreeDepth` uses the named constant `EmptyTreeHeight` (−1) for null. This is a pure tree-traversal unit, not a search/matching operation, so the repository suffix tree is **not applicable** (no substring/occurrence search is involved).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Leaf = node with no children [1][2][3].
- Total tree length = sum of all branch lengths [3][4].
- Tree height = longest root→leaf path in edges; leaf-only tree = 0 [1][2].
- Empty tree height = −1 [1][2].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Branch-length depth (`Biopython depths()`); **users should rely on:** `PatristicDistance` / `CalculateTreeLength` for branch-length quantities.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Null root → empty tree | Assumption | `GetTreeDepth(null) = -1` per convention | accepted | ASM-01; [1][2] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null root | leaves ∅, length 0, height −1 | Empty-tree convention [1][2] |
| Single leaf | one leaf, length = its `BranchLength`, height 0 | Single node has height 0 [1][2] |
| All branch lengths default (0) | length 0 | Sum of zeros [4] |

### 6.2 Limitations

Binary tree model only; multifurcating trees are not represented. `GetTreeDepth` is topological (edges), not patristic (branch-length) depth.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var t = PhylogeneticAnalyzer.ParseNewick("((A:1,B:1):1,(C:1,D:1):1);");
int leaves = PhylogeneticAnalyzer.GetLeaves(t).Count();      // 4
double len = PhylogeneticAnalyzer.CalculateTreeLength(t);    // 6.0
int height = PhylogeneticAnalyzer.GetTreeDepth(t);           // 2
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PhylogeneticAnalyzer_TreeStatistics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_TreeStatistics_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [PHYLO-STATS-001-Evidence.md](../../../docs/Evidence/PHYLO-STATS-001-Evidence.md)
- Related algorithms: [Tree_Construction](Tree_Construction.md), [Newick_Format](Newick_Format.md)

## 8. References

1. Wikipedia. 2026. Tree (graph theory). https://en.wikipedia.org/wiki/Tree_(graph_theory)
2. Wikipedia. 2026. Tree (abstract data type). https://en.wikipedia.org/wiki/Tree_(abstract_data_type)
3. Cock PJA et al. Biopython `Bio.Phylo.BaseTree` (get_terminals, count_terminals, is_terminal, total_branch_length). https://biopython.org/docs/latest/api/Bio.Phylo.BaseTree.html
4. Sukumaran J, Holder MT. DendroPy `Tree.length()`. https://dendropy.org/library/treemodel.html
5. Rzhetsky A, Nei M. 1992. A simple method for estimating and testing minimum-evolution trees. (via) Wikipedia. Minimum evolution. https://en.wikipedia.org/wiki/Minimum_evolution
