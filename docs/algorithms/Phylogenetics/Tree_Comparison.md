# Tree Comparison

| Field | Value |
|-------|-------|
| Algorithm Group | Phylogenetics |
| Test Unit ID | PHYLO-COMP-001 |
| Related Projects | Seqeron.Genomics; Seqeron.Genomics.Tests |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

This document covers three complementary tree-comparison routines implemented by the repository: Robinson-Foulds distance for rooted topology comparison, most recent common ancestor (MRCA) lookup, and patristic distance for branch-length path sums. All three operate on rooted binary `PhyloNode` trees and are deterministic for a fixed tree structure and set of leaf names. Robinson-Foulds distance is reported as a raw rooted distance, while MRCA and patristic distance depend on exact taxon-name matches at leaf nodes. These routines analyze already-constructed trees; tree construction itself is documented separately in [Tree_Construction.md](Tree_Construction.md).

## 2. Scientific / Formal Basis

> A = Robinson-Foulds Distance, B = Most Recent Common Ancestor, C = Patristic Distance

### 2.A Robinson-Foulds Distance

#### Domain Context

Robinson-Foulds (RF) distance measures how much two phylogenetic tree topologies differ by comparing the bipartitions or clades induced by their internal edges (References 1, 4). In rooted trees, the same idea is commonly expressed in terms of rooted clades or clusters rather than unrooted splits (References 1, 4).

#### Core Model

For rooted trees $T_1$ and $T_2$, RF distance is the size of the symmetric difference between the clade sets induced by the two trees (References 1, 4).

$$
RF(T_1, T_2) = |S_1 \triangle S_2| = |S_1 \setminus S_2| + |S_2 \setminus S_1|
$$

Here, $S_1$ and $S_2$ are the clade or split sets extracted from the two trees. In this repository's rooted interpretation, only non-trivial clades are relevant: leaves are excluded because they contain a single taxon, and the full-tree root is excluded because it contains all taxa.

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-RF-01 | $RF(T, T) = 0$ | Identity is a defining metric property (References 1, 4). |
| INV-RF-02 | $RF(T_1, T_2) = RF(T_2, T_1)$ | Symmetric difference is symmetric (References 1, 4). |
| INV-RF-03 | $RF(T_1, T_2) \ge 0$ | Set-cardinality counts are non-negative (References 1, 4). |
| INV-RF-04 | RF is even for binary trees with the same leaf set | Rooted binary trees contribute paired clade differences in the symmetric difference (References 1, 4). |
| INV-RF-05 | For rooted binary trees with $n$ taxa, $RF(T_1, T_2) \le 2(n-2)$ | Each rooted binary tree has $n-2$ non-trivial clades, so the largest possible symmetric difference is twice that count (References 1, 4). |

### 2.B Most Recent Common Ancestor

#### Domain Context

In a rooted phylogenetic tree, each internal node represents the common ancestor of the taxa descending from that node (Reference 2). The MRCA of two taxa is the deepest node that is an ancestor of both taxa (Reference 2).

#### Core Model

For taxa $x$ and $y$ in a rooted tree, $MRCA(x, y)$ is the deepest shared ancestor on the root-to-leaf paths for $x$ and $y$ (Reference 2). When the same taxon is queried twice, the taxon itself is the MRCA (Reference 2).

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-MRCA-01 | $MRCA(x, x)$ is the node for $x$ | A taxon is trivially its own most recent common ancestor (Reference 2). |
| INV-MRCA-02 | $MRCA(x, y) = MRCA(y, x)$ | The common-ancestor relation is symmetric in the queried taxa (Reference 2). |
| INV-MRCA-03 | The MRCA is an ancestor of both queried taxa | This is the defining property of an MRCA (Reference 2). |
| INV-MRCA-04 | The MRCA is unique for a fixed rooted tree and taxon pair | A rooted tree has a unique deepest shared ancestor for a given pair of descendants (Reference 2). |

### 2.C Patristic Distance

#### Domain Context

Patristic distance measures the total branch length along the unique tree path connecting two taxa through their MRCA (Reference 3). It is a path-based quantity, so its biological interpretation depends on the meaning of the branch lengths stored on the tree (Reference 3).

#### Core Model

If $MRCA(x, y)$ is the most recent common ancestor of taxa $x$ and $y$, then patristic distance is the sum of the branch lengths from that MRCA to each taxon (Reference 3).

$$
PD(x, y) = d(x, MRCA(x, y)) + d(y, MRCA(x, y))
$$

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-PD-01 | Stored branch lengths are meaningful edge weights | The returned value is still a path sum, but it may not be biologically interpretable as divergence (Reference 3). |
| ASM-PD-02 | Interpreting branch length as elapsed time requires an ultrametric or molecular-clock context | The path sum should not be read as elapsed time unless the upstream tree construction justifies that interpretation (Reference 3). |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-PD-01 | $PD(x, x) = 0$ | No path needs to be traversed when the same taxon is queried twice (Reference 3). |
| INV-PD-02 | $PD(x, y) = PD(y, x)$ | The underlying tree path is undirected even when the tree itself is rooted (Reference 3). |
| INV-PD-03 | $PD(x, y) \ge 0$ when branch lengths are non-negative | The distance is a sum of path lengths (Reference 3). |
| INV-PD-04 | $PD(x, y) = d(x, MRCA) + d(y, MRCA)$ | This is the defining decomposition of patristic distance (Reference 3). |
| INV-PD-05 | Patristic distance satisfies triangle inequality only when stored branch lengths define a non-negative tree metric | Tree-path distance is a metric on a fixed weighted tree only under non-negative edge weights (References 2, 3). |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| [RF] `tree1` | `PhyloNode` | required | First rooted tree to compare | Expected to be a rooted binary `PhyloNode` tree; branch lengths are ignored by RF comparison. |
| [RF] `tree2` | `PhyloNode` | required | Second rooted tree to compare | Expected to be a rooted binary `PhyloNode` tree over the same conceptual taxon namespace as `tree1`. |
| [MRCA] `root` | `PhyloNode` | required | Root of the tree to search | Taxa are matched by exact leaf `Name` values. |
| [MRCA] `taxon1` | `string` | required | First queried taxon | Compared by exact string match against leaf names. |
| [MRCA] `taxon2` | `string` | required | Second queried taxon | Compared by exact string match against leaf names. |
| [PD] `root` | `PhyloNode` | required | Root of the weighted tree to search | The traversed child nodes must carry the `BranchLength` values to be summed. |
| [PD] `taxon1` | `string` | required | First queried taxon | Compared by exact string match against leaf names. |
| [PD] `taxon2` | `string` | required | Second queried taxon | Compared by exact string match against leaf names. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| [RF] Return value | `int` | Raw rooted Robinson-Foulds distance, computed as the symmetric-difference count of non-trivial clades. |
| [MRCA] Return value | `PhyloNode?` | The matching leaf when the same taxon is queried twice, the deepest common ancestor when both taxa exist, or `null` when no MRCA can be resolved. |
| [PD] Return value | `double` | Sum of branch lengths from the MRCA to each queried taxon, or `double.NaN` when the MRCA cannot be resolved. |

### 3.3 Preconditions and Validation

- [RF] `RobinsonFouldsDistance` compares rooted tree topology only; it does not use branch lengths.
- [RF] The implementation derives clades from internal subtrees with more than one taxon and fewer than all taxa.
- [MRCA] `FindMRCA` returns `null` immediately when `root` is `null`.
- [MRCA] Taxa are matched by exact `PhyloNode.Name` comparison at leaf nodes.
- [MRCA] When `taxon1 != taxon2`, a leaf result is interpreted as "only one taxon found", and the public method returns `null`.
- [PD] `PatristicDistance` first resolves the MRCA with `FindMRCA`; if no MRCA is found, it returns `double.NaN`.
- [PD] Path lengths are accumulated from the traversed child node's `BranchLength` values.

## 4. Algorithm

### 4.A Robinson-Foulds Distance

#### High-Level Steps

1. Traverse the first rooted tree and collect every non-trivial clade.
2. Traverse the second rooted tree and collect every non-trivial clade.
3. Compute the symmetric difference of the two clade sets.
4. Return the total number of clades present in exactly one of the two sets.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `RobinsonFouldsDistance` | `O(n^2 log n)` worst case | `O(n^2)` worst case | The current implementation repeatedly materializes and sorts subtree taxon lists while building clade keys, so it is not linear in the number of taxa or visited nodes. |

### 4.B Most Recent Common Ancestor

#### High-Level Steps

1. Recursively descend the left and right subtrees.
2. If the current node is a matching leaf, return that leaf.
3. If both subtrees report a match, return the current node as the MRCA.
4. Otherwise, propagate the non-null subtree result upward.
5. After recursion, reject a leaf result when two distinct taxa were requested, because that means only one taxon was found.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindMRCA` | `O(n)` | `O(h)` | The traversal may visit the full tree; stack usage follows the tree height `h`. |

### 4.C Patristic Distance

#### High-Level Steps

1. Resolve the MRCA of the two queried taxa.
2. Starting at that MRCA, recursively accumulate branch lengths to the first taxon.
3. Starting at that MRCA, recursively accumulate branch lengths to the second taxon.
4. Return the sum of the two path lengths.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PatristicDistance` | `O(n)` | `O(h)` | MRCA resolution dominates the cost; each path sum then follows at most the tree height `h`. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PhylogeneticAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs)

- `PhylogeneticAnalyzer.PhyloNode`: binary tree node type with `Name`, `BranchLength`, `Left`, `Right`, `IsLeaf`, and `Taxa`.
- `PhylogeneticAnalyzer.RobinsonFouldsDistance(PhyloNode, PhyloNode)`: computes raw rooted RF distance from the symmetric difference of clade sets built by `GetClades` and `CollectClades`.
- `PhylogeneticAnalyzer.FindMRCA(PhyloNode, string, string)`: resolves the MRCA via `FindMRCAInternal` and applies a public post-check for missing taxa.
- `PhylogeneticAnalyzer.PatristicDistance(PhyloNode, string, string)`: finds the MRCA and sums both MRCA-to-taxon paths through `DistanceToTaxon`.

### 5.2 Current Behavior

- [RF] The repository implements rooted RF through non-trivial clade comparison, not through branch lengths or a normalized RF variant. Each clade is materialized as a sorted taxon-name string, and leaves plus the full-tree root are excluded from the comparison set.
- [MRCA] The recursive helper returns matching leaves upward; the public method then converts a leaf result for two distinct taxa into `null`, which prevents a single found taxon from masquerading as a common ancestor.
- [PD] `PatristicDistance` delegates missing-taxon handling to `FindMRCA`, then adds the traversed child `BranchLength` values on the two paths descending from the MRCA without enforcing non-negative edge weights.

### 5.3 Conformance to Theory / Spec

#### 5.3.A Robinson-Foulds Distance

**Implemented (verbatim from the cited theory/spec):**

- Raw rooted RF distance is returned as the symmetric-difference count of non-trivial clade sets.
- Identical trees produce `0`, and the tests pin exact rooted values of `2` for differing three-taxa topologies and `4` for maximally different four-taxa examples.
- The implementation treats RF as a rooted clade comparison, which matches the rooted interpretation documented for RF in the evidence and test specification.

**Intentionally simplified:**

- (none)

**Not implemented:**

- Unrooted RF variants and multifurcating-tree handling; **users should rely on:** no current in-repo alternative.
- Normalized or generalized RF variants, including information-theoretic extensions (Reference 5); **users should rely on:** external tooling or a separate metric implementation.

#### 5.3.B Most Recent Common Ancestor

**Implemented (verbatim from the cited theory/spec):**

- Recursive deepest-common-ancestor search on a rooted binary tree.
- Querying the same taxon twice returns that taxon's leaf node.
- Cross-clade pairs return the root when that root is the deepest shared ancestor.
- Missing taxa return `null` rather than a partial match.

**Intentionally simplified:**

- (none)

**Not implemented:**

- Multifurcating-tree support beyond the repository's binary `PhyloNode` model; **users should rely on:** no current in-repo alternative.

#### 5.3.C Patristic Distance

**Implemented (verbatim from the cited theory/spec):**

- Patristic distance is computed as the sum of the two MRCA-to-taxon path lengths.
- Querying the same taxon twice returns `0.0`.
- Exact rooted examples in the tests verify documented path sums such as `PD(A,B) = 1.0`, `PD(A,C) = 5.0`, and `PD(C,D) = 2.0` on the four-taxa reference tree.

**Intentionally simplified:**

- Biological interpretation of branch lengths is not validated; **consequence:** the method always returns the stored path sum, but that value is only as meaningful as the caller-supplied `BranchLength` data.

**Not implemented:**

- Automatic branch-length calibration or alternate distance conventions for trees without meaningful edge weights; **users should rely on:** upstream tree construction or external calibration.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | [RF] Rooted binary-tree scope | Assumption | RF values describe rooted `PhyloNode` trees only; unrooted or multifurcating comparisons are out of scope. | accepted | `PhyloNode` exposes `Left` and `Right` children only, and the evidence doc resolves RF scope as rooted and binary. |
| 2 | [RF] Raw RF count | Assumption | Results differ from normalized-RF software unless callers rescale externally. | accepted | The implementation returns the raw symmetric-difference count with no normalization step. |
| 3 | [PD] Caller-supplied branch-length semantics | Assumption | The numerical result is always a path sum, but its biological interpretation depends on how the tree was built. | accepted | `PatristicDistance` sums stored edge lengths and does not infer time calibration or confidence. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| [RF] A rooted tree with no non-trivial clades | RF contributes `0` clades from that tree | Leaves and the full-tree root are excluded from RF clade sets, so there is nothing to compare. |
| [MRCA] `root` is `null` | `FindMRCA` returns `null` | The public method guards `root == null` before recursion. |
| [MRCA] One or both taxa are missing | `FindMRCA` returns `null` | Tests cover both one-missing and both-missing cases, and the public post-check rejects partial matches. |
| [PD] Same taxon queried twice | `PatristicDistance` returns `0.0` | The MRCA is the queried leaf itself, so both path components are zero. |
| [PD] One or both taxa are missing | `PatristicDistance` returns `double.NaN` | The method returns `NaN` when `FindMRCA` cannot resolve an ancestor. |
| [PD] Traversed branch lengths are zero | The returned distance includes those zero-valued edges and can therefore be `0.0` even for different taxa | Patristic distance is defined as a path-length sum, so zero-valued edges contribute nothing. |

### 6.2 Limitations

- The repository models trees as rooted binary `PhyloNode` structures; none of the three routines handles multifurcating trees.
- `RobinsonFouldsDistance` reports a raw rooted clade distance only; there is no normalized or unrooted RF mode.
- MRCA and patristic distance compare taxa by exact leaf-name equality and do not normalize labels.
- Patristic distance inherits whatever branch-length semantics the caller or upstream tree-construction routine supplied.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PhylogeneticAnalyzer_TreeComparison_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_TreeComparison_Tests.cs) - covers `INV-RF-01` through `INV-RF-05`, `INV-MRCA-01` through `INV-MRCA-04`, and `INV-PD-01` through `INV-PD-05`.
- Test specification: [PHYLO-COMP-001.md](../../../tests/TestSpecs/PHYLO-COMP-001.md) - records the canonical scope, expected complexities, and the RF/MRCA implementation fixes adopted by the current code.
- Evidence: [PHYLO-COMP-001-Evidence.md](../../Evidence/PHYLO-COMP-001-Evidence.md) - captures the external references and resolved scope decisions for rooted binary trees.
- Related algorithms: [Tree_Construction.md](Tree_Construction.md)
- Related algorithms: [Distance_Matrix.md](Distance_Matrix.md)

## 8. References

1. Wikipedia contributors. Robinson-Foulds metric. Wikipedia. https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric
2. Wikipedia contributors. Most recent common ancestor. Wikipedia. https://en.wikipedia.org/wiki/Most_recent_common_ancestor
3. Wikipedia contributors. Phylogenetic tree. Wikipedia. https://en.wikipedia.org/wiki/Phylogenetic_tree
4. Robinson, D. F., and L. R. Foulds. 1981. Comparison of phylogenetic trees. Mathematical Biosciences 53(1-2): 131-147. doi:10.1016/0025-5564(81)90043-2
5. Smith, M. R. 2020. Information theoretic Generalized Robinson-Foulds metrics. Bioinformatics 36(20): 5007-5013. doi:10.1093/bioinformatics/btaa614
