# Tree Construction Algorithms

> **Baseline / reference method.** UPGMA (documented here alongside Neighbor-Joining) is retained as a standard baseline — it assumes a molecular clock and can bias topology/branch lengths, so Neighbor-Joining is the default choice. See [Legacy / Baseline Methods](../CANONICAL_MAP.md).

| Field | Value |
|-------|-------|
| Algorithm Group | Phylogenetics |
| Test Unit ID | PHYLO-TREE-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Distance-based phylogenetic tree construction converts a distance matrix into a binary tree over the input taxa. The repository supports two classical methods: UPGMA and Neighbor-Joining. Both are deterministic agglomerative procedures, but they make different evolutionary assumptions and therefore suit different data regimes. The current implementation accepts either aligned sequences or a precomputed distance matrix and returns a `PhylogeneticTree` whose `Root` is a binary `PhyloNode` structure.

## 2. Scientific / Formal Basis

> A = UPGMA, B = Neighbor-Joining

### 2.A UPGMA

#### Domain Context

UPGMA is a hierarchical clustering method that builds a rooted ultrametric tree from a distance matrix. It is appropriate only when the data satisfy a molecular-clock assumption so that all leaves should end up equidistant from the root.

#### Core Model

UPGMA repeatedly merges the closest two clusters and updates distances by weighted averaging:

$$
d(u, k) = \frac{|i| d(i, k) + |j| d(j, k)}{|i| + |j|}
$$

with new-cluster height:

$$
h(u) = \frac{d(i, j)}{2}
$$

and child branch lengths computed from the height difference to the new cluster.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-UPGMA-01 | Evolution follows a molecular clock | Non-ultrametric data can lead to incorrect topology and branch lengths |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-UPGMA-01 | Output is rooted and binary | The algorithm repeatedly merges two active clusters into one internal node |
| INV-UPGMA-02 | All tips are equidistant from the root in the ideal model | UPGMA assigns branch lengths from cluster heights |
| INV-UPGMA-03 | Branch lengths are non-negative in the implementation | Child branch lengths are computed as `Math.Max(0, newHeight - childHeight)` |

### 2.B Neighbor-Joining

#### Domain Context

Neighbor-Joining (NJ) builds a tree without assuming a molecular clock. It is designed to recover good topology from additive distance matrices even when evolutionary rates vary among lineages (Saitou and Nei, 1987).

#### Core Model

NJ chooses the next pair to join by minimizing the $Q$-matrix:

$$
Q(i, j) = (n - 2) d(i, j) - \sum_k d(i, k) - \sum_k d(j, k)
$$

with branch lengths:

$$
\delta(i, u) = \frac{d(i, j)}{2} + \frac{\sum_k d(i, k) - \sum_k d(j, k)}{2(n - 2)}
$$

$$
\delta(j, u) = d(i, j) - \delta(i, u)
$$

and updated distances:

$$
d(u, k) = \frac{d(i, k) + d(j, k) - d(i, j)}{2}
$$

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-NJ-01 | The distance matrix is approximately additive | Tree topology can become less reliable when the matrix deviates strongly from additivity |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-NJ-01 | Output is binary | Each join creates a new internal node with two children |
| INV-NJ-02 | Negative branch lengths are possible | The NJ branch-length formula is used directly without clamping |
| INV-NJ-03 | No molecular clock is assumed | The algorithm is formulated for non-ultrametric distance matrices |

#### Comparison with Related Methods

| Aspect | UPGMA | Neighbor-Joining |
|--------|-------|------------------|
| Tree type in theory | Rooted, ultrametric | Unrooted |
| Molecular-clock assumption | Required | Not required |
| Branch lengths | Non-negative in this implementation | May be negative |
| Best fit | Similar rates across lineages | Variable rates across lineages |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[BuildTree] sequences` | `IReadOnlyDictionary<string, string>` | required | Named aligned sequences | At least two sequences are required and all lengths must match |
| `[BuildTree] distanceMethod` | `DistanceMethod` | `JukesCantor` | Distance method used to compute the matrix from aligned sequences | Passed to `CalculateDistanceMatrix` |
| `[BuildTree] treeMethod` | `TreeMethod` | `UPGMA` | Tree-building method | One of `UPGMA` or `NeighborJoining` |
| `[BuildTreeFromMatrix] taxa` | `IReadOnlyList<string>` | required | Taxon names in matrix order | Must contain at least two taxa |
| `[BuildTreeFromMatrix] distanceMatrix` | `double[,]` | required | Precomputed symmetric distance matrix | Dimensions must match the taxon count |
| `[BuildTreeFromMatrix] treeMethod` | `TreeMethod` | `UPGMA` | Tree-building method | One of `UPGMA` or `NeighborJoining` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Root` | `PhyloNode` | Root node of the constructed binary tree |
| `Taxa` | `IReadOnlyList<string>` | Taxa in input order |
| `DistanceMatrix` | `double[,]` | Matrix used by the construction routine |
| `Method` | `string` | Name of the tree-building method |

### 3.3 Preconditions and Validation

`BuildTree` throws `ArgumentException` when fewer than two sequences are supplied or when the sequences are not aligned to equal length. `BuildTreeFromMatrix` throws `ArgumentException` when fewer than two taxa are supplied, when the matrix is missing, or when matrix dimensions do not match the taxon list. `BuildTree` computes the distance matrix first and then dispatches to the chosen tree builder.

## 4. Algorithm

### 4.A UPGMA

#### High-Level Steps

1. Initialize each taxon as a singleton cluster with height `0`.
2. Find the active cluster pair with minimum distance.
3. Create a new parent node joining that pair.
4. Set the new cluster height to `d(i, j) / 2`.
5. Assign child branch lengths from the new height minus child heights.
6. Update distances by weighted averaging and repeat until one cluster remains.

#### Decision Rules / Reference Tables

The implementation stores working distances in a dictionary keyed by cluster indices and keeps cluster size and cluster height maps to compute weighted averages and ultrametric branch lengths.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| UPGMA tree construction | `O(n^3)` | `O(n^2)` | `n` = number of taxa |

### 4.B Neighbor-Joining

#### High-Level Steps

1. Initialize one leaf node per taxon.
2. Compute the `r` sums and the `Q` matrix over active taxa.
3. Join the pair with minimum `Q` value.
4. Compute the two child branch lengths.
5. Update distances from the new node to the remaining active taxa.
6. Repeat until two active taxa remain, then join them under a final root node.

#### Decision Rules / Reference Tables

The repository returns a rooted binary `PhyloNode` even for NJ by joining the last two active nodes under a final root. This is a representation convention for the in-memory tree structure rather than a theoretical claim that NJ is rooted.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Neighbor-Joining tree construction | `O(n^3)` | `O(n^2)` | `n` = number of taxa |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PhylogeneticAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs)

- `PhylogeneticAnalyzer.BuildTree(IReadOnlyDictionary<string, string>, DistanceMethod, TreeMethod)`: Computes a distance matrix from aligned sequences and builds a tree.
- `PhylogeneticAnalyzer.BuildTreeFromMatrix(IReadOnlyList<string>, double[,], TreeMethod)`: Builds a tree from a precomputed matrix.

### 5.2 Current Behavior

The repository first computes a distance matrix unless one is supplied directly. UPGMA uses incremental cluster heights and explicitly clamps child branch lengths to non-negative values. Neighbor-Joining uses the standard branch-length formula directly, so negative branch lengths can appear. For NJ, the in-memory representation becomes rooted only at the final join step when the last two active nodes are attached to a new parent node.

### 5.3 Conformance to Theory / Spec

#### 5.3.A UPGMA

**Implemented (verbatim from the cited theory/spec):**

- Repeated joining of the minimum-distance cluster pair.
- Weighted-average distance updates.
- Branch-length assignment from cluster heights.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

#### 5.3.B Neighbor-Joining

**Implemented (verbatim from the cited theory/spec):**

- Q-matrix selection rule.
- Standard NJ branch-length formulas.
- Standard distance update rule for newly joined nodes.

**Intentionally simplified:**

- The final in-memory result is rooted by convention so it fits the binary `PhyloNode` API; **consequence:** the returned object is a rooted representation of an algorithm that is theoretically unrooted.

**Not implemented:**

- (none)

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Rooted `PhyloNode` representation for NJ output | Deviation | Consumers see a rooted object even though NJ itself is an unrooted method | accepted | Representation choice made at the final join step |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Fewer than two sequences or taxa | `ArgumentException` | A binary tree requires at least two leaves |
| Unequal aligned-sequence lengths | `ArgumentException` | Distance calculation requires aligned inputs |
| Exactly two taxa | Simple binary tree | One join is sufficient |
| Identical sequences | Valid tree with zero distances and zero or equal branch lengths | Distances collapse to zero |
| Saturated nucleotide distance from the chosen model | Matrix may contain positive infinity | The upstream distance method can return infinity |

### 6.2 Limitations

These builders are classical `O(n^3)` implementations with a binary `PhyloNode` output model. They do not include bootstrap support values, multifurcations, or advanced optimization strategies. Interpretation quality still depends on the suitability of the chosen distance model and on the molecular-clock/additivity assumptions of the selected method.

## 7. Examples and Related Material

- [PHYLO-TREE-001](../../../tests/TestSpecs/PHYLO-TREE-001.md) documents the repository's tree-construction test specification.
- [Distance_Matrix.md](./Distance_Matrix.md) documents the distance models used to create the input matrix.

## 8. References

1. Saitou, N., and M. Nei. 1987. The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4(4):406-425.
2. Sokal, R. R., and C. D. Michener. 1958. A statistical method for evaluating systematic relationships. University of Kansas Science Bulletin 38:1409-1438.
3. Wikipedia contributors. UPGMA. Wikipedia. https://en.wikipedia.org/wiki/UPGMA
4. Wikipedia contributors. Neighbor joining. Wikipedia. https://en.wikipedia.org/wiki/Neighbor_joining
