# Tree Construction Algorithms

**Category:** Phylogenetics  
**Test Unit:** PHYLO-TREE-001  
**Last Updated:** 2026-02-01

---

## Overview

Phylogenetic tree construction algorithms build evolutionary trees from distance matrices calculated from aligned sequences. The implementation supports two distance-based methods: UPGMA and Neighbor-Joining.

## Algorithms

### UPGMA (Unweighted Pair Group Method with Arithmetic Mean)

UPGMA is a bottom-up hierarchical clustering algorithm that produces rooted ultrametric trees.

**Properties:**
- Time complexity: O(n³) naive, O(n²) optimized
- Space complexity: O(n²)
- Produces rooted trees
- Assumes molecular clock (constant evolutionary rate)
- All tips equidistant from root (ultrametric)

**Algorithm:**
1. Initialize each taxon as a singleton cluster
2. Find the pair of clusters (i, j) with minimum distance
3. Merge into a new cluster u
4. Set branch lengths: δ(i,u) = δ(j,u) = d(i,j)/2
5. Update distances using weighted average:
   ```
   d(u,k) = (|i|·d(i,k) + |j|·d(j,k)) / (|i| + |j|)
   ```
6. Repeat until single cluster remains

**Source:** Sokal & Michener (1958)

### Neighbor-Joining

Neighbor-Joining is a bottom-up algorithm that does not assume a molecular clock, making it suitable for data with varying evolutionary rates.

**Properties:**
- Time complexity: O(n³)
- Space complexity: O(n²)
- Produces unrooted trees (rooted by convention)
- No clock assumption required
- Guarantees correct topology for additive distance matrices
- May produce negative branch lengths

**Algorithm:**
1. Calculate Q-matrix:
   ```
   Q(i,j) = (n-2)·d(i,j) - Σd(i,k) - Σd(j,k)
   ```
2. Find pair (f, g) with minimum Q value
3. Create new node u joining f and g
4. Calculate branch lengths:
   ```
   δ(f,u) = d(f,g)/2 + (Σd(f,k) - Σd(g,k)) / (2(n-2))
   δ(g,u) = d(f,g) - δ(f,u)
   ```
5. Update distances:
   ```
   d(u,k) = (d(f,k) + d(g,k) - d(f,g)) / 2
   ```
6. Repeat until fully resolved

**Source:** Saitou & Nei (1987)

## Implementation

### Entry Point

```csharp
public static PhylogeneticTree BuildTree(
    IReadOnlyDictionary<string, string> sequences,
    DistanceMethod distanceMethod = DistanceMethod.JukesCantor,
    TreeMethod treeMethod = TreeMethod.UPGMA)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| sequences | IReadOnlyDictionary<string, string> | required | Named aligned sequences |
| distanceMethod | DistanceMethod | JukesCantor | Distance calculation method |
| treeMethod | TreeMethod | UPGMA | Tree building algorithm |

### Return Value

```csharp
public readonly record struct PhylogeneticTree(
    PhyloNode Root,
    IReadOnlyList<string> Taxa,
    double[,] DistanceMatrix,
    string Method);
```

### Tree Node Structure

```csharp
public class PhyloNode
{
    public string Name { get; set; }
    public double BranchLength { get; set; }
    public PhyloNode? Left { get; set; }
    public PhyloNode? Right { get; set; }
    public bool IsLeaf => Left == null && Right == null;
    public List<string> Taxa { get; set; }
}
```

## Invariants

| ID | Invariant | Description |
|----|-----------|-------------|
| INV-01 | All taxa present | Tree contains all input taxa as leaves |
| INV-02 | Binary structure | Each internal node has exactly 2 children |
| INV-03 | Non-negative lengths | Branch lengths ≥ 0 (UPGMA guarantee) |
| INV-04 | Ultrametric (UPGMA) | All tips equidistant from root |
| INV-05 | Correct topology | NJ produces correct tree for additive matrices |

## Edge Cases

| Case | Behavior |
|------|----------|
| < 2 sequences | ArgumentException |
| Unequal lengths | ArgumentException |
| 2 sequences | Simple binary tree |
| Identical sequences | Valid tree, zero branch lengths |
| Saturated distance | JC returns infinity (high divergence) |

## Comparison

| Feature | UPGMA | Neighbor-Joining |
|---------|-------|------------------|
| Tree type | Rooted | Unrooted |
| Clock assumption | Required | Not required |
| Branch lengths | Always positive | May be negative |
| Best for | Similar rates | Variable rates |
| Accuracy | Lower | Higher |

## References

1. Saitou N, Nei M (1987). "The neighbor-joining method: a new method for reconstructing phylogenetic trees." Molecular Biology and Evolution 4(4):406-425.
2. Sokal RR, Michener CD (1958). "A statistical method for evaluating systematic relationships." University of Kansas Science Bulletin 38:1409-1438.
3. Wikipedia: [UPGMA](https://en.wikipedia.org/wiki/UPGMA)
4. Wikipedia: [Neighbor joining](https://en.wikipedia.org/wiki/Neighbor_joining)
