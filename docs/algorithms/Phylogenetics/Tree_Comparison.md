# Tree Comparison Algorithms

**Test Unit:** PHYLO-COMP-001
**Category:** Phylogenetics
**Implementation:** `PhylogeneticAnalyzer.cs`

---

## Overview

Tree comparison algorithms measure similarity or distance between phylogenetic trees. This document covers three complementary methods:

1. **Robinson-Foulds Distance** - Topological comparison based on splits
2. **Most Recent Common Ancestor (MRCA)** - Finding shared ancestors
3. **Patristic Distance** - Path-based distance between taxa

---

## 1. Robinson-Foulds Distance

### Definition

The Robinson-Foulds (RF) metric measures the distance between two phylogenetic trees by counting the number of splits (bipartitions) that differ between them.

**Formula:**
$$RF(T_1, T_2) = |S_1 \triangle S_2| = |S_1 \setminus S_2| + |S_2 \setminus S_1|$$

Where:
- $S_1$, $S_2$ are sets of splits induced by trees $T_1$, $T_2$
- $\triangle$ denotes symmetric difference

### Theory (Wikipedia, Robinson & Foulds 1981)

A **split** (or bipartition) is created by removing an internal edge from a tree, which divides the taxa into two disjoint sets. Each internal edge induces exactly one split.

**Properties:**
- RF distance is a true metric (identity, symmetry, triangle inequality)
- For binary trees with n taxa: maximum RF = 2(n-3)
- Identical trees: RF = 0
- Computable in O(n) time (Day 1985)

### Implementation

```csharp
public static int RobinsonFouldsDistance(PhyloNode tree1, PhyloNode tree2)
```

**Algorithm:**
1. Collect all splits from tree1 into set S1
2. Collect all splits from tree2 into set S2
3. Represent splits canonically (smaller partition side, lexicographically sorted)
4. Return |S1 - S2| + |S2 - S1|

### Complexity

- **Time:** O(n) where n is number of taxa
- **Space:** O(n) for storing splits

---

## 2. Most Recent Common Ancestor (MRCA)

### Definition

The MRCA of two taxa is the deepest node in a rooted phylogenetic tree that is an ancestor of both taxa.

### Theory (Wikipedia)

In a rooted phylogenetic tree:
- Each internal node represents the MRCA of all its descendants
- MRCA(x, y) is the node where paths from root to x and y diverge
- MRCA is symmetric: MRCA(x, y) = MRCA(y, x)

### Implementation

```csharp
public static PhyloNode? FindMRCA(PhyloNode root, string taxon1, string taxon2)
```

**Algorithm (recursive):**
1. If node is null, return null
2. If node is leaf and matches either taxon, return node
3. Recursively find MRCA in left and right subtrees
4. If both subtrees return non-null, current node is MRCA
5. Otherwise, return whichever subtree result is non-null

### Complexity

- **Time:** O(n) worst case traversal
- **Space:** O(h) stack space where h is tree height

---

## 3. Patristic Distance

### Definition

The patristic distance between two taxa is the sum of branch lengths along the path connecting them through their MRCA.

**Formula:**
$$PD(x, y) = d(x, MRCA(x,y)) + d(y, MRCA(x,y))$$

### Theory

Patristic distance reflects:
- Evolutionary divergence along tree branches
- Proportional to time if tree is ultrametric (molecular clock)
- Requires branch lengths to be meaningful

### Implementation

```csharp
public static double PatristicDistance(PhyloNode root, string taxon1, string taxon2)
```

**Algorithm:**
1. Find MRCA of taxon1 and taxon2
2. Sum branch lengths from MRCA to taxon1
3. Sum branch lengths from MRCA to taxon2
4. Return total sum

### Complexity

- **Time:** O(n) for MRCA + O(h) for distance calculation
- **Space:** O(h) stack space

---

## Properties and Invariants

### Robinson-Foulds
| Property | Invariant |
|----------|-----------|
| Identity | RF(T, T) = 0 |
| Symmetry | RF(T1, T2) = RF(T2, T1) |
| Non-negativity | RF(T1, T2) ≥ 0 |
| Even | RF is always even |
| Bounded | RF ≤ 2(n-3) for binary trees |

### MRCA
| Property | Invariant |
|----------|-----------|
| Self-MRCA | MRCA(x, x) = x |
| Symmetry | MRCA(x, y) = MRCA(y, x) |
| Uniqueness | Exactly one MRCA per taxon pair |
| Ancestry | MRCA is ancestor of both taxa |

### Patristic Distance
| Property | Invariant |
|----------|-----------|
| Identity | PD(x, x) = 0 |
| Symmetry | PD(x, y) = PD(y, x) |
| Non-negativity | PD(x, y) ≥ 0 |
| Triangle inequality | PD(x, z) ≤ PD(x, y) + PD(y, z) |

---

## Edge Cases

| Case | RobinsonFoulds | FindMRCA | PatristicDistance |
|------|----------------|----------|-------------------|
| Null tree/root | Handle gracefully | Returns null | Returns NaN |
| Same taxon | N/A | Returns taxon | Returns 0 |
| Non-existent taxon | N/A | Returns null | Returns NaN |
| Single taxon tree | RF = 0 | Returns taxon | Returns 0 |
| Zero branch lengths | N/A | N/A | Returns 0 |

---

## References

1. Robinson, D.F.; Foulds, L.R. (1981). "Comparison of phylogenetic trees". Mathematical Biosciences. 53(1-2): 131-147.
2. Day, W.H.E. (1985). "Optimal algorithms for comparing trees with labeled leaves". Journal of Classification. 2(1): 7-28.
3. Wikipedia contributors. "Robinson–Foulds metric". Wikipedia.
4. Wikipedia contributors. "Most recent common ancestor". Wikipedia.
5. Smith, M.R. (2020). "Information theoretic Generalized Robinson-Foulds metrics". Bioinformatics. 36(20): 5007-5013.
