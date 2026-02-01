# PHYLO-TREE-001: Tree Construction - Evidence Document

**Test Unit ID:** PHYLO-TREE-001  
**Algorithm:** Phylogenetic Tree Construction (UPGMA, Neighbor-Joining)  
**Date:** 2026-02-01  
**Status:** Complete

---

## 1. Authoritative Sources

### Primary Sources

| Source | Type | URL |
|--------|------|-----|
| Wikipedia: UPGMA | Encyclopedia | https://en.wikipedia.org/wiki/UPGMA |
| Wikipedia: Neighbor joining | Encyclopedia | https://en.wikipedia.org/wiki/Neighbor_joining |
| Wikipedia: Phylogenetic tree | Encyclopedia | https://en.wikipedia.org/wiki/Phylogenetic_tree |
| Saitou & Nei (1987) | Original Paper | Molecular Biology and Evolution 4(4):406-425 |
| Sokal & Michener (1958) | Original Paper | University of Kansas Science Bulletin 38:1409-1438 |

### Secondary Sources

| Source | Type | Relevance |
|--------|------|-----------|
| Felsenstein (2004) | Textbook | Inferring Phylogenies, Sinauer Associates |

---

## 2. Algorithm Definition

### 2.1 UPGMA (Unweighted Pair Group Method with Arithmetic Mean)

From Wikipedia (UPGMA):
> "UPGMA is a simple agglomerative (bottom-up) hierarchical clustering method... At each step, the nearest two clusters are combined into a higher-level cluster."

**Key properties:**
- Produces rooted ultrametric trees
- Assumes molecular clock (constant rate of evolution)
- Time complexity: O(n³) naive, O(n²) optimized

**Algorithm steps:**
1. Initialize each taxon as its own cluster
2. Find pair of clusters with minimum distance
3. Join into new cluster; compute new distances as weighted average
4. Repeat until single cluster remains

**Branch length formula:**
- Height = distance / 2 (ultrametric property)

### 2.2 Neighbor-Joining Algorithm

From Wikipedia (Neighbor joining):
> "Neighbor joining is a bottom-up (agglomerative) clustering method for the creation of phylogenetic trees, created by Naruya Saitou and Masatoshi Nei in 1987."

**Key properties:**
- Produces unrooted trees (rooted by convention)
- Does NOT assume molecular clock
- Time complexity: O(n³)
- Guarantees correct topology for additive distance matrices

**Algorithm steps:**
1. Calculate Q-matrix: Q(i,j) = (n-2)·d(i,j) - Σd(i,k) - Σd(j,k)
2. Find pair (i,j) with minimum Q value
3. Calculate branch lengths using formula from source
4. Update distance matrix
5. Repeat until tree is resolved

**Branch length formulas (from Wikipedia):**
```
δ(f,u) = d(f,g)/2 + (Σd(f,k) - Σd(g,k)) / (2(n-2))
δ(g,u) = d(f,g) - δ(f,u)
d(u,k) = (d(f,k) + d(g,k) - d(f,g)) / 2
```

---

## 3. Key Invariants

### 3.1 Tree Structure Invariants

| ID | Invariant | Source |
|----|-----------|--------|
| INV-01 | Tree contains all input taxa as leaves | Definition |
| INV-02 | Binary tree: each internal node has exactly 2 children | UPGMA/NJ definition |
| INV-03 | Number of leaves = n (input sequences) | Definition |
| INV-04 | Number of internal nodes = n-1 (for binary tree) | Graph theory |
| INV-05 | All branch lengths ≥ 0 | UPGMA definition |

### 3.2 UPGMA-Specific Invariants

| ID | Invariant | Source |
|----|-----------|--------|
| INV-U01 | Ultrametric: all tips equidistant from root | Wikipedia UPGMA |
| INV-U02 | Produces rooted tree | Wikipedia UPGMA |
| INV-U03 | Height = distance/2 | UPGMA algorithm |

### 3.3 Neighbor-Joining Invariants

| ID | Invariant | Source |
|----|-----------|--------|
| INV-N01 | Correct topology for additive matrices | Wikipedia NJ |
| INV-N02 | May produce negative branch lengths | Wikipedia NJ (limitation) |
| INV-N03 | Does not require clock assumption | Wikipedia NJ |

---

## 4. Test Data from Sources

### 4.1 UPGMA Example (Wikipedia)

From Wikipedia UPGMA working example (5S rRNA):
```
Initial distance matrix D1:
    a    b    c    d    e
a   0   17   21   31   23
b  17    0   30   34   21
c  21   30    0   28   39
d  31   34   28    0   43
e  23   21   39   43    0
```

**Expected clustering order:**
1. (a,b) merged first (d=17)
2. ((a,b),e) merged (d=22)
3. (c,d) merged (d=28)
4. Final join at d=33

**Expected branch lengths:**
- δ(a,u) = δ(b,u) = 8.5
- δ(e,v) = 11
- δ(c,w) = δ(d,w) = 14
- Root height = 16.5

### 4.2 Neighbor-Joining Example (Wikipedia)

From Wikipedia NJ working example (5 taxa):
```
Initial distance matrix:
    a   b   c   d   e
a   0   5   9   9   8
b   5   0  10  10   9
c   9  10   0   8   7
d   9  10   8   0   3
e   8   9   7   3   0
```

**Expected branch lengths:**
- δ(a,u) = 2
- δ(b,u) = 3
- δ(u,v) = 3
- δ(c,v) = 4
- δ(v,w) = 2
- δ(d,w) = 2
- δ(e,w) = 1

---

## 5. Edge Cases

### 5.1 Input Validation

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| < 2 sequences | Exception | Definition (need pairs) |
| Unequal lengths | Exception | Alignment required |
| 2 sequences | Trivial binary tree | Minimum case |
| Identical sequences | Zero distance, arbitrary join order | Algorithm |

### 5.2 Distance Matrix Properties

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| All zeros | All taxa identical | Algorithm |
| All equal distances | Star topology (arbitrary resolution) | Algorithm |
| Saturated distances (p > 0.75) | JC returns infinity | JC formula limit |

### 5.3 Degenerate Cases

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Single nucleotide sequences | Valid but trivial | Algorithm |
| All gaps | Zero comparable sites | Implementation |

---

## 6. Testing Methodology

Based on the evidence, tests should verify:

1. **Correctness tests** - Use known examples from Wikipedia
2. **Invariant tests** - Verify tree structure properties
3. **Edge case tests** - Validate error handling
4. **Property tests** - All taxa present, binary structure
5. **Method selection tests** - UPGMA vs NJ produce valid trees

---

## 7. Implementation Notes

The implementation in `PhylogeneticAnalyzer.cs`:
- `BuildTree()` is the canonical method
- Supports both UPGMA and NeighborJoining methods
- Returns `PhylogeneticTree` with Root, Taxa, DistanceMatrix, Method
- Uses distance matrix internally
- Branch lengths are computed according to standard formulas

---

## References

1. Sokal RR, Michener CD (1958). "A statistical method for evaluating systematic relationships." University of Kansas Science Bulletin 38:1409-1438.
2. Saitou N, Nei M (1987). "The neighbor-joining method: a new method for reconstructing phylogenetic trees." Molecular Biology and Evolution 4(4):406-425.
3. Wikipedia contributors. "UPGMA." Wikipedia, The Free Encyclopedia.
4. Wikipedia contributors. "Neighbor joining." Wikipedia, The Free Encyclopedia.
5. Wikipedia contributors. "Phylogenetic tree." Wikipedia, The Free Encyclopedia.
