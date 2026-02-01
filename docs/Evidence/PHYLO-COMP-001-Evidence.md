# Evidence Document: PHYLO-COMP-001 (Tree Comparison)

**Test Unit ID:** PHYLO-COMP-001
**Algorithm:** Tree Comparison (Robinson-Foulds Distance, MRCA, Patristic Distance)
**Date:** 2026-02-01
**Status:** Evidence Gathered

---

## 1. Sources Consulted

### Primary Sources

| # | Source | Type | URL | Accessed |
|---|--------|------|-----|----------|
| 1 | Wikipedia: Robinson–Foulds metric | Encyclopedia | https://en.wikipedia.org/wiki/Robinson–Foulds_metric | 2026-02-01 |
| 2 | Wikipedia: Most recent common ancestor | Encyclopedia | https://en.wikipedia.org/wiki/Most_recent_common_ancestor | 2026-02-01 |
| 3 | Wikipedia: Phylogenetic tree | Encyclopedia | https://en.wikipedia.org/wiki/Phylogenetic_tree | 2026-02-01 |
| 4 | Robinson & Foulds (1981) | Original Paper | doi:10.1016/0025-5564(81)90043-2 | 2026-02-01 |
| 5 | Smith (2020) - Information theoretic Generalized RF metrics | Paper | doi:10.1093/bioinformatics/btaa614 | 2026-02-01 |

---

## 2. Algorithm Definitions

### 2.1 Robinson-Foulds Distance

**Definition (Wikipedia):**
> "The Robinson–Foulds or symmetric difference metric [...] is a simple way to calculate the distance between phylogenetic trees. It is defined as (A + B) where A is the number of partitions of data implied by the first tree but not the second tree and B is the number of partitions of data implied by the second tree but not the first tree."

**Key Properties (Robinson & Foulds 1981):**
- RF distance is a proper metric (satisfies identity, symmetry, triangle inequality)
- Computed via splits (bipartitions) induced by removing edges
- For n taxa, maximum RF distance = 2(n-3) for unrooted binary trees
- Identical trees → RF distance = 0
- Linear time algorithm exists (Day 1985)

**Implementation Approach:**
1. Extract splits (bipartitions) from each tree
2. Represent splits canonically (using smaller partition side)
3. Count symmetric difference between split sets

### 2.2 Most Recent Common Ancestor (MRCA)

**Definition (Wikipedia):**
> "A most recent common ancestor (MRCA), also known as a last common ancestor (LCA) or concestor, is the most recent individual from which all organisms of a set are inferred to have descended."

**Properties:**
- In a rooted tree, each node represents MRCA of its descendants
- For two taxa: MRCA is the deepest node that is an ancestor of both
- Finding MRCA on a tree: O(n) via recursive traversal
- Same taxon queried twice: MRCA is the taxon itself

### 2.3 Patristic Distance

**Definition (derived from phylogenetic tree literature):**
Patristic distance is the sum of branch lengths along the path connecting two taxa through their MRCA.

**Properties:**
- Patristic distance = distance(taxon1 → MRCA) + distance(MRCA → taxon2)
- Same taxon: patristic distance = 0
- Requires branch lengths to be meaningful
- Reflects evolutionary divergence along tree path

---

## 3. Test Cases from Sources

### 3.1 Robinson-Foulds Distance

| # | Test Case | Source | Expected Behavior |
|---|-----------|--------|-------------------|
| RF-1 | Identical trees | Wikipedia, RF metric definition | RF distance = 0 |
| RF-2 | Completely different topologies | Wikipedia | Maximum RF distance |
| RF-3 | Trees with shared splits | Wikipedia | Count symmetric difference |
| RF-4 | Single taxon difference in split | Wikipedia | RF increases by number of differing splits |

### 3.2 MRCA

| # | Test Case | Source | Expected Behavior |
|---|-----------|--------|-------------------|
| MRCA-1 | Two sibling taxa | Wikipedia, MRCA definition | Returns their parent node |
| MRCA-2 | Same taxon queried twice | Wikipedia | Returns the taxon node itself |
| MRCA-3 | Distant taxa | Wikipedia | Returns deepest common ancestor |
| MRCA-4 | Root is MRCA | Wikipedia | All taxa share root as MRCA |
| MRCA-5 | Non-existent taxon | Edge case | Returns null or handles gracefully |

### 3.3 Patristic Distance

| # | Test Case | Source | Expected Behavior |
|---|-----------|--------|-------------------|
| PD-1 | Same taxon | Definition | Distance = 0 |
| PD-2 | Sibling taxa | Definition | Sum of branch lengths to parent |
| PD-3 | Distant taxa | Definition | Sum through path via MRCA |
| PD-4 | Non-existent taxon | Edge case | Returns NaN or handles gracefully |

---

## 4. Edge Cases and Corner Cases

### 4.1 Robinson-Foulds Edge Cases
- **Null tree:** Should handle gracefully (return 0 or error)
- **Single taxon trees:** No internal edges, RF = 0
- **Two taxa trees:** Only one possible topology, RF = 0 if same
- **Star topology vs binary tree:** Different split counts

### 4.2 MRCA Edge Cases
- **Null root:** Returns null
- **Empty Taxa list:** Implementation-dependent
- **Taxon not in tree:** Returns null

### 4.3 Patristic Distance Edge Cases
- **Zero branch lengths:** Distance = 0 even for different taxa
- **Taxon not in tree:** Returns NaN
- **Single taxon tree:** Only valid query is same taxon → 0

---

## 5. Invariants

### 5.1 Robinson-Foulds Invariants
1. RF(T, T) = 0 (identity)
2. RF(T1, T2) = RF(T2, T1) (symmetry)
3. RF(T1, T2) ≥ 0 (non-negativity)
4. RF is even (symmetric difference of two sets)

### 5.2 MRCA Invariants
1. MRCA(x, x) returns node for x
2. MRCA(x, y) = MRCA(y, x) (symmetry)
3. MRCA is always an ancestor of both taxa
4. MRCA is unique for a given tree and taxon pair

### 5.3 Patristic Distance Invariants
1. PD(x, x) = 0
2. PD(x, y) = PD(y, x) (symmetry)
3. PD(x, y) ≥ 0 (non-negativity)
4. PD(x, y) = dist(x, MRCA) + dist(y, MRCA)

---

## 6. Implementation Notes

### Current Implementation Analysis

**RobinsonFouldsDistance:**
- Uses split-based comparison
- Canonical representation: smaller partition side
- Returns count of symmetric difference

**FindMRCA:**
- Recursive tree traversal
- Returns node when both taxa found in subtrees
- Handles leaf nodes by name matching

**PatristicDistance:**
- Uses FindMRCA to locate common ancestor
- Sums branch lengths from MRCA to each taxon
- Returns NaN for non-existent taxa

---

## 7. Open Questions

1. **Multifurcating trees:** Implementation assumes binary trees; behavior with polytomies?
2. **Unrooted trees:** RF implementation assumes rooted; unrooted support?
3. **Normalized RF:** Should normalized version (0-1 range) be provided?

---

## 8. Summary

Evidence gathered from Wikipedia and original literature provides clear definitions and expected behaviors for all three tree comparison methods. Test cases are derived from documented properties and mathematical invariants of these metrics.
