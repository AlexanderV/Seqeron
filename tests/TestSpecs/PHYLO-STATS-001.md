# Test Specification: PHYLO-STATS-001

**Test Unit ID:** PHYLO-STATS-001
**Area:** Phylogenetic
**Algorithm:** Tree Statistics (leaves, total tree length, tree height/depth)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wikipedia — Tree (graph theory) | 4 | https://en.wikipedia.org/wiki/Tree_(graph_theory) | 2026-06-13 |
| 2 | Wikipedia — Tree (abstract data type) | 4 | https://en.wikipedia.org/wiki/Tree_(abstract_data_type) | 2026-06-13 |
| 3 | Biopython `Bio.Phylo.BaseTree` | 3 | https://biopython.org/docs/latest/api/Bio.Phylo.BaseTree.html | 2026-06-13 |
| 4 | DendroPy `Tree` model | 3 | https://dendropy.org/library/treemodel.html | 2026-06-13 |
| 5 | Wikipedia — Minimum evolution (Rzhetsky & Nei 1992) | 4 | https://en.wikipedia.org/wiki/Minimum_evolution | 2026-06-13 |

### 1.2 Key Evidence Points

1. A leaf is "a vertex with no children" — Source 1; Biopython `get_terminals` returns "all of this tree's terminal (leaf) nodes" — Source 3.
2. Tree length = "sum of edge lengths" (DendroPy `Tree.length()` — Source 4) = "sum of all the branch lengths in this tree" (Biopython `total_branch_length` — Source 3); it is the quantity minimized by minimum evolution (Rzhetsky & Nei 1992 — Source 5).
3. Tree height = "the length of the longest downward path to a leaf from [the root]" (edges) — Sources 1, 2.
4. Single-node tree has height/depth 0; empty tree has height −1 — Sources 1, 2.
5. Undefined/missing edge length is treated as 0 when summing tree length — Source 4.

### 1.3 Documented Corner Cases

- Single-node (leaf-only) tree → height 0, one leaf (Sources 1, 2).
- Empty tree → height −1 (Sources 1, 2).
- Edge with no defined length → counted as 0 in the length sum (Source 4).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing height in **edges** with height in **branch-length distance** — Biopython `depths()` is "by branch length", a different metric; this unit's `GetTreeDepth` is the topological (edge-count) height — Sources 2, 3.
2. Returning 0 for both a leaf and a null/empty tree conflates two distinct cases; the empty tree is −1 — Sources 1, 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetLeaves(PhyloNode root)` | PhylogeneticAnalyzer | Canonical | Returns terminal nodes (no children). |
| `CalculateTreeLength(PhyloNode root)` | PhylogeneticAnalyzer | Canonical | Sum of all branch lengths. |
| `GetTreeDepth(PhyloNode root)` | PhylogeneticAnalyzer | Canonical | Topological height in edges to deepest leaf. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every node returned by `GetLeaves` satisfies `IsLeaf` (no children). | Yes | Source 1, 3 |
| INV-2 | For an N-leaf fully-bifurcating tree, `GetLeaves(root).Count() == N`. | Yes | Source 1, 3 |
| INV-3 | `CalculateTreeLength(root)` equals Σ over all nodes of `node.BranchLength`. | Yes | Source 3, 4 |
| INV-4 | `CalculateTreeLength(root) >= 0` when all branch lengths are ≥ 0. | Yes | Source 4 |
| INV-5 | `GetTreeDepth` = number of edges on the longest root-to-leaf path; leaf-only tree = 0. | Yes | Source 1, 2 |
| INV-6 | `GetTreeDepth(null) == -1` (empty-tree convention); `GetLeaves(null)` empty; `CalculateTreeLength(null) == 0`. | Yes | Source 1, 2; ASSUMPTION-1 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Leaves of balanced 4-taxon tree | `GetLeaves` on `((A,B),(C,D))` | Exactly {A,B,C,D}, all `IsLeaf`, in order A,B,C,D | Source 1, 3 |
| M2 | Leaf count = N | `GetLeaves(...).Count()` on 4-leaf tree | 4 | Source 3 (count_terminals) |
| M3 | Single-leaf tree leaves | `GetLeaves` on a lone leaf node | The single node itself | Source 1 |
| M4 | Length of balanced tree | `CalculateTreeLength` of `((A:1,B:1):1,(C:1,D:1):1)` | 6.0 | Source 3, 4 |
| M5 | Length of caterpillar tree | `CalculateTreeLength` of `(A:1,(B:1,(C:1,D:1):0.5):0.5)` | 5.0 | Source 3, 4 |
| M6 | Length sums root branch too | Length includes the root node's own `BranchLength` | Σ all branch lengths | Source 3, 4 |
| M7 | Height of balanced tree (edges) | `GetTreeDepth` of `((A,B),(C,D))` | 2 | Source 1, 2 |
| M8 | Height of caterpillar tree | `GetTreeDepth` of `(A,(B,(C,D)))` | 3 | Source 1, 2 |
| M9 | Single-leaf tree height | `GetTreeDepth` of a lone leaf | 0 | Source 1, 2 |
| M10 | Null tree leaves/length/height | `GetLeaves(null)`, `CalculateTreeLength(null)`, `GetTreeDepth(null)` | ∅, 0.0, −1 | Source 1, 2 (ASSUMPTION-1) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Default (zero) branch lengths | Length of a tree whose branch lengths are unset (0) | 0.0 | DendroPy: missing length → 0 |
| S2 | Two-leaf tree height | `GetTreeDepth` of `(A,B)` | 1 | One edge root→leaf |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Leaves order is left-to-right | Pre-order leaf traversal order | A,B,C,D for `((A,B),(C,D))` | Traversal contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for existing coverage of `GetLeaves`, `CalculateTreeLength`, `GetTreeDepth`. No dedicated test file existed for PHYLO-STATS-001 prior to this unit; sibling file `PhylogeneticAnalyzer_Bootstrap_Tests.cs` covers a different method.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New unit |
| M2 | ❌ Missing | New unit |
| M3 | ❌ Missing | New unit |
| M4 | ❌ Missing | New unit |
| M5 | ❌ Missing | New unit |
| M6 | ❌ Missing | New unit |
| M7 | ❌ Missing | New unit |
| M8 | ❌ Missing | New unit |
| M9 | ❌ Missing | New unit |
| M10 | ❌ Missing | New unit |
| S1 | ❌ Missing | New unit |
| S2 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_TreeStatistics_Tests.cs` — all PHYLO-STATS-001 cases.
- **Remove:** nothing (no prior tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PhylogeneticAnalyzer_TreeStatistics_Tests.cs` | Canonical for PHYLO-STATS-001 | 13 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | S1 | ❌ Missing | Implemented | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `GetLeaves_BalancedTree_ReturnsAllTerminalNodesInOrder` |
| M2 | ✅ Covered | `GetLeaves_FourLeafTree_CountEqualsFour` |
| M3 | ✅ Covered | `GetLeaves_SingleLeaf_ReturnsThatNode` |
| M4 | ✅ Covered | `CalculateTreeLength_BalancedTree_SumsAllBranchLengths` |
| M5 | ✅ Covered | `CalculateTreeLength_CaterpillarTree_SumsAllBranchLengths` |
| M6 | ✅ Covered | `CalculateTreeLength_IncludesRootBranchLength` |
| M7 | ✅ Covered | `GetTreeDepth_BalancedTree_ReturnsEdgeCountToDeepestLeaf` |
| M8 | ✅ Covered | `GetTreeDepth_CaterpillarTree_ReturnsThree` |
| M9 | ✅ Covered | `GetTreeDepth_SingleLeaf_ReturnsZero` |
| M10 | ✅ Covered | `TreeStatistics_NullRoot_ReturnEmptyZeroAndMinusOne` |
| S1 | ✅ Covered | `CalculateTreeLength_DefaultBranchLengths_ReturnsZero` |
| S2 | ✅ Covered | `GetTreeDepth_TwoLeafTree_ReturnsOne` |
| C1 | ✅ Covered | `GetLeaves_BalancedTree_PreOrderTraversalOrder` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Null `PhyloNode` maps to the graph-theory "empty tree" → `GetTreeDepth(null)` = −1, `GetLeaves(null)` = ∅, `CalculateTreeLength(null)` = 0. | INV-6, M10 |

---

## 7. Open Questions / Decisions

1. Decision: `GetTreeDepth` reports topological height in **edges** (not branch-length distance), per the graph-theory/ADT definition (Sources 1, 2). Biopython's branch-length `depths()` is a distinct metric and out of scope for this unit.
