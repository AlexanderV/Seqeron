# Validation Report: PHYLO-STATS-001 — Tree Statistics (leaves, total tree length, tree height/depth)

- **Validated:** 2026-06-15   **Re-validated (Phase-3, N-ary refactor):** 2026-06-17   **Area:** Phylogenetics
- **Canonical method(s):** `PhylogeneticAnalyzer.GetLeaves(PhyloNode)`, `PhylogeneticAnalyzer.CalculateTreeLength(PhyloNode)`, `PhylogeneticAnalyzer.GetTreeDepth(PhyloNode)`, `PhylogeneticAnalyzer.PatristicDistance(PhyloNode, string, string)` (via `FindMRCA`/`DistanceToTaxon`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** CLEAN (no defect found; N-ary traversal mutation-verified 2026-06-17)

## Stage A — Description

### Sources opened this session (retrieved, not trusted by label)

1. **Wikipedia — Tree (graph theory)** (WebFetch, 2026-06-15)
   - Leaf: "A *leaf* is a vertex with no children." / "An *external vertex* (or outer vertex, terminal vertex or leaf) is a vertex of degree 1."
   - Depth: "The *depth* of a vertex is the length of the path to its root (*root path*)." "The root has depth zero…"
   - Height: "The *height* of a vertex in a rooted tree is the length of the longest downward path to a leaf from that vertex. The *height* of the tree is the height of the root."
   - Single vertex: "…a tree with only a single vertex (hence both a root and leaf) has depth and height zero."
   - Empty tree: "Conventionally, an empty tree (a tree with no vertices, if such are allowed) has depth and height −1."

2. **Wikipedia — Tree (abstract data type)** (WebFetch, 2026-06-15)
   - Leaf: "An **external node** (also known as an **outer node**, **leaf node**, or **terminal node**) is any node that does not have child nodes."
   - Height: "The **height** of a node is the length of the longest downward path to a leaf from that node. The height of the root is the height of the tree."
   - "…the root node has depth zero, leaf nodes have height zero, and a tree with only a single node (hence both a root and leaf) has depth and height zero."
   - Empty tree: "Conventionally, an empty tree (tree with no nodes, if such are allowed) has height −1."

3. **Biopython `Bio.Phylo.BaseTree`** (WebFetch, 2026-06-15)
   - `total_branch_length`: "Calculate the sum of all the branch lengths in this tree."
   - `get_terminals`: "Get a list of all of this tree's terminal (leaf) nodes."
   - `count_terminals`: "Count the number of terminal (leaf) nodes within this tree."
   - `depths()`: distance from root to each clade, by cumulative branch length (default) **or** by branch count when `unit_branch_lengths=True` — confirms branch-length depth is a *distinct* metric from this unit's topological height.

4. **DendroPy `Tree.length()`** (WebFetch, 2026-06-15)
   - "Returns sum of edge lengths of self. Edges with no lengths defined (None) will be considered to have a length of 0." — confirms tree length = Σ edge lengths and the missing-length → 0 convention.

### Formula check

- **Leaf = node with no children** — matches Sources 1, 2, 3 verbatim.
- **Total tree length = Σ branch lengths** — matches DendroPy "sum of edge lengths" and Biopython "sum of all the branch lengths".
- **Tree height = longest root→leaf path in edges; height(tree) = height(root)** — matches Sources 1, 2 verbatim.
- **Single node → height 0; empty tree → height −1** — matches Sources 1, 2 verbatim.
- **Missing/default branch length → 0** — matches DendroPy.

### Edge-case semantics

- Single leaf → height 0, one leaf: sourced (Sources 1, 2).
- Empty tree → height −1: sourced convention (Sources 1, 2). The null-`PhyloNode` ↔ empty-tree mapping is ASSUMPTION-1 (documented); reasonable and standard.
- All-default (0) branch lengths → length 0: sourced (DendroPy).

### Independent cross-check (Biopython 1.85, executed this session)

| Newick | leaves | count | total_branch_length | height (edges, `unit_branch_lengths=True`) |
|--------|--------|-------|---------------------|---------------------------------------------|
| `((A:1,B:1):1,(C:1,D:1):1);` | A,B,C,D | 4 | **6.0** | **2** |
| `(A:1,(B:1,(C:1,D:1):0.5):0.5);` | A,B,C,D | 4 | **5.0** | **3** (C,D at depth 3) |
| `(C:1,D:1):0.5;` | C,D | 2 | **2.5** | 1 |
| `(A:1,B:1);` | A,B | 2 | **2.0** | **1** |

Every Stage-A / TestSpec expected value is reproduced by an **independent reference implementation** (Biopython), not merely hand-computed or echoed from the repo. Notably, Biopython's `total_branch_length` for `(C:1,D:1):0.5` = **2.5**, confirming the root node's own branch length is counted (M6).

### Findings / divergences

None. Description matches all retrieved sources exactly.

**Stage A verdict: PASS.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs`:
- `GetLeaves` (lines 675–692): null → `yield break`; leaf → yield self; else recurse Left then Right (pre-order). Realises leaf-def + pre-order order.
- `CalculateTreeLength` (lines 707–716): null → 0; else `root.BranchLength + len(Left) + len(Right)`. Sums every node including root → matches Biopython `total_branch_length` (verified 2.5 for the subtree case).
- `GetTreeDepth` (lines 723–748): `EmptyTreeHeight = -1` for null; leaf → 0; else `1 + max(depth(Left), depth(Right))`. Realises edge-count height; single node 0, empty −1.

### Formula realised correctly?

Yes. Each method is a single linear traversal computing exactly the validated formula. Verified branch-by-branch:

| Method | null | leaf | internal |
|--------|------|------|----------|
| GetLeaves | `yield break` (∅) ✓ | yield self ✓ | recurse L,R pre-order ✓ |
| CalculateTreeLength | 0 ✓ | `BranchLength` ✓ | `BL + Σ children` ✓ |
| GetTreeDepth | −1 ✓ | 0 ✓ | `1 + max(L,R)` ✓ |

### Cross-verification table recomputed vs code (full suite run)

All test datasets match the Biopython cross-check above. M6 (`(C:1,D:1):0.5` → 2.5) confirms root-edge inclusion against the external reference.

### Variant/delegate consistency

No `*Fast` variants or instance-property delegates exist for these three methods; canonical methods only. `GetLeaves` is also consumed internally by `RobinsonFouldsDistance` (out of scope) consistently.

### Test quality audit (HARD gate)

- **Sourced expectations, not code echoes:** M4=6.0, M5=5.0, M6=2.5, counts=4, heights 2/3/1 are all independently reproduced by Biopython 1.85 this session. A deliberately-wrong implementation (e.g. omitting the root edge, or counting nodes instead of edges) would fail M6 / M7 / M8.
- **No green-washing:** all assertions use exact `Is.EqualTo` (with `1e-10` numeric tolerance on doubles, appropriate for exact small sums). No ranges, no `Greater`/`AtLeast`, no skipped/ignored tests, no widened tolerances.
- **Cover all logic:** every branch of all three methods is exercised — null (M10), leaf (M3, M9), internal multi-level (M1, M2, C1, M4, M5, M7, M8), root-edge inclusion (M6), default-0 lengths (S1), two-leaf (S2), pre-order ordering (C1), and the empty-tree convention triple (M10). 13 tests, all 13 TestSpec cases covered.
- **Honest green:** full **unfiltered** suite = **Failed: 0, Passed: 6561** (`dotnet test … --no-build`). Build = 0 errors (4 pre-existing NUnit2007 warnings in an unrelated file, `ApproximateMatcher_EditDistance_Tests.cs`; none in PHYLO-STATS-001's file).

**Result of test-quality gate: PASS.**

### Findings / defects

None. No code or test change was required.

**Stage B verdict: PASS.**

## Verdict & follow-ups

- **Stage A: PASS** — description matches Wikipedia (graph theory & ADT), Biopython, and DendroPy verbatim.
- **Stage B: PASS** — implementation realises the validated formulas; tests assert exact externally-sourced values covering every branch.
- **End-state: CLEAN** — no defect found; full suite green (6561 passed, 0 failed).
- **No logged defects.**

## Phase-3 independent re-validation (2026-06-17) — N-ary (multifurcating) refactor (ledger row #11)

### Scope and concern

`PhyloNode` was refactored from the binary `Left`/`Right` model to **N-ary** (commit `c4f0190`):
canonical storage is now an ordered `List<PhyloNode> Children` (`Left`/`Right` retained as accessors
over the first two children). The traversals in this unit — `GetLeaves`, `CalculateTreeLength`,
`GetTreeDepth`, and **`PatristicDistance`** (via `FindMRCAInternal` + `DistanceToTaxon`) — were
rewritten to iterate N children. The load-bearing concern this session: **does every child get
visited, or could a "first-two-children-only" regression have slipped in?**

### Stage A — definitions re-sourced independently this session

Re-fetched the reference docs (not trusted by label):

- **Biopython `Bio.Phylo.BaseTree`** (WebFetch 2026-06-17): `total_branch_length` = "Calculate the sum
  of all the branch lengths in this tree."; `count_terminals` = "Count the number of terminal (leaf)
  nodes within this tree."; `get_terminals` = "Get a list of all of this tree's terminal (leaf)
  nodes."; **`distance()`** = "Calculate the sum of the branch lengths between two targets." — i.e.
  **patristic distance** = sum of branch lengths along the path between two leaves through their MRCA.
  URL: https://biopython.org/docs/latest/api/Bio.Phylo.BaseTree.html
- **DendroPy `Tree.length()`** (WebFetch 2026-06-17): "Returns sum of edge lengths of self. Edges with
  no lengths defined (None) will be considered to have a length of 0."
  URL: https://dendropy.org/library/treemodel.html
- **Wikipedia — Tree (graph theory) / Tree (abstract data type)**: leaf = node with no children;
  height = length of longest downward path to a leaf in **edges**; single node 0; empty tree −1.
  (Biopython's branch-length `depths()` is a *distinct* metric — confirmed not what this unit reports.)

All four definitions (leaf count, total length, edge-height, patristic distance) match the code's
behaviour. **Stage A: PASS.**

### Stage B — polytomy hand-derivation (load-bearing regression guard)

Tree: `((A:0.1,B:0.2,C:0.3):0.0,(D:0.4,E:0.5):0.6);` — internal node **P** (BranchLength 0.0) is a
**3-child polytomy** over leaves A(0.1), B(0.2), **C(0.3)**; node Q(0.6) over D(0.4), E(0.5).

| Quantity | Hand-derived from definition | Code |
|----------|------------------------------|------|
| Leaf count | A,B,C,D,E = **5** | 5 ✓ |
| Total length | 0.1+0.2+0.3+0.0+0.4+0.5+0.6 = **2.1** (incl. 3rd child C's 0.3) | 2.1 ✓ |
| Depth (edges) | root→P→leaf = root→Q→leaf = **2** | 2 ✓ |
| Patristic A→B (MRCA = P) | 0.1 + 0.2 = **0.3** | 0.3 ✓ |
| **Patristic A→C (3rd-child guard, MRCA = P)** | 0.1 + 0.3 = **0.4** | 0.4 ✓ |

The A→C path crosses the polytomy to its **3rd** child (`Children[2]`): a binary-only traversal would
never reach C. Code reviewed: `GetLeaves`/`CalculateTreeLength`/`GetTreeDepth` (lines ~767/796/827)
and `FindMRCAInternal`/`DistanceToTaxon` (lines ~1057/1102) all iterate `foreach (… in
node.Children)` — no `Take(2)`/`Left`+`Right`-only shortcut. **Stage B: PASS.**

### HARD test-quality gate — mutation-verified this session

The polytomy test `TreeStatistics_OverMultifurcatingNode_TraverseAllChildren` lives in
`PhylogeneticAnalyzer_TreeComparison_Tests.cs` (§"Multifurcation (N-ary)" region). To prove it is a
genuine guard and not green-washing, the production source was mutated to `Children.Take(2)` in
`CalculateTreeLength` and `DistanceToTaxon`:

- `CalculateTreeLength` returned **1.8** (drops 3rd child C's 0.3) ≠ 2.1 → **FAIL**.
- `PatristicDistance(A,C)` returned **NaN** (C, the 3rd child, never found) ≠ 0.4 → **FAIL**.

Both assertions failed exactly as required; the source was then restored byte-for-byte (working tree
clean, `git status` empty after restore). Expectations are exact sourced values (`Within 1e-9`), no
ranges, no skips, no widened tolerances. **Test-quality gate: PASS.**

### End-state

**CLEAN** — no code or test defect; all four N-ary traversals visit every child, confirmed by hand
and by mutation. No code/test change required. Full **unfiltered** suite = **Passed: 6779, Failed:
0**; build 0 errors (4 pre-existing unrelated NUnit2007 warnings in `ApproximateMatcher_EditDistance_Tests.cs`).
See FINDINGS_REGISTER §D (PHYLO-STATS-001) and VALIDATION_LEDGER Phase-3 row #11.
