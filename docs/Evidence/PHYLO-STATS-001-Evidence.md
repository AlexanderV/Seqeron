# Evidence Artifact: PHYLO-STATS-001

**Test Unit ID:** PHYLO-STATS-001
**Algorithm:** Tree Statistics (leaves, total tree length, tree height/depth)
**Date Collected:** 2026-06-13

---

## Online Sources

### Tree (graph theory) — Wikipedia (cites primary graph-theory definitions)

**URL:** https://en.wikipedia.org/wiki/Tree_(graph_theory)
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 4 (encyclopedia of standard graph-theory definitions)

**Key Extracted Points:**

1. **Leaf:** "A _leaf_ is a vertex with no children." (also: "an external vertex … terminal vertex or leaf … is a vertex of degree 1.")
2. **Depth:** "The _depth_ of a vertex is the length of the path to its root (root path)." "The root has depth zero."
3. **Height:** "The _height_ of a vertex in a rooted tree is the length of the longest downward path to a leaf from that vertex. The _height_ of the tree is the height of the root."
4. **Single-node tree:** "a tree with only a single vertex (hence both a root and leaf) has depth and height zero."
5. **Empty tree:** "Conventionally, an empty tree (a tree with no vertices, if such are allowed) has depth and height −1."

### Tree (abstract data type) — Wikipedia

**URL:** https://en.wikipedia.org/wiki/Tree_(abstract_data_type)
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 4

**Key Extracted Points:**

1. **Leaf / external node:** "An external node (also known as an outer node, leaf node, or terminal node) is any node that does not have child nodes."
2. **Depth:** "The depth of a node is the length of the path to its root (i.e., its root path)." The root has depth zero.
3. **Height:** "The height of a node is the length of the longest downward path to a leaf from that node. The height of the root is the height of the tree." "Leaf nodes have height zero, and a tree with only a single node … has depth and height zero."
4. **Empty tree:** "Conventionally, an empty tree (tree with no nodes, if such are allowed) has height −1."

### Biopython `Bio.Phylo.BaseTree` API documentation (reference implementation)

**URL:** https://biopython.org/docs/latest/api/Bio.Phylo.BaseTree.html
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (established bioinformatics library)

**Key Extracted Points:**

1. **get_terminals:** "Get a list of all of this tree's terminal (leaf) nodes."
2. **count_terminals:** "Count the number of terminal (leaf) nodes within this tree."
3. **is_terminal:** "Check if this is a terminal (leaf) node." (leaf = node with no children)
4. **total_branch_length:** "Calculate the sum of all the branch lengths in this tree."

### DendroPy `Tree` model documentation (reference implementation)

**URL:** https://dendropy.org/library/treemodel.html
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (established phylogenetics library)

**Key Extracted Points:**

1. **Tree.length():** "Returns sum of edge lengths of self. Edges with no lengths defined (None) will be considered to have a length of 0." — confirms tree length = sum of branch lengths, and the "missing length → 0" convention.

### Minimum evolution — Wikipedia (context for tree length as an optimized quantity)

**URL:** https://en.wikipedia.org/wiki/Minimum_evolution
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 4

**Key Extracted Points:**

1. **Total branch length:** Minimum evolution "selects the branching pattern with the smallest total branch length"; the modern framework is attributed to Rzhetsky & Nei (1992). This establishes "total tree length = sum of branch lengths" as the quantity used by the ME criterion. (No standalone formula L = Σ branch lengths is quoted in the intro; the explicit "sum of edge lengths" wording comes from DendroPy/Biopython above.)

---

## Documented Corner Cases and Failure Modes

### From Tree (graph theory) / Tree (abstract data type)

1. **Single-node tree:** root that is also a leaf has height and depth 0.
2. **Empty tree:** by convention has height −1.

### From DendroPy

1. **Undefined branch length:** edges with no length are treated as length 0 when summing tree length.

---

## Test Datasets

### Dataset: Balanced 4-taxon tree `((A:1,B:1):1,(C:1,D:1):1)`

**Source:** Constructed from the cited tree-statistics definitions (graph-theory height/leaf; DendroPy/Biopython tree length). Worked by hand below.

| Parameter | Value |
|-----------|-------|
| Topology | `((A,B),(C,D))` rooted, fully bifurcating |
| Branch lengths | A=1, B=1, C=1, D=1, internal(AB)=1, internal(CD)=1 |
| Leaves | A, B, C, D → count = 4 |
| Total tree length | 1+1+1+1+1+1 = 6 |
| Height (edges to deepest leaf) | root→internal→leaf = 2 |

### Dataset: Caterpillar (ladder) tree `(A:1,(B:1,(C:1,D:1)))`

**Source:** Constructed from the cited definitions; worked by hand.

| Parameter | Value |
|-----------|-------|
| Topology | `(A,(B,(C,D)))` rooted |
| Leaves | A, B, C, D → count = 4 |
| Branch lengths (leaves) | A=1, B=1, C=1, D=1 |
| Branch lengths (internal) | (CD)=0.5, (B,(CD))=0.5 |
| Total tree length | 1+1+1+1+0.5+0.5 = 5 |
| Height (edges to C or D) | root→·→·→leaf = 3 |

### Dataset: Single-leaf tree and empty tree

**Source:** Graph-theory / ADT conventions (single node height 0; empty tree height −1).

| Parameter | Value |
|-----------|-------|
| Single leaf `A` | leaves = {A} (count 1); length = A.BranchLength; height = 0 |
| Null root | leaves = ∅; length = 0; height = −1 |

---

## Assumptions

1. **ASSUMPTION: Null root maps to the "empty tree" convention.** The cited sources define the *empty tree* (no vertices) to have height −1. A `null` `PhyloNode` reference is the repository's representation of "no tree", so `GetTreeDepth(null)` returns −1 and `GetLeaves(null)`/`CalculateTreeLength(null)` return the empty/zero result. This is the natural mapping of the graph-theory convention onto the nullable node model; no source addresses a C# null reference directly.

---

## Recommendations for Test Coverage

1. **MUST Test:** `GetLeaves` returns exactly the terminal nodes (no children), in order, for a multi-leaf tree — Evidence: graph-theory leaf definition; Biopython `get_terminals`.
2. **MUST Test:** `CalculateTreeLength` equals the sum of every branch length (balanced tree = 6, caterpillar = 5) — Evidence: DendroPy `length()`, Biopython `total_branch_length`.
3. **MUST Test:** `GetTreeDepth` equals the number of edges to the deepest leaf (balanced = 2, caterpillar = 3, single leaf = 0) — Evidence: graph-theory/ADT height definition.
4. **MUST Test:** Single-leaf tree → height 0, one leaf, length = its branch length — Evidence: "single node has height zero".
5. **MUST Test:** Null root → no leaves, length 0, height −1 — Evidence: empty-tree convention.
6. **SHOULD Test:** Tree length treats default (0) branch lengths as 0 — Rationale: DendroPy "missing length → 0".
7. **COULD Test:** Leaf count via `GetLeaves(...).Count()` matches Biopython `count_terminals` semantics — Rationale: cross-library agreement.

---

## References

1. Wikipedia. 2026. Tree (graph theory). https://en.wikipedia.org/wiki/Tree_(graph_theory)
2. Wikipedia. 2026. Tree (abstract data type). https://en.wikipedia.org/wiki/Tree_(abstract_data_type)
3. Cock PJA et al. Biopython `Bio.Phylo.BaseTree` (get_terminals, count_terminals, is_terminal, total_branch_length). https://biopython.org/docs/latest/api/Bio.Phylo.BaseTree.html
4. Sukumaran J, Holder MT. DendroPy `Tree` model (`Tree.length`). https://dendropy.org/library/treemodel.html
5. Rzhetsky A, Nei M. 1992. A simple method for estimating and testing minimum-evolution trees. (cited via) Wikipedia. Minimum evolution. https://en.wikipedia.org/wiki/Minimum_evolution

---

## Change History

- **2026-06-13**: Initial documentation.
