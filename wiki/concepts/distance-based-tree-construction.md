---
type: concept
title: "Distance-based phylogenetic tree construction (UPGMA & Neighbor-Joining — the PHYLO tree-building core)"
tags: [phylogenetics, algorithm]
sources:
  - docs/Evidence/PHYLO-TREE-001-Evidence.md
  - docs/algorithms/Phylogenetics/Newick_Format.md
source_commit: f9b89f022474226460c4ed9124c6062cd89b4815
created: 2026-07-10
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: phylo-tree-001-evidence
      evidence: "Test Unit ID: PHYLO-TREE-001 ... Algorithm: Phylogenetic Tree Construction (UPGMA, Neighbor-Joining)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:evolutionary-distance-matrix
      source: phylo-tree-001-evidence
      evidence: "BuildTreeFromMatrix() accepts pre-computed distance matrices; UPGMA/NJ consume the symmetric n×n distance matrix (PHYLO-DIST-001) as input and turn it into a PhyloNode tree"
      confidence: high
      status: current
---

# Distance-based phylogenetic tree construction (UPGMA & Neighbor-Joining)

The **tree-building core** of the phylogenetics (`PHYLO-*`) family, PHYLO-TREE-001 — the step that
**consumes an [[evolutionary-distance-matrix]]** (the symmetric *n×n* pairwise-distance matrix from
PHYLO-DIST-001) and **produces a `PhyloNode` tree**. Two agglomerative (bottom-up) methods share this
unit: **UPGMA** (Unweighted Pair Group Method with Arithmetic Mean) and **Neighbor-Joining (NJ)**. This
is the piece the rest of the family sits *on top of*: [[phylogenetic-bootstrap-support]] rebuilds a tree
here per replicate, and [[tree-comparison-metrics]] / [[tree-statistics]] compare, query, and summarize
the trees it emits. Validated under test unit **PHYLO-TREE-001**; the literature-traced record is
[[phylo-tree-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade correctness
reference ([[scientific-rigor|research-grade]]), not for clinical use.

The public surface is two methods on `PhylogeneticAnalyzer`: **`BuildTree(sequences, method)`** — the
canonical path that computes distances then builds — and **`BuildTreeFromMatrix(matrix, taxa, method)`**
— which takes a **pre-computed** distance matrix directly (this is what lets the reference Wikipedia
matrices below be tested exactly). Both return a `PhylogeneticTree` carrying `Root`, `Taxa`,
`DistanceMatrix`, and `Method`.

## UPGMA (Unweighted Pair Group Method with Arithmetic Mean)

A simple agglomerative hierarchical clustering method (Sokal & Michener 1958).

- **Output:** a **rooted, ultrametric** tree — all tips are equidistant from the root.
- **Assumption:** a **molecular clock** (constant rate of evolution across lineages).
- **Steps:** (1) each taxon is its own cluster; (2) find the pair of clusters at **minimum distance**;
  (3) merge them, computing new distances as a **weighted (arithmetic-mean) average**; (4) repeat until
  one cluster remains.
- **Branch length:** node **height = distance / 2** (the ultrametric property); the implementation
  tracks cluster heights and emits **incremental** branch lengths (`height_new − height_child`).
- **Complexity:** O(n³) naive, O(n²) optimized.

## Neighbor-Joining (Saitou & Nei 1987)

The workhorse distance method that does **not** assume a clock.

- **Output:** an **unrooted** tree (rooted here by convention — see the final-join rule below).
- **Guarantee:** recovers the **correct topology for additive distance matrices** (INV-N01).
- **Q-matrix step:** compute `Q(i,j) = (n−2)·d(i,j) − Σ_k d(i,k) − Σ_k d(j,k)`, then join the pair
  `(i,j)` with the **minimum Q value** (not the minimum raw distance — NJ's key departure from UPGMA).
- **Branch-length formulas** (from Wikipedia NJ):
  - `δ(f,u) = d(f,g)/2 + (Σ_k d(f,k) − Σ_k d(g,k)) / (2(n−2))`
  - `δ(g,u) = d(f,g) − δ(f,u)`
  - `d(u,k) = (d(f,k) + d(g,k) − d(f,g)) / 2` (distance from the new node *u* to every remaining node)
- **Negative branches preserved:** NJ **may produce negative branch lengths** (INV-N02); the
  implementation **does not clamp** them — the algorithm specification is followed verbatim.
- **Final join = midpoint rooting:** the last join splits the remaining distance **d/2 each**, which
  **preserves all patristic distances** (the additive-matrix guarantee).
- **Complexity:** O(n³).

## Invariants (§3)

**Shared tree structure:**

| ID | Invariant |
|----|-----------|
| INV-01 | Every input taxon appears as a **leaf** |
| INV-02 | **Binary** tree — each internal node has exactly 2 children |
| INV-03 | Number of leaves = *n* (input sequences) |
| INV-04 | Number of internal nodes = *n−1* (binary-tree graph theory) |
| INV-05 | All UPGMA branch lengths ≥ 0 |

**UPGMA-specific:** ultrametric (all tips equidistant from root, INV-U01), rooted (INV-U02),
height = distance/2 (INV-U03).

**NJ-specific:** correct topology for additive matrices (INV-N01), **may produce negative branch
lengths** (INV-N02), no clock assumption (INV-N03).

## Worked oracles (from the Wikipedia reference matrices, §4)

- **UPGMA 5S-rRNA example** (taxa a–e; input matrix `d(a,b)=17`, etc.): clustering order **(a,b)@17 →
  ((a,b),e)@22 → (c,d)@28 → final@33**; branch lengths **δ(a,u)=δ(b,u)=8.5, δ(e,v)=11,
  δ(c,w)=δ(d,w)=14, root height 16.5**; ultrametric — every tip is **16.5** from the root
  (verified to tolerance 1e-10 with exact integer arithmetic).
- **NJ 5-taxon example** (input `d(a,b)=5, d(a,c)=9, …`): first join **(a,b)** with
  **Q₁(a,b) = −50** (minimum Q), branch lengths **δ(a,u)=2, δ(b,u)=3, δ(u,v)=3, δ(c,v)=4, δ(v,w)=2,
  δ(d,w)=2, δ(e,w)=1**; patristic distances recovered from the tree **match the input matrix exactly**
  (additive-matrix guarantee, INV-N01).

## Edge cases and validation (§5)

| Case | Behavior |
|------|----------|
| < 2 sequences | throw (a tree needs ≥ 2 taxa / a pair) |
| Unequal-length sequences | throw (aligned input is a precondition) |
| 2 sequences | trivial binary tree (minimum case) |
| Identical sequences | zero distance → **arbitrary join order** among ties |
| All-zero matrix | all taxa identical |
| All-equal distances | **star topology** resolved arbitrarily |
| Saturated distances (p > 0.75) | inherited from the JC distance step → **+∞** (see [[evolutionary-distance-matrix]] saturation) |
| Single-nucleotide sequences | valid but trivial |
| All-gap columns | zero comparable sites; gap-only columns skipped, identical non-gap sites → distance 0 |

## Newick serialization I/O layer (PHYLO-NEWICK-001)

The `PhyloNode` tree this unit emits is round-tripped to and from **Newick** (New Hampshire) text by the
family's **I/O layer**, PHYLO-NEWICK-001 (`PhylogeneticAnalyzer.ToNewick(PhyloNode, bool)` /
`ParseNewick(string)`, both in `PhylogeneticAnalyzer.cs`). Newick is **grammatical, not biological**: a
tree is recursively nested subtrees with optional labels and branch lengths. This is a **format
serializer, not a separate algorithm** — its literature-traced Evidence record is
[[phylo-newick-001-evidence]]. The `docs/algorithms/Phylogenetics/Newick_Format.md` spec is
reconciled here because the format serializes *this* concept's output rather than defining a new
inference step.

**Core grammar** (Olsen-style, the subset implemented): `Tree → Subtree ";"`,
`Internal → "(" BranchSet ")" Name`, `Branch → Subtree Length`, `Length → empty | ":" number`. Parentheses
delimit internal nodes, commas separate siblings, an optional label follows a subtree, and a branch length
is introduced by `:`.

**Contract.** `ToNewick(node, includeBranchLengths=true)` returns a string; a **`null` node → empty
string**. `ParseNewick(newick)` returns a `PhyloNode`; **null/empty/whitespace-only input throws
`ArgumentException`**. The parser trims the input, strips a trailing `;`, then recurses over `(`, `,`, `)`,
labels, and numeric branch lengths. Both are **O(n)** time, **O(h)** space (h = recursion/parse depth).

**Invariants.**

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `ToNewick` output ends with `;` | serializer appends `;` after the recursive traversal |
| INV-02 | Branch lengths render with `.` as decimal separator | serializer formats via `ToString("F4", CultureInfo.InvariantCulture)` — locale-independent |
| INV-03 | Internal-node names emitted only when they are **valid unquoted labels** | serializer suppresses names containing Newick metacharacters (blanks, `()`, `[]`, `'`, `:`, `;`, `,`) instead of quoting them |
| INV-04 | `ParseNewick` accepts an optional **root branch length** after the main subtree | dedicated post-root `:` handling; parsed into `root.BranchLength` (the Olsen/TreeAlign `:0.0` convention) |

**Serializer/parser asymmetry** (a deliberate limitation): serialization filters *internal* labels to the
conservative unquoted subset above (metachar-bearing names are **omitted**, not quoted), while **leaf
labels are emitted verbatim** (the unquoted check is applied only to internal names). Parsing is more
**permissive** — it reads labels as raw character runs up to the next structural delimiter and does not
enforce the unquoted restriction. The parser also accepts **scientific-notation** branch lengths
(`e`/`E`/`+`/`-`) in addition to digits, decimal points, and signs.

**Documented scope (out-of-scope, not bugs).** Only **binary** trees are supported (`PhyloNode` has just
`Left`/`Right`, so the grammar's multifurcating `BranchSet` is out of scope); **no quoted `'…'` labels**,
**no `[]` comments**, and **no underscore→blank** rewrite. Float precision is ±0.00005 (F4 — adequate for
UPGMA/NJ output; the spec imposes no precision limit). Appropriate for the binary trees this unit's UPGMA
and NJ builders produce, but not for general-purpose Newick interoperability across all phylogenetics
tools. No source contradictions.

## Relationship to the rest of the PHYLO family

This unit is the **hinge** of the family. Upstream, it **consumes** the
[[evolutionary-distance-matrix]] (PHYLO-DIST-001) — `BuildTreeFromMatrix` even takes a pre-built matrix
directly. Downstream, everything else operates on the `PhyloNode` tree it emits:
[[phylogenetic-bootstrap-support]] (PHYLO-BOOT-001) re-runs **this exact distance-matrix → UPGMA/NJ
machinery per resampled replicate** (bootstrap *wraps* tree construction — its oracles use UPGMA +
JukesCantor); [[tree-comparison-metrics]] (PHYLO-COMP-001) compares two such trees (RF / MRCA /
patristic); and [[tree-statistics]] (PHYLO-STATS-001) reads whole-tree summaries off one. The trees are
serialized to/from **Newick** text by the I/O layer, PHYLO-NEWICK-001 ([[phylo-newick-001-evidence]]).

It is **distinct from** [[tumor-phylogeny-clonal-tree-reconstruction]] (ONCO-PHYLO-001), the oncology
CCF-constraint / perfect-phylogeny builder that computes **no distance matrix** and runs **no UPGMA/NJ**,
and from [[phylogenetic-marker-selection]] (PANGEN-MARKER-001), which selects the informative *columns*
a distance matrix would be built over but builds no tree itself. UPGMA and NJ are **alternatives** within
this one unit: UPGMA is rooted/ultrametric and clock-assuming; NJ is unrooted (rooted by midpoint
convention) and clock-free with the additive-topology guarantee.

## Implementation notes and assumptions (§7, §8)

The implementation lives in `PhylogeneticAnalyzer.cs`: `BuildTree()` is canonical (from sequences),
`BuildTreeFromMatrix()` accepts pre-computed matrices for reference-example testing. UPGMA tracks cluster
heights and emits incremental branch lengths; NJ preserves negative branch lengths (no clamping) and
midpoint-roots the final join to preserve patristic distances. The source records **no deviations and no
assumptions** — the implementation strictly follows UPGMA (Sokal & Michener 1958) and NJ (Saitou & Nei
1987) as described in the authoritative sources.

## Reference tools

Definitions trace to **Sokal & Michener (1958)** *"A statistical method for evaluating systematic
relationships"* (University of Kansas Science Bulletin 38:1409–1438 — UPGMA), **Saitou & Nei (1987)**
*"The neighbor-joining method: a new method for reconstructing phylogenetic trees"* (Molecular Biology
and Evolution 4(4):406–425 — NJ, the Q-matrix and branch-length formulas), **Felsenstein (2004)**
*Inferring Phylogenies* (Sinauer), and the Wikipedia articles **UPGMA**, **Neighbor joining**, and
**Phylogenetic tree** (the worked 5S-rRNA and 5-taxon reference matrices). No source contradictions.
