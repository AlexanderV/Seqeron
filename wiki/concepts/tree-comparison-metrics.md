---
type: concept
title: "Tree comparison metrics (Robinson–Foulds distance, MRCA, patristic distance)"
tags: [phylogenetics, algorithm]
sources:
  - docs/Evidence/PHYLO-COMP-001-Evidence.md
source_commit: 3f492b584e4bfe5aee958659ec2f15a8fabed25a
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: phylo-comp-001-evidence
      evidence: "Test Unit ID: PHYLO-COMP-001 ... Algorithm: Tree Comparison (Robinson-Foulds Distance, MRCA, Patristic Distance)"
      confidence: high
      status: current
---

# Tree comparison metrics (Robinson–Foulds distance, MRCA, patristic distance)

The **second phylogenetics-family (`PHYLO-*`) unit**, PHYLO-COMP-001 — three distinct **operations on an
already-built phylogenetic tree** that compare or query it rather than infer it:
**Robinson–Foulds (RF) distance** (how topologically different two trees are), **most recent common
ancestor (MRCA)** (the deepest node ancestral to a taxon pair), and **patristic distance** (branch-length
path between two taxa through their MRCA). This is genuinely separate from
[[phylogenetic-bootstrap-support]] (PHYLO-BOOT-001): bootstrap attaches a *confidence* to a single tree by
resampling; this unit *compares topologies* and *reads distances off* an existing tree — no resampling, no
tree building. Validated under test unit **PHYLO-COMP-001**; the literature-traced record is
[[phylo-comp-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade correctness
reference ([[scientific-rigor|research-grade]]), not for clinical use.

All three operate on the same rooted binary `PhyloNode` (`Left`/`Right` children) that the
UPGMA/NJ builder produces from an [[evolutionary-distance-matrix]] (PHYLO-DIST-001, the family's
pairwise-distance substrate) — the same machinery [[phylogenetic-bootstrap-support]] resamples; the
scope is therefore **rooted, bifurcating trees** (§7 resolved design decisions).

## 1. Robinson–Foulds distance (`RobinsonFouldsDistance`)

The **symmetric-difference metric** (Robinson & Foulds 1981): `RF = A + B`, where **A** is the number of
splits (bipartitions) implied by the first tree but not the second and **B** the reverse. A split is the
partition of taxa induced by **removing one internal edge**; each tree's split set is compared as a set.

- **Implementation:** extract splits (bipartitions) from each tree, represent each **canonically** (using
  the smaller partition side), then count the **symmetric difference** of the two split sets.
- **Rooted-clade approach:** this unit computes RF on **rooted** trees using clades (subtree leaf-sets),
  equivalent to the Wikipedia dummy-leaf construction ("Rooted trees can be examined by attaching a dummy
  leaf to the root node"). A rooted binary tree with *n* taxa has *n−2* non-trivial clades.
- **Raw count, not normalized:** returns the Robinson & Foulds (1981) raw count — not divided by 2, not
  scaled to [0,1] (some software does both; the original definition does not).
- **Properties:** a **proper metric** (identity, symmetry, triangle inequality); linear-time algorithm
  exists (Day 1985).

### Max-RF formulas (rooted vs unrooted — a genuine difference)

- **Rooted** binary, *n* taxa: **max RF = 2(n−2)** (each tree has *n−2* non-trivial clades).
- **Unrooted** binary, *n* taxa: **max RF = 2(n−3)**.
- The two are reconciled by the dummy-leaf equivalence: an unrooted tree on *n+1* leaves gives
  `2((n+1)−3) = 2(n−2)`, matching the rooted formula.

## 2. Most recent common ancestor (`FindMRCA`)

The deepest node that is an ancestor of **both** queried taxa (a.k.a. LCA / concestor). In a rooted tree
every node is the MRCA of its descendant leaves.

- **Implementation:** recursive tree traversal; returns the node when **both** taxa are found in its
  subtrees; leaf nodes matched **by name**.
- Complexity **O(n)** via one recursive pass.
- `MRCA(x, x)` = the taxon node itself; the root is the MRCA of all taxa.

## 3. Patristic distance (`PatristicDistance`)

The **sum of branch lengths** along the path connecting two taxa through their MRCA:
`PD(x, y) = dist(x → MRCA) + dist(MRCA → y)`. It reflects evolutionary divergence along the tree and
**requires meaningful branch lengths** (zero-length branches → distance 0 even between distinct taxa).

- **Implementation:** uses `FindMRCA` to locate the common ancestor, then sums branch lengths from the
  MRCA down to each taxon.

## Invariants and test oracles

| Metric | Invariants (from §5) |
|--------|----------------------|
| **RF** | `RF(T,T)=0` (identity) · `RF(T1,T2)=RF(T2,T1)` (symmetry) · `RF ≥ 0` · **RF is even** (symmetric difference of two sets) |
| **MRCA** | `MRCA(x,x)` = node for x · `MRCA(x,y)=MRCA(y,x)` · always an ancestor of both · unique per (tree, pair) |
| **Patristic** | `PD(x,x)=0` · `PD(x,y)=PD(y,x)` · `PD ≥ 0` · `PD(x,y)=dist(x,MRCA)+dist(y,MRCA)` |

Behavioural oracles (§3): RF identical trees → 0, completely different topologies → maximum, shared splits
→ symmetric-difference count; MRCA two siblings → their parent, distant taxa → deepest common ancestor,
root is MRCA of all taxa; patristic same taxon → 0, siblings → sum of branch lengths to parent.

## Edge cases and API contract (§4)

- **RF:** null / single-taxon / two-taxon trees → 0 (no internal edges to differ on); star topology vs
  binary tree → different split counts.
- **MRCA:** null root → null; **taxon not in tree → null**; empty taxa list is
  implementation-dependent.
- **Patristic:** **taxon not in tree → NaN**; zero branch lengths → 0 even for different taxa; a
  single-taxon tree admits only the same-taxon query (→ 0).

## Two documented scope decisions (source-backed, §7)

1. **Binary trees only** — `PhyloNode` has `Left`/`Right` by design; multifurcating trees are out of
   scope (Wikipedia notes both bifurcating and multifurcating trees exist; scope is limited to
   bifurcating, the standard phylogenetics case).
2. **Rooted trees only** — UPGMA and NJ here produce rooted trees; RF is computed via the rooted-clade
   approach, equivalent to Wikipedia's dummy-leaf-at-root construction for unrooted trees.

Both are modeling/scope choices, not correctness gaps; the source records **no deviations** and states no
assumptions remain unresolved.

## Relationship to the rest of the PHYLO family

RF's split/bipartition comparison is the same primitive that [[phylogenetic-bootstrap-support]] uses when
it scores a reference clade by presence across replicate trees — bootstrap counts clade agreement,
RF counts clade *disagreement* — but the two are separate units: this one is a deterministic tree-vs-tree
(or tree-query) computation with **no resampling and no tree inference**. Distinct again from
[[tumor-phylogeny-clonal-tree-reconstruction]] (ONCO-PHYLO-001), which builds an oncology-specific tree
from cancer-cell-fraction constraints and computes none of these metrics.

## Reference tools

Definitions trace to **Robinson & Foulds (1981)** *"Comparison of phylogenetic trees"*
(Mathematical Biosciences 53:131–147, doi:10.1016/0025-5564(81)90043-2 — RF = symmetric difference of
split sets, a proper metric), **Day (1985)** (linear-time RF), **Smith (2020)**
*"Information theoretic Generalized Robinson–Foulds metrics"* (Bioinformatics, doi:10.1093/bioinformatics/btaa614),
and the Wikipedia articles on the **Robinson–Foulds metric**, **Most recent common ancestor**, and
**Phylogenetic tree**. No source contradictions.
