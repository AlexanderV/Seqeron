---
type: source
title: "Evidence: PHYLO-COMP-001 (Tree Comparison — Robinson–Foulds distance, MRCA, patristic distance)"
tags: [validation, phylogenetics]
doc_path: docs/Evidence/PHYLO-COMP-001-Evidence.md
sources:
  - docs/Evidence/PHYLO-COMP-001-Evidence.md
source_commit: 3f492b584e4bfe5aee958659ec2f15a8fabed25a
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PHYLO-COMP-001

The validation-evidence artifact for test unit **PHYLO-COMP-001** — **Tree Comparison**: three
operations on an already-built rooted binary phylogenetic tree — **Robinson–Foulds (RF) distance**
(symmetric difference of the two trees' split/bipartition sets), **MRCA** (deepest common ancestor of a
taxon pair), and **patristic distance** (branch-length path between two taxa through their MRCA). This is
the **second phylogenetics-family (`PHYLO-*`) Evidence file** (after PHYLO-BOOT-001) and one instance of
the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the definitions,
formulas, invariants, and scope decisions are synthesized in the dedicated concept
[[tree-comparison-metrics]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:** Wikipedia **Robinson–Foulds metric** (RF = `A + B`, splits implied by one tree but
  not the other; rooted trees examined by attaching a dummy leaf to the root; some tools divide by 2 or
  scale to max 1) + Wikipedia **Most recent common ancestor** + Wikipedia **Phylogenetic tree**;
  **Robinson & Foulds (1981)** original paper (doi:10.1016/0025-5564(81)90043-2 — RF is a proper metric,
  computed via splits from edge removal); **Smith (2020)** information-theoretic generalized RF metrics
  (doi:10.1093/bioinformatics/btaa614). Day (1985) cited for the linear-time RF algorithm.
- **Key properties:** RF is a proper metric (identity/symmetry/triangle inequality) and **even**; raw
  count returned (not normalized). Max RF = **2(n−2)** for **rooted** binary trees vs **2(n−3)** for
  **unrooted** — reconciled by dummy-leaf equivalence (unrooted on n+1 leaves → 2((n+1)−3)=2(n−2)). A
  rooted binary tree has n−2 non-trivial clades. MRCA found in O(n) by recursive traversal (leaf
  name-matching); patristic distance uses `FindMRCA` then sums branch lengths from MRCA to each taxon.
- **Test cases (§3) / corner cases (§4):** RF identical→0, completely-different→max, shared-splits→
  symmetric-difference count; null/single-taxon/two-taxon trees→0; star vs binary→differing split counts.
  MRCA siblings→parent, `MRCA(x,x)`→x, distant→deepest ancestor, root→MRCA of all; null root or
  **taxon-not-in-tree → null**. Patristic same-taxon→0, siblings→sum-to-parent, distant→sum via MRCA;
  **taxon-not-in-tree → NaN**; zero branch lengths → 0 even for distinct taxa.
- **Invariants (§5):** RF `RF(T,T)=0` / symmetric / `≥0` / even; MRCA symmetric / always-an-ancestor /
  unique; patristic `PD(x,x)=0` / symmetric / `≥0` / `= dist(x,MRCA)+dist(y,MRCA)`.

## Deviations and assumptions

No deviations from the literature; **all §7 design decisions resolved with external sources, no
assumptions remain**. Two source-backed **scope** decisions: (1) **binary trees only** — `PhyloNode` has
`Left`/`Right` children by design; multifurcating trees out of scope (the standard phylogenetics case).
(2) **Rooted trees only** — UPGMA/NJ produce rooted trees; RF is computed via the rooted-clade approach,
equivalent to Wikipedia's dummy-leaf-at-root construction. Raw RF count returned (consistent with Robinson
& Foulds 1981, not the divide-by-2 / scale-to-1 variants). No source contradictions.
