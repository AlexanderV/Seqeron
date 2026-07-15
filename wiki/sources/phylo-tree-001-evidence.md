---
type: source
title: "Evidence: PHYLO-TREE-001 (Phylogenetic Tree Construction — UPGMA & Neighbor-Joining)"
tags: [validation, phylogenetics]
doc_path: docs/Evidence/PHYLO-TREE-001-Evidence.md
sources:
  - docs/Evidence/PHYLO-TREE-001-Evidence.md
source_commit: 6e5907cbbe510e53e4afa483dd33991739bb93fa
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PHYLO-TREE-001

The validation-evidence artifact for test unit **PHYLO-TREE-001** — **Phylogenetic Tree Construction**:
the distance-based tree-building step that **consumes a distance matrix and emits a `PhyloNode` tree**
via two agglomerative methods, **UPGMA** (rooted, ultrametric, clock-assuming; height = d/2) and
**Neighbor-Joining** (unrooted/midpoint-rooted, clock-free, Q-matrix-driven, additive-topology
guarantee, negative branches preserved). This is the **tree-building core** the phylogenetics family
sits on; the definitions, formulas, invariants, worked oracles, and edge cases are synthesized in the
dedicated concept [[distance-based-tree-construction]]. One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; see [[test-unit-registry]] for how units
are tracked.

## What this file records

- **Authoritative sources (§1):** Wikipedia **UPGMA**, **Neighbor joining**, **Phylogenetic tree** +
  **Saitou & Nei (1987)** (MBE 4(4):406–425 — NJ) + **Sokal & Michener (1958)** (Univ. Kansas Sci.
  Bull. 38:1409–1438 — UPGMA) + **Felsenstein (2004)** *Inferring Phylogenies*.
- **Algorithm definitions (§2):** UPGMA — nearest-pair merge, arithmetic-mean distance update, rooted
  ultrametric, **height = d/2**, O(n³)/O(n²). NJ — **Q(i,j) = (n−2)·d(i,j) − Σd(i,k) − Σd(j,k)**, join
  minimum-Q pair, branch lengths `δ(f,u)=d(f,g)/2 + (Σd(f,k)−Σd(g,k))/(2(n−2))`, `δ(g,u)=d(f,g)−δ(f,u)`,
  `d(u,k)=(d(f,k)+d(g,k)−d(f,g))/2`; O(n³); correct topology for additive matrices.
- **Invariants (§3):** shared — all taxa are leaves, binary tree, *n* leaves / *n−1* internal nodes,
  UPGMA branches ≥ 0. UPGMA — ultrametric (tips equidistant from root), rooted, height = d/2. NJ —
  additive-matrix topology correct, **may produce negative branch lengths**, no clock.
- **Worked oracles (§4):** UPGMA 5S-rRNA matrix → clustering order (a,b)@17 → ((a,b),e)@22 → (c,d)@28 →
  final@33, δ(a,u)=8.5 / δ(e,v)=11 / δ(c,w)=14 / root 16.5, all tips 16.5 from root (tol 1e-10). NJ
  5-taxon matrix → first join (a,b) with **Q₁(a,b)=−50**, δ(a,u)=2 / δ(b,u)=3 / δ(u,v)=3 / δ(c,v)=4 /
  δ(v,w)=2 / δ(d,w)=2 / δ(e,w)=1, patristic distances recover the input matrix exactly.
- **Edge cases (§5):** <2 seqs / unequal lengths → throw; 2 seqs → trivial tree; identical → zero
  distance, arbitrary join order; all-zero matrix → identical taxa; all-equal → star topology; saturated
  p>0.75 → JC +∞; single-nucleotide → trivial; all-gap columns → zero comparable sites.
- **Methodology (§6):** 32 test runs (M01–M13, S01–S05, C01–C05) verifying the Wikipedia reference
  examples, ultrametric property, non-negative UPGMA branches, binary structure, taxa preservation,
  validation errors, performance (50 seqs × 100bp UPGMA < 30s), gap handling. **S06 removed** (duplicate
  of S02). `BuildTreeFromMatrix` public API enables direct testing against known reference matrices.
- **Implementation (§7):** `PhylogeneticAnalyzer.cs` — `BuildTree()` canonical, `BuildTreeFromMatrix()`
  for pre-computed matrices; UPGMA incremental branch lengths (height_new − height_child); NJ preserves
  negative branches (no clamping) and midpoint-roots the final join; returns `PhylogeneticTree` with
  Root / Taxa / DistanceMatrix / Method.

## Deviations and assumptions

**None (§8).** The implementation strictly follows UPGMA (Sokal & Michener 1958) and NJ (Saitou & Nei
1987) as described in the authoritative sources. No source contradictions.
