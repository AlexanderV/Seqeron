---
type: concept
title: "Phylogenetic bootstrap support (Felsenstein bootstrap proportions — column resampling → tree replicates → clade support)"
tags: [phylogenetics, algorithm]
mcp_tools:
  - bootstrap_support
sources:
  - docs/Evidence/PHYLO-BOOT-001-Evidence.md
  - docs/algorithms/Phylogenetics/Bootstrap_Analysis.md
source_commit: c9a7fd983fa1cc3da6231edf3531528dd8a35a6f
created: 2026-07-10
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: phylo-boot-001-evidence
      evidence: "Test Unit ID: PHYLO-BOOT-001 ... Algorithm: Phylogenetic Bootstrap Analysis (Felsenstein's Bootstrap Proportions)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:distance-based-tree-construction
      source: phylo-boot-001-evidence
      evidence: "rebuilds a tree from each pseudo-alignment with the same method as the reference tree (here a distance matrix + UPGMA/NJ) — bootstrap re-runs the tree-construction step per replicate"
      confidence: high
      status: current
---

# Phylogenetic bootstrap support (Felsenstein bootstrap proportions)

The **first phylogenetics-family (`PHYLO-*`) unit** and the anchor concept for it: **Felsenstein's
bootstrap** (Felsenstein 1985), the standard way to attach a **confidence value to each clade** of an
inferred phylogenetic tree. Given a multiple-sequence **alignment**, it **resamples the alignment
columns (sites) with replacement** to make many pseudo-alignments, **rebuilds a tree** from each, and
scores every internal clade of the reference tree by the **fraction of replicate trees that contain
that same clade** — the clade's **bootstrap support**. Validated under test unit **PHYLO-BOOT-001**;
the literature-traced record is [[phylo-boot-001-evidence]], [[test-unit-registry]] tracks the unit,
and [[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade
correctness reference ([[scientific-rigor|research-grade]]), not for clinical use.

Bootstrap is the first `PHYLO-*` unit ingested, so this page doubles as the family anchor. The
per-replicate tree here is built by the same [[distance-based-tree-construction|distance-matrix →
UPGMA/NJ]] machinery (PHYLO-TREE-001) the datasets exercise (`UPGMA` + `JukesCantor` distance below) —
bootstrap **wraps and re-runs** that tree-construction step per replicate rather than replacing it. The
pairwise evolutionary distances it consumes are the dedicated [[evolutionary-distance-matrix]] unit
(PHYLO-DIST-001: p-distance / JC69 / K2P / Hamming). The second `PHYLO-*` unit,
[[tree-comparison-metrics]] (PHYLO-COMP-001 — Robinson–Foulds distance, MRCA, patristic distance),
operates on the *same* rooted `PhyloNode` tree but *compares* topologies instead of inferring them:
bootstrap counts clade **agreement** across replicates, Robinson–Foulds counts clade **disagreement**
between two trees — the same split/clade primitive, opposite direction, and no resampling on the RF
side.

## The FBP procedure (Felsenstein 1985 / Lemoine 2018 / Biopython)

Three sources agree on the same four steps:

1. **Resample sites, keep taxa.** "keep all of the original species while sampling characters with
   replacement" (Felsenstein). For an alignment of length *L*, draw *L* column indices uniformly in
   `[0, L-1]` **with replacement** and assemble the columns into a **pseudo-alignment of the same
   length** (Lemoine: "pseudo-alignments of the same length"). Taxa (rows) are never resampled —
   only characters/columns (sites). Biopython's reference code is literally
   `for j in range(length): col = random.randint(0, length - 1)`.
2. **Rebuild a tree** from each pseudo-alignment with the same method as the reference tree (here a
   distance matrix + UPGMA/NJ).
3. **Score clades by presence.** For each **non-terminal** clade of the **reference** (original-data)
   tree, count the replicate trees that contain a clade with the **identical set of terminal (leaf)
   names**. "measure the support of every branch in the reference tree as the proportion of
   pseudo-trees containing that branch" (Lemoine).
4. **Support = count / replicates.** Felsenstein/Lemoine report it as a **percentage** and Biopython
   multiplies by 100 (`c.confidence = (t + 1) * 100.0 / size`); practitioners read thresholds like
   **70%** or **95%** ("group shows up 95% of the time or more → statistically significant").

## What counts as "the same clade"

- **Clade identity is by terminal (leaf) name set** — two clades match iff their leaf-name sets are
  identical; **branch lengths, internal labels, and child order are irrelevant** (Biopython
  `get_support`). A branch of a tree "defines a bipartition of X" (Lemoine); the clade is the cluster
  of taxa on one side.
- **Only non-trivial clades are scored** — single-leaf (terminal) clades are excluded
  (`find_clades(terminal=False)`); trivial one-taxon groups are not bootstrap entities.
- **The reference tree provides the entity set** — support is measured only for clades that appear in
  the original-data tree. Clades that surface **only** in replicates are not reported.
- **Binary per-replicate scoring** — a clade in a replicate either matches exactly (counted) or does
  not; there is **no partial credit**.

## Worked oracles (analytically derived, not copied from a run)

- **Two-group deterministic alignment (support = 1.0).** Taxa A,B,C,D with
  `A=B=AAAAAAAAAA`, `C=D=GGGGGGGGGG` (length 10), method UPGMA + JukesCantor distance, 100 replicates,
  fixed seed 42. Every column is an all-`A` or all-`G` site, so **any** resampled multiset leaves the
  pairwise distances unchanged: `d(A,B)=d(C,D)=0` and the cross-group pairs are JC-saturated
  (`p=1 → +∞`). Every replicate therefore yields the **same** UPGMA topology grouping {A,B} and {C,D},
  so `support({A,B}) = support({C,D}) = 100/100 = 1.0` for **all** seeds.
- **All-identical alignment (degenerate).** Taxa A,B,C each `ACGTACGT` → every replicate produces an
  identical zero-distance matrix → identical reference and bootstrap topology → **every reported clade
  has support 1.0**.

## Invariants and test oracles

- **Support ∈ [0,1]** always (`count/replicates`, count ∈ `[0, replicates]`).
- **Determinism** — same seed + inputs give bit-identical results across runs (resampling is
  RNG-driven; reproducibility requires a fixed seed).
- **Keys = the non-trivial clades of the reference tree** built from the original data — nothing more,
  nothing less.
- **Denominator = replicate count** — every support is `k/replicates` for integer `k`; changing the
  replicate count changes the denominator.
- **Input validation** — null sequences, fewer than 2 sequences, or `replicates < 1` throw
  `ArgumentException`/`ArgumentNullException` (a tree needs ≥2 taxa; <1 replicate is undefined).

## Two documented assumptions (source-consistent)

1. **Rooted-clade scoring rather than unrooted bipartitions.** Felsenstein/Lemoine describe
   bipartitions of **unrooted** trees; this unit scores **rooted clades** (subtree leaf-sets) of the
   UPGMA/NJ tree, matching Biopython's `get_support` (compares clades by terminal set). For rooted
   ultrametric UPGMA trees this is the conventional, consistent representation — a modeling choice,
   not a correctness gap.
2. **Support reported as a proportion in [0,1], not a percentage.** The literature and Biopython
   express *P* as a percent; this unit returns the raw `count/replicates ∈ [0,1]` (multiply by 100 to
   recover the published percentage). A units/labeling choice that changes neither which clades are
   reported nor their ranking.

Thresholds such as 70% or 95% are **descriptive** interpretation conventions, **not parameters** of
the computation.

## Contract, complexity, and method boundary (from the algorithm spec)

The `PhylogeneticAnalyzer.Bootstrap(sequences, replicates, distanceMethod, treeMethod, seed)` entry
point (in `PhylogeneticAnalyzer.cs`) fixes the following defaults, which double as the reproducibility
contract:

| Parameter | Type | Default | Constraint |
|-----------|------|---------|------------|
| `sequences` | `IReadOnlyDictionary<string,string>` | required | non-null, ≥2 entries, equal length |
| `replicates` (B) | `int` | 100 | ≥ 1 |
| `distanceMethod` | `DistanceMethod` | `JukesCantor` | — |
| `treeMethod` | `TreeMethod` | `UPGMA` | UPGMA or NeighborJoining |
| `seed` | `int` | 42 | fixed ⇒ reproducible |

The return is `IReadOnlyDictionary<string,double>` mapping each clade key (sorted, `|`-joined leaf
names) to its support proportion in `[0,1]`. Distance computation uppercases bases and ignores
gaps/ambiguous characters; clade comparison is case-sensitive on taxon names. The `seed` was a later
deviation — added to a previously seed-hardcoded method to enable deterministic tests; the default 42
preserves prior behavior for existing callers.

**Complexity:** `O(B · (n·L + n³))` time, `O(n² + n·L)` space for B replicates, n taxa, L sites —
each replicate resamples L columns (`n·L`) and rebuilds an `O(n³)` distance tree; clade matching is
`O(n)` per internal node. This is fundamentally an `O(B·n³)` procedure (B independent tree builds).

**Explicitly not implemented** (rely on external tooling): (1) **majority-rule consensus tree**
construction from the replicates — this unit returns only the per-clade support map over the reference
tree, with no in-repo consensus-tree builder; (2) **Transfer Bootstrap Expectation (TBE)** — the
gradual transfer-distance support of Lemoine 2018 (use `booster`). The Felsenstein vs transfer contrast:
FBP matches branches by **exact leaf-set identity (binary)** and its support drops with a single
misplaced taxon on large trees, whereas TBE uses a **gradual transfer distance** and is more robust to
rogue taxa. The repository suffix tree was evaluated and is **not used** here — bootstrap does RNG
resampling, distance-matrix tree building, and leaf-set clade matching, with no substring/pattern
search.

## Not the same as tumor-phylogeny reconstruction

[[tumor-phylogeny-clonal-tree-reconstruction]] (ONCO-PHYLO-001) also builds a rooted tree, but it is a
**different, oncology-specific method**: a constraint-satisfaction / perfect-phylogeny builder ordered
by per-sample **cancer-cell-fraction** inequalities (LICHeE/PICTograph) — it computes **no distance
matrix**, runs **no UPGMA/NJ**, and does **no site resampling**. This unit, by contrast, is
**distance-based sequence phylogenetics** and its whole content is the **resampling** confidence
assessment. The two share the word "phylogeny" and nothing of the algorithm. Bootstrap support is a
confidence layer over the distance-matrix tree builders; the upstream sequence-selection step for such
core-genome trees is [[phylogenetic-marker-selection]] (which picks the informative columns a tree —
and its bootstrap — would be built from, but builds no tree and computes no support itself).

## Reference tools

Definitions trace to **Felsenstein, J. (1985)** *"Confidence Limits on Phylogenies: An Approach Using
the Bootstrap"* (Evolution 39(4):783–791 — origin of the method: resample characters with replacement,
same sample size, *P* = fraction of replicates containing a group, 95% ≈ significant), **Lemoine et
al. (2018)** *"Renewing Felsenstein's Phylogenetic Bootstrap in the Era of Big Data"* (Nature 556:452
— restates the formal FBP procedure: same-length pseudo-alignments, branch = bipartition, support =
proportion of pseudo-trees containing the branch), and the **Biopython** `Bio.Phylo.Consensus` module
(`bootstrap`, `bootstrap_trees`, `get_support` — the reference column-resampling + terminal-set
clade-matching implementation). No source contradictions: the three describe the identical procedure,
differing only in unrooted-bipartition vs rooted-clade phrasing (the assumption reconciled above).
