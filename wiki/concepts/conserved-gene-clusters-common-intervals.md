---
type: concept
title: "Conserved gene clusters (common intervals)"
tags: [comparative-genomics, algorithm]
sources:
  - docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md
  - docs/algorithms/Comparative_Genomics/Conserved_Gene_Clusters.md
  - docs/Validation/reports/COMPGEN-CLUSTER-001.md
source_commit: 665dc3361ce2789ca8ede9ad2e88ea718c20310e
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: compgen-cluster-001-evidence
      evidence: "Test Unit ID: COMPGEN-CLUSTER-001 ... Algorithm: Conserved Gene Clusters (common intervals of permutations)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:synteny-and-rearrangement-detection
      source: compgen-cluster-001-evidence
      evidence: "Conserved gene clusters (common intervals, COMPGEN-CLUSTER-001) and synteny (COMPGEN-SYNTENY) are sibling Comparative-genomics units — clusters ask for a gene SET contiguous in every genome (order-free), synteny for a collinear ordered block"
      confidence: medium
      status: current
---

# Conserved gene clusters (common intervals)

A **conserved gene cluster** is a set of genes that appears as a **contiguous block in
every genome** being compared, regardless of the internal gene order or orientation inside
that block. Formally this is the **common interval** model of permutations (Uno & Yagiura
2000; Heber & Stoye 2001). It is the sixth-ingested unit of the **Comparative-genomics**
family (`COMPGEN-*`) and a sibling of [[average-nucleotide-identity]] and the shared synteny
anchor [[synteny-and-rearrangement-detection]]: where ANI measures *nucleotide* identity and
synteny requires a *collinear ordered* block, a common interval only requires the same **gene
set** to be contiguous in each genome — order and strand inside the window are free. The
end-to-end pipeline [[genome-comparison-core-dispensable]] composes ortholog detection and
synteny into a core/dispensable gene partition. Validated
under test unit **COMPGEN-CLUSTER-001**; the validation record is
[[compgen-cluster-001-evidence]], the independent two-stage re-validation verdict is
[[compgen-cluster-001-report]] (Stage A PASS / Stage B PASS-WITH-NOTES / End state CLEAN — no code
defect; three weak test assertions strengthened to exact sourced sets), [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

## The common-interval model

Map each gene to its **ortholog-group label**, so each genome becomes a permutation (or, with
paralogs, a sequence) of labels. The definitions trace to the literature:

- **Interval of a permutation** (`P_k`): the set of elements between positions `i` and `j`
  (inclusive), defined only for `1 ≤ i < j ≤ n` — so intervals have **size ≥ 2**; singletons
  are excluded by the `i < j` constraint (Bui-Xuan, Habib & Paul 2013, §2).
- **Common interval** (Uno & Yagiura 2000, Definition 1): a set of integers that is an
  interval of **every** `P_k`, `k ∈ [K]` — i.e. a label set that is contiguous in *all* K
  genomes simultaneously.
- **k-permutation generalisation** (Heber & Stoye 2001): the model extends to a family of `k`
  permutations; a cluster is a set contiguous in **all** `k` genomes. Optimal `O(kn + z)` time
  / `O(n)` space for `z` common intervals.
- **Sequence extension with duplicates** (Didier et al. 2013): when ortholog groups recur
  (paralogs/duplications) a genome is a *sequence*, not a permutation. A set `I` is a common
  interval iff **some** contiguous window in each genome has exactly the label set `I` (any
  *location* suffices). This is the model that handles repeated genes.

## Algorithm and parameters

Seqeron finds all conserved clusters by the simple strict-model check — for each candidate
label set, verify it forms a contiguous window (exact label set, no foreign groups inside) in
every genome. This is the `O(n² · K_genomes)` variant, adequate for the small gene-cluster
inputs in scope; the output-sensitive `O(n + K)` (Uno & Yagiura RC) and `O(kn + z)` (Heber &
Stoye) algorithms are the reference upper bounds, not what is implemented.

| Parameter | Meaning |
|-----------|---------|
| `minClusterSize` | discard clusters smaller than this threshold (filters trivial short clusters) |
| `maxGap` | retained for **API/MCP backward compatibility only** — does **not** relax the strict model (see below) |

Invariants / behaviour:

- The **whole set** `(1..n)` is always a common interval of any family (trivial).
- A common interval is defined over a **family of K ≥ 2** genomes — with a single genome
  every interval is trivially "common", so the conserved-cluster question is meaningful only
  for **K ≥ 2**.
- Result is a **set** of clusters — order-independent, deterministic, reproducible.

## Strict (gap-free) model — the `maxGap` assumption

The single documented **assumption** is API-shape only, not a correctness gap: the public
method keeps a `maxGap` parameter, but the validated/tested behaviour is the **strict,
gap-free common-interval model** (Uno & Yagiura 2000; Heber & Stoye 2001) — a cluster must
occupy a *contiguous* window in every genome (its group set equals the set of all groups in
some window, with no foreign groups inside). `maxGap` does **not** relax this in the validated
path. The **gene-teams** gapped extension (Bergeron, Corteel & Raffinot 2002) is **not**
implemented (its source was not retrievable in the evidence session). The correctness contract
is the strict common-interval definition.

## Documented oracles

- **Two-permutation golden vector** (Bui-Xuan, Habib & Paul, Example 1): genome 1 = `1 2 3 4 5
  6 7` (Id₇), genome 2 = `7 2 1 3 6 4 5`. All **non-trivial** common intervals: `{1,2}`,
  `{1,2,3}`, `{3,4,5,6}`, `{4,5}`, `{4,5,6}`, `{1,2,3,4,5,6}`; plus the trivial whole set
  `{1,…,7}`. (Independently reproduced by brute force over all subsets of size ≥ 2.)
- **Split-in-one-genome negative** — a set contiguous in genome 1 but split in genome 2 is
  **not** a cluster. E.g. `{2,3}`: contiguous in Id₇ but positions of 2 and 3 in genome 2 are
  2 and 4 (non-adjacent) → not common. (Note `{1,2}` *is* common: 1,2 sit at positions 3,2 in
  genome 2 → adjacent.)
- **Sequence with duplicates** (Didier et al., Example 1): `T = 1 2 5 2 1 4 3 1 2 6 5`,
  `S = 5 6 4 2 3 4 1 5`. `{1,2}` is an interval of `T` but **not** a common interval (not an
  interval of `S`); `{1,2,3,4}` **is** a common interval (a contiguous window with that label
  set exists in both, five locations on `T`, two on `S`).

## Edge cases

- **`minClusterSize`** filters out clusters below the threshold.
- **Fewer than two genomes** → no clusters (common interval undefined for K < 2).
- **Identical gene order across all genomes** → every window of size ≥ `minClusterSize` is
  conserved.
- **Repeated ortholog-group labels** (paralogs) → a set is reported iff a contiguous window in
  each genome has exactly that label set (any location suffices).

## Reference tools

The definitions trace to **Uno & Yagiura 2000** (Algorithmica, the originating common-interval
model, `O(n²)` LHP + `O(n+K)` RC algorithms), **Heber & Stoye 2001** (CPM, the k-permutation
generalisation, `O(kn+z)`), **Bui-Xuan, Habib & Paul 2013** (MinMax-Profiles, the unifying
view + golden-vector Example 1), and **Didier et al. 2013** (extension from permutations to
sequences with duplicates). No deviations from the sources are recorded; the one assumption is
the strict-model / `maxGap`-inert API-compat note above.
