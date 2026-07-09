---
type: source
title: "Evidence: COMPGEN-CLUSTER-001 (Conserved gene clusters — common intervals)"
tags: [validation, comparative-genomics]
doc_path: docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md
sources:
  - docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md
source_commit: 2ef49f21aacb36551908473c37944c54ebe55323
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: COMPGEN-CLUSTER-001

The validation-evidence artifact for test unit **COMPGEN-CLUSTER-001** — Conserved Gene
Clusters under the **common-interval** model of permutations: a gene set that is contiguous in
*every* genome, order- and strand-free inside the window. This is a **Comparative-genomics**
family Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm, its parameters,
invariants, worked oracles, and corner cases are summarized in
[[conserved-gene-clusters-common-intervals]]. Its sibling COMPGEN units are
[[average-nucleotide-identity]] and the shared synteny anchor
[[synteny-and-rearrangement-detection]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (all authority rank 1):
  - **Bui-Xuan, Habib & Paul 2013** (MinMax-Profiles, arXiv:1304.5140) — interval `[i,j]`
    defined only for `1 ≤ i < j ≤ n` (size ≥ 2, singletons excluded); common interval
    (Definition 1, attributed to Uno & Yagiura) = a set that is an interval of each `P_k`;
    the whole set `(1..n)` is always common. Worked **Example 1** (verbatim golden vector).
  - **Uno & Yagiura 2000** (Algorithmica 26(2):290–309) — the originating common-interval
    model: "a pair of intervals of these permutations consisting of the same set of elements".
    Complexity `O(n²)` (LHP) + output-sensitive `O(n+K)` (RC). (Full text behind a Springer
    auth redirect; only citation + abstract facts used.)
  - **Didier, Schmidt, Stoye & Tsur 2013** (arXiv:1310.4290) — extension from permutations to
    **sequences with duplicates**: a common interval is a label set that is an interval (some
    contiguous window / *location*) of both sequences; handles paralogs. Worked Example 1.
  - **Heber & Stoye 2001** (CPM, LNCS 2089:207–218) — the k-permutation generalisation: all
    `z` common intervals of `k` permutations in optimal `O(kn+z)` time, `O(n)` space.
- **Algorithm behaviour (from the artifact):** map genes to ortholog-group labels; a cluster
  is a label set contiguous in **every** genome (K ≥ 2). The repository implements the simple
  `O(n² · K_genomes)` strict check (adequate for small gene-cluster inputs). `minClusterSize`
  filters short clusters; `maxGap` is retained for API/MCP back-compat only.
- **Datasets (documented oracles):**
  - *Two-permutation golden vector* — genome 1 `1 2 3 4 5 6 7` (Id₇), genome 2 `7 2 1 3 6 4 5`:
    non-trivial common intervals `{1,2}`, `{1,2,3}`, `{3,4,5,6}`, `{4,5}`, `{4,5,6}`,
    `{1,2,3,4,5,6}`; trivial whole set `{1,…,7}`. (Brute-force recomputation matches the paper.)
    Split-negative: `{2,3}` (positions 2,4 in genome 2 → non-adjacent) is not common.
  - *Sequence with duplicates* — `T = 1 2 5 2 1 4 3 1 2 6 5`, `S = 5 6 4 2 3 4 1 5`: `{1,2}` not
    a common interval (interval of T, not S); `{1,2,3,4}` is (a contiguous window in both).

## Deviations and assumptions

The artifact records **one API-shape assumption, not a correctness gap**: the public method
keeps a `maxGap` parameter but the validated/tested behaviour is the **strict, gap-free**
common-interval model (Uno & Yagiura 2000; Heber & Stoye 2001) — the cluster's group set must
equal the set of all groups in some contiguous window in every genome, no foreign groups
inside. `maxGap` does **not** relax this in the validated path, and the **gene-teams** gapped
extension (Bergeron, Corteel & Raffinot 2002) is **not** implemented (its source was not
retrievable this session). No contradictions among sources — Bui-Xuan/Uno & Yagiura/Heber &
Stoye/Didier et al. agree on the interval definition, the "contiguous in every genome" cluster
rule, the size-≥2 constraint, and the sequence-with-duplicates generalisation.
