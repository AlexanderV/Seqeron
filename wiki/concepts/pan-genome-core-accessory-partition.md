---
type: concept
title: "Pan-genome partition (core/accessory/unique + genomic fluidity + open/closed)"
tags: [comparative-genomics, pan-genome, algorithm]
sources:
  - docs/Evidence/PANGEN-CORE-001-Evidence.md
  - docs/Validation/reports/PANGEN-CORE-001.md
source_commit: 9957b475387e91865d553279d914fddef30ae85b
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pangen-core-001-evidence
      evidence: "Test Unit ID: PANGEN-CORE-001 ... Algorithm: Core / Accessory / Unique genome construction (pan-genome partitioning), genome fluidity, open/closed classification"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:pan-genome-gene-clustering
      source: pangen-core-001-evidence
      evidence: "ASSUMPTION: this repository's ConstructPanGenome delegates to the in-repo ClusterGenes ... The partitioning logic (core/accessory/unique by cluster occupancy, fluidity, openness) operates on the gene clusters that ClusterGenes (PANGEN-CLUSTER-001) produces"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:genome-comparison-core-dispensable
      source: pangen-core-001-evidence
      evidence: "Both realise the Tettelin 2005 core/dispensable pan-genome model; this unit is the N-genome occupancy-based partition (+ fluidity + Heaps openness), COMPGEN-COMPARE-001 is the pairwise two-genome RBH-based CompareGenomes"
      confidence: medium
      status: current
---

# Pan-genome partition (core/accessory/unique + genomic fluidity + open/closed)

`ConstructPanGenome` is the **pan-genome family**'s (`PANGEN-*`) partitioning unit: given gene
**clusters** across **N genomes**, it splits them into **core / accessory / unique** by cluster
**occupancy**, computes the pairwise **genomic fluidity** of the collection, and classifies the
pan-genome as **open or closed** by the **Heaps'-law** decay exponent. It sits **downstream** of the
gene-family clustering step [[pan-genome-gene-clustering]] (PANGEN-CLUSTER-001, which produces the
clusters this unit counts) and is the **N-genome, occupancy-based** counterpart of the pairwise
reciprocal-best-hit pipeline [[genome-comparison-core-dispensable]] (COMPGEN-COMPARE-001) — both
realise the Tettelin et al. 2005 core/dispensable model, one over many genomes by occupancy, the
other over two genomes by RBH. The **single-copy subset** of the core it defines is what
[[phylogenetic-marker-selection]] (PANGEN-MARKER-001) picks as phylogenetic markers. Validated under
test unit **PANGEN-CORE-001**; the validation record
is [[pangen-core-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The core / accessory / unique partition (occupancy)

Tettelin et al. 2005 (PNAS, the paper that coined "pan-genome") define the **pan-genome** as the
union of all gene families across a set of genomes, the **core genome** as families present in
*every* genome (the intersection), and the **accessory / dispensable genome** as families present in
some but not all — with **strain-specific / unique** families found in a single genome. Each cluster
is assigned by its **occupancy** = the number of distinct genomes it appears in:

- **Core** — occupancy meets the core threshold (present in ~all genomes).
- **Accessory / shell** — present in more than one but fewer than the core threshold.
- **Unique / cloud** — occupancy 1 (a single genome).

The **core threshold is fractional** (Roary, Page et al. 2015): a family is core when it is in **"at
least 99% of samples"**, i.e. `occupancy / N ≥ coreFraction` — **not** `floor(coreFraction·N)`. The
99% default tolerates assembly error in large datasets. Because the rule is fractional it is strict
for small N: with N = 3 and coreFraction 0.99, only 3/3 is core; **2/3 (66.7%) is accessory**, not
core. Roary's four-tier refinement is **hard core >99% / soft core 95–99% / shell 15–95% / cloud
<15%**.

This fractional rule was itself the outcome of validation: the two-stage verdict is recorded in
[[pangen-core-001-report]] — **Stage A/B both PASS-WITH-NOTES, State CLEAN**. The original code and
description used the unsourced floor `floor(coreFraction·N)`, which wrongly made 2/3 core at N=3;
the defect was fixed to the fractional `occupancy/N ≥ coreFraction` (shared `IsCoreOccupancy`) across
code, tests, and all three description artifacts. Fluidity and Heaps openness were already correct
(the single-order α fit is a disclosed simplification). Full suite 6509/0.

Worked partition (occupancy oracle): 3 genomes; c1 in all 3, c2 in 2, c3/c4/c5 each in 1. With
`coreFraction = 1.0` (coreThreshold 3) → **core {c1}**, **accessory {c2}**, **unique {c3,c4,c5}**.

## Genomic fluidity (Kislyuk 2011)

**Genomic fluidity** `φ` measures gene-content dissimilarity averaged over genome *pairs*:

```
φ = [2 / (N(N−1))] · Σ_{k<l} (U_k + U_l) / (M_k + M_l)
```

where `N` = number of genomes, `U_k` = gene families **unique** to genome k, `M_k` = **total**
families in genome k. `φ ∈ [0, 1]`: **0** = identical gene content, **1** = complete dissimilarity;
a fluidity of 0.1 means a pair shares 90% of genes on average. The jackknife variance is
`σ² = ((N−1)/N)·Σ_i(φ̂_(i) − φ̂)²`, where `φ̂_(i)` drops genome i.

Worked fluidity oracle (3 genomes): A = {c1,c2,c3}, B = {c1,c2,c4}, C = {c1,c5,c6}. Pair terms are
2/6, 4/6, 4/6, so `φ = (1/3)(2/6 + 4/6 + 4/6) = 10/18 = `**`0.5̄`**. Bounds: identical content → 0,
disjoint content → 1.

## Open vs closed (Heaps' law)

The pan-genome's **growth regime** is classified by fitting **Heaps' law** `N = k·n^(-alpha)` to the
number of *new* gene families observed as genomes are added in random order (micropan `heaps()`),
returning the **decay exponent alpha**:

- **Open** ⟺ **alpha < 1** (Wikipedia: alpha ≤ 1) — new genes keep accumulating, no asymptote
  (e.g. *E. coli*).
- **Closed** ⟺ **alpha > 1** — few new genes per added genome, the pan-genome size approaches a
  limit (e.g. *S. pneumoniae*).

micropan's verbatim rule: *"if alpha<1.0 the pan-genome is open, if alpha>1.0 it is closed."* The
decay-exponent fit itself — presence/absence binarization, the first-appearance new-gene curve, and
the bounded power-law least-squares — is the dedicated growth-model unit [[pan-genome-heaps-law-fit]]
(PANGEN-HEAP-001); this partition only *reports* its open/closed verdict as one output alongside the
occupancy partition and fluidity.

## Edge cases

- **Empty input** → empty partition.
- **Single genome (N < 2)** → no genome pairs, **fluidity undefined → taken as 0**; every occupancy-1
  cluster is unique.
- **< 3 genomes** → the Heaps decay exponent is degenerate (the new-gene curve needs several genomes),
  so open/closed is not meaningful.
- **Empty genome pair** (`M_k + M_l = 0`) → the fluidity term is undefined and contributes **0**.

## Assumptions (source-backed)

1. **Clustering identity metric (assumption).** Roary/Tettelin cluster gene families by BLASTP
   percentage identity; `ConstructPanGenome` delegates clustering to the in-repo `ClusterGenes` (the
   k-mer Jaccard heuristic of [[pan-genome-gene-clustering]], PANGEN-CLUSTER-001). The partition logic
   under test (occupancy core/accessory/unique, fluidity, openness) is **independent** of the upstream
   identity metric, so test inputs use identical or fully-disjoint sequences where occupancy is
   unambiguous.
2. **Empty-pair fluidity convention (assumption).** Pairs with `M_k + M_l = 0` contribute 0 (the
   neutral element for the otherwise-undefined term); only arises for empty genomes.

## Reference tools

Definitions trace to **Tettelin et al. 2005** (PNAS, pan-genome core/dispensable model) and
**Tettelin et al. 2008** (Curr Opin Microbiol, open/closed), **Kislyuk et al. 2011** (BMC Genomics
12:32, the fluidity formula and jackknife variance), **Page et al. 2015** (Roary, the fractional 99%
core rule and hard-core/soft-core/shell/cloud tiers), and **micropan** (`heaps()` alpha openness
criterion, `fluidity()` corroboration), with the **Wikipedia Pan-genome** article corroborating the
Heaps'-law classification. No source contradictions — the occupancy partition, fluidity formula, and
openness criterion are mutually consistent across the primaries and reference implementations.
