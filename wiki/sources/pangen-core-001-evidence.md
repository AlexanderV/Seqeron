---
type: source
title: "Evidence: PANGEN-CORE-001 (Pan-genome partition — core/accessory/unique + fluidity + open/closed)"
tags: [validation, comparative-genomics, pan-genome]
doc_path: docs/Evidence/PANGEN-CORE-001-Evidence.md
sources:
  - docs/Evidence/PANGEN-CORE-001-Evidence.md
source_commit: 9957b475387e91865d553279d914fddef30ae85b
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PANGEN-CORE-001

The validation-evidence artifact for test unit **PANGEN-CORE-001** — **pan-genome partitioning**:
`ConstructPanGenome` partitions gene clusters across N genomes into **core / accessory / unique**
by cluster **occupancy**, computes **genome fluidity** (Kislyuk 2011), and classifies the pan-genome
**open vs closed** by the **Heaps'-law** decay exponent. This is a **pan-genome family** (`PANGEN-*`)
Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm, its formulas, invariants,
worked oracles, and assumptions are summarized in [[pan-genome-core-accessory-partition]]. Its PANGEN
sibling is the upstream clustering step [[pan-genome-gene-clustering]] (PANGEN-CLUSTER-001, which
supplies the gene families this unit partitions); the pairwise two-genome relative is
[[genome-comparison-core-dispensable]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Tettelin et al. 2005** (PNAS, authority 1, the paper that coined "pan-genome") — **pan-genome**
    = union of all genes across the genomes; **core** = genes in *every* genome (intersection);
    **accessory / dispensable** = genes in some but not all; **strain-specific / unique** = genes in a
    single genome.
  - **Tettelin, Riley, Cattuto, Medini 2008** (Curr Opin Microbiol, authority 1) — core + dispensable
    decomposition; **open** pan-genome when the Heaps'-law power-law exponent indicates unbounded
    growth (**alpha < 1**), many bacteria are open (small conserved core + large dynamic accessory).
  - **Wikipedia "Pan-genome"** (authority 4, corroboration) — component vocabulary core / shell
    (accessory, ≥2 strains) / cloud (strain-specific); **Heaps' law** `N = k·n^(-alpha)`, **open ⟺
    alpha ≤ 1**, **closed ⟺ alpha > 1**; examples *E. coli* open, *S. pneumoniae* closed.
  - **Kislyuk, Haegeman, Bergman, Weitz 2011** (BMC Genomics 12:32, authority 1) — the **genomic
    fluidity** formula `φ = [2/(N(N−1))]·Σ_{k<l}(U_k+U_l)/(M_k+M_l)` (U = unique gene families per
    genome, M = total families per genome), range 0..1 (0 = identical gene content, 1 = fully
    dissimilar; 0.1 ⇒ pairs share 90% on average), plus the jackknife variance
    `σ² = ((N−1)/N)·Σ_i(φ̂_(i)−φ̂)²`.
  - **Page et al. 2015** (Roary, *Bioinformatics*, authority 3) — core = "a gene being in **at least
    99% of samples**" (default `-cd 99`), tolerating assembly error; clustering identity default 95%
    BLASTP; the four-tier vocabulary **hard core >99% / soft core 95–99% / shell 15–95% / cloud <15%**.
  - **micropan `heaps()` / `fluidity()`** (CRAN, authority 3, cites Tettelin 2008) — Heaps model fitted
    to the number of *new* gene clusters over random genome orderings, returns intercept + **decay
    parameter alpha**; verbatim openness criterion "**if alpha<1.0 the pan-genome is open, if alpha>1.0
    it is closed**"; corroborates the Kislyuk fluidity definition.
- **Corner cases / failure modes:** < 3 genomes → Heaps decay exponent degenerate (openness fit not
  meaningful); N < 2 → no genome pairs, fluidity undefined (taken as **0** by convention); a pair with
  `M_k + M_l = 0` (empty genomes) contributes **0**; the Roary 99% core rule is **fractional**
  (`occupancy / N ≥ coreFraction`), **not** `floor(coreFraction·N)` — with N=3, coreFraction 0.99 only
  3/3 is core, 2/3 (66.7%) is shell/accessory (corrected 2026-06-15 during validation).
- **Datasets (documented oracles):**
  - *Fluidity bounds* — two genomes 10% unique each → φ = 0.1; identical content → 0; disjoint → 1.
  - *Hand-derived 3-genome fluidity* — A={c1,c2,c3}, B={c1,c2,c4}, C={c1,c5,c6}: pair terms 2/6, 4/6,
    4/6 → φ = (1/3)(10/6) = **10/18 = 0.5̄**.
  - *Core/accessory/unique partition* — 3 genomes, c1 in all 3, c2 in 2, c3/c4/c5 each in 1;
    coreFraction 1.0 (coreThreshold 3) → core {c1}, unique {c3,c4,c5}, accessory {c2}.
- **Coverage recommendations:** MUST-test the fractional occupancy core rule, the closed-form fluidity
  on the 0.5̄ example, the fluidity 0/1 bounds, and open-vs-closed by alpha (open ⟺ alpha < 1);
  SHOULD-test empty→empty, single-genome (no pairs → fluidity 0, occupancy-1 clusters are unique);
  COULD-test the coreFraction = CoreGeneCount/TotalClusters invariant.

## Deviations and assumptions

Two assumptions, both source-backed, no deviations from the literature. (1) **Clustering identity
metric (assumption).** Roary/Tettelin cluster by BLASTP percentage identity; `ConstructPanGenome`
delegates clustering to the in-repo `ClusterGenes` (k-mer Jaccard heuristic, the separate
PANGEN-CLUSTER-001 unit / [[pan-genome-gene-clustering]]). The **partition logic under test**
(occupancy-based core/accessory/unique, fluidity, openness) is **independent** of the upstream
identity metric, so test inputs use identical or fully-disjoint sequences where occupancy is
unambiguous. (2) **Empty-pair convention (assumption).** Fluidity pairs with `M_k + M_l = 0`
contribute 0 (the neutral element for the undefined term); not stated by Kislyuk but only arises for
empty genomes. No source contradictions.
