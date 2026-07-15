---
type: source
title: "Evidence: ONCO-CLONAL-001 (clonal vs subclonal classification, CCF posterior)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-CLONAL-001-Evidence.md
sources:
  - docs/Evidence/ONCO-CLONAL-001-Evidence.md
source_commit: 730939482dd663a30a7dfb71f56e1f47e1bf1dd9
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-CLONAL-001

The validation-evidence artifact for test unit **ONCO-CLONAL-001** — **Clonal vs Subclonal Mutation
Classification (cancer cell fraction posterior)**. The **seventh ingested unit of the Oncology
family** and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern. The distinct method is synthesized in
[[clonal-subclonal-classification-ccf-posterior]] — the Bayesian CCF-posterior classifier, an
`alternative_to` the point-estimate + k-means [[cancer-cell-fraction-clonal-clustering]].
[[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (two primary papers, mutually consistent — no contradictions):**
  - **Landau et al. (2013)** *Cell* 152(4):714 (rank 1) — the ABSOLUTE-style CCF posterior: expected
    allele fraction **f(c) = αc/(2(1−α)+αq)** (α purity, c CCF, q per-locus absolute somatic copy
    number); posterior **P(c) ∝ Binom(a|N,f(c))** with a **uniform prior**, evaluated over a
    **100-point grid c∈[0.01,1]** and normalised; **classification rule (verbatim): clonal iff the
    posterior probability that CCF > 0.95 exceeds 0.5, subclonal otherwise**. Inputs: a mutation seen
    in `a` of `N` reads at a locus of copy number `q`, sample purity `α`. Subclonal-CCF sanity
    reference: median subclonal CCF 0.49.
  - **Satas, Zaccaria, El-Kebir, Raphael (2021) DeCiFering** *Cell Systems* 12(10):1004,
    PMC8542635 (rank 1) — the **multiplicity-general** CCF (Eq. 1) inverting **f(c) =
    ρ·M·c/(2(1−ρ)+ρ·N_tot)** to arbitrary SNV multiplicity M; clonal (CCF≈1) vs subclonal (CCF≪1)
    definition.

- **Documented corner cases / failure modes:** CCF lower bound > 0 (grid c∈[0.01,1] — a detected
  mutation sits in ≥1 cancer cell); **probabilistic, not point** classification (a point estimate
  near 1 with wide uncertainty from shallow coverage can still be subclonal); **multiplicity
  matters** — for M>1 the same VAF maps to a lower CCF, so ignoring M overestimates CCF.

- **Datasets (deterministic, grid-evaluated independently of the implementation):**
  - **Posterior classification cases** A1 (a=N=300, q=2, M=1, α=1.0 → mean 0.9995, P=1.0, Clonal),
    B2 (400/1000, α=0.8 → 0.9725, P=0.864, Clonal), C1 (240/1000, α=0.8 → 0.601, P=0, Subclonal),
    D (200/1000, α=1.0 → 0.401, P≈0, Subclonal), **E (100/100, M=2, α=1.0 → 0.994, P=0.998,
    Clonal)** — the multiplicity lift.
  - **Point-estimate `IdentifyClonalMutations`** (clonal iff CCF > 0.95, boundary excluded):
    {0.96, 0.95, 1.00, 0.50, 0.951} → clonal indices **{0, 2, 4}**.

- **Coverage recommendations (7 items):** MUST — clonal calls for deep-coverage variants with
  posterior mass > 0.5 above 0.95 (A1/B2/E); subclonal for CCF≪1 (C1/D); M=2 raises CCF for the same
  VAF (E); invariant ClonalCount+SubclonalCount=total & ClonalFraction=ClonalCount/total;
  `IdentifyClonalMutations` strict CCF>0.95 (0.95 excluded). SHOULD — null / purity∉(0,1] / invalid
  read counts / copy number / multiplicity / CCF∉[0,1] throw. COULD — empty variant set → empty
  calls, counts 0, ClonalFraction 0.

## Deviations and assumptions

- **ASSUMPTION — registry `ploidy` parameter is the per-variant local copy number `q`.** The registry
  stub was `ClassifyClonality(variants, purity, ploidy)`, but Landau's model uses the **per-locus**
  absolute copy number q, not a genome-wide ploidy scalar. The canonical method carries q per variant
  (`ClonalityVariant.LocalCopyNumber`) and takes `(variants, purity)`. Mirrors the prior ONCO-WGD
  decision where a registry scalar was superseded by per-segment data to match the authoritative
  definition. **Non-correctness-affecting** (API shape only): the numerical rule and outputs are
  exactly Landau's.

No source contradictions — Landau (posterior + 0.95/0.5 rule) and Satas (multiplicity generalisation)
are complementary, and the single assumption is API-shape only, not source-contradicting.
