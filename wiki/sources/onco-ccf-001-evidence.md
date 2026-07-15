---
type: source
title: "Evidence: ONCO-CCF-001 (cancer cell fraction estimation + 1D clonal/subclonal clustering)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-CCF-001-Evidence.md
sources:
  - docs/Evidence/ONCO-CCF-001-Evidence.md
source_commit: 280c1e8efb03103002e158f47f67a1e48267d3b0
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-CCF-001

The validation-evidence artifact for test unit **ONCO-CCF-001** — **Cancer Cell Fraction (CCF) point
estimation and 1D CCF clustering into clones/subclones**. The **fifth ingested unit of the Oncology
family** and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern. The distinct method is synthesized in
[[cancer-cell-fraction-clonal-clustering]] (the standalone estimator + the deterministic
clustering); the shared CCF closed form is also carried by the upstream
[[allele-specific-copy-number-ascat]] §4. [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent — the CCF formula corroborated three ways, no contradictions):**
  - **Tarabichi et al. (2021)** *Nature Methods* 18:144, subclonal-reconstruction guide, PMC7867630
    (rank 1) — Box 1 definitions: **CCF = CP / purity** (CP = cellular prevalence); the closed form
    **CCF = f·(ρ·N_T + 2(1−ρ)) / (ρ·m)** (VAF f, purity ρ, tumor total CN N_T, multiplicity m);
    **multiplicity m = f·(ρ·N_T + 2(1−ρ)) / ρ** rounded to nearest non-zero integer for clonal CN
    regions; **clonal-cluster rule** — "the cluster with the highest CP can be deemed clonal."
  - **Zheng et al. (2022) PICTograph** *Bioinformatics* 38(15):3677 (rank 1) — the VAF generative
    model **VAF = (m·CCF·p)/(c·p + 2(1−p))**; solving for CCF gives the identical closed form; no
    explicit [0,1] cap in the generative model (CCF treated as a continuous fraction).
  - **McGranahan et al. (2016)** *Science* 351:1463 (rank 1) — observed mutation copy number
    **n_mut = VAF·(1/p)·[p·CN_t + CN_n·(1−p)]** (CN_n = 2); **CCF = n_mut / multiplicity** — the
    same formula, plus the clonal definition (a mutation in all cancer cells has CCF = 1).
  - **CNAqc vignette (Caravagna lab)** (rank 3, reference implementation) — worked outputs
    **VAF 0.883 / m 2 → CCF 0.993** and **VAF 0.471 / m 1 → CCF 1.06**; demonstrates the raw formula
    can exceed 1 from sampling noise (does not intrinsically bound CCF ≤ 1).
  - **Lloyd 1982 / k-means (Wikipedia)** (rank 4, for the cited primary Lloyd 1982) — assignment
    (nearest-centroid, least squared Euclidean distance), update (recompute means), objective
    (minimise within-cluster sum of squares WCSS).

- **Documented corner cases / failure modes:** **CCF > 1 from noise** (CNAqc 1.06) — reported CCF
  capped at 1, raw exposed; **multi-copy loci (N_T > 2)** — multiplicity m ≠ 1 must be supplied,
  else CCF ambiguous (m ≥ 1, ≤ tumor CN); **purity ρ divides** → must be in (0,1].

- **Datasets (deterministic, hand-derived):**
  - **CNAqc worked outputs:** VAF 0.471 / m 1 → CCF ~1.06 (uncapped).
  - **Exact CCF cases A–F** from CCF = VAF·(ρ·N_T + 2(1−ρ))/(ρ·m): A clonal-diploid (0.40, 0.80, 2,
    1)→1.0; B subclonal (0.20, 0.80, 2, 1)→0.5; C multi-copy m=2 (0.50, 1.0, 4, 2)→1.0; D purity-0.5
    (0.25, 0.50, 2, 1)→1.0; E noise (0.471, 1.0, 2, 1)→0.942; **F overshoot→cap (0.60, 0.80, 2, 1)→
    raw 1.50, reported 1.0**.
  - **1D clustering (deterministic Lloyd k=2):** CCFs {1.0, 0.98, 0.96, 0.50, 0.48, 0.52} → centroids
    {0.50, 0.98}, low = indices {3,4,5} / high = {0,1,2}, **clonal = high (centroid 0.98)**.

- **Coverage recommendations (8 items):** MUST — `EstimateCCF` returns exact A–D CCFs within 1e-10;
  caps reported CCF at 1.0 when raw > 1 (case F) and exposes raw 1.5; enforces m ∈ [1, tumor CN];
  rejects ρ ∉ (0,1], VAF ∉ [0,1], CN < 1; `ClusterCCFValues` returns the exact centroids/assignments
  and identifies the highest-centroid cluster as clonal; and is **deterministic** across repeated /
  shuffled runs. SHOULD — k=1 → single mean cluster; null/empty handling on both methods.

## Deviations and assumptions

- **ASSUMPTION — reported CCF capped to [0,1].** Raw can exceed 1 (CNAqc 1.06); the unit reports
  min(1, raw) as bounded CCF while exposing the uncapped raw value. Justification: registry invariant
  + McGranahan clonal definition (mutation in all cancer cells → CCF = 1); no source forbids capping.
- **ASSUMPTION — 1D clustering = deterministic Lloyd k-means with quantile seeding.** Sources name
  CCF clustering broadly (Dirichlet process, variational beta mixtures) but the unit requires a
  *deterministic, well-defined* 1D method (no fabricated Dirichlet process). Lloyd 1982 is fully
  specified; determinism via sorting + evenly-spaced-quantile centroid seeding (no RNG). The
  clonal-cluster rule (highest centroid) is source-backed (Tarabichi 2021).

No source contradictions — the CCF formula is corroborated three independent ways (Tarabichi / Zheng
/ McGranahan) and the two assumptions are source-consistent, not source-contradicting.
