---
type: source
title: "Evidence: ONCO-MSI-001 (Microsatellite instability — unstable-loci fraction + status)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-MSI-001-Evidence.md
sources:
  - docs/Evidence/ONCO-MSI-001-Evidence.md
source_commit: ea6d7a9858e7a4c0541dacc8d3df2a8b227021d9
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-MSI-001

The validation-evidence artifact for test unit **ONCO-MSI-001** — **Microsatellite Instability (MSI)
detection**: the fraction-of-unstable-loci MSI score and MSI status classification. The **twenty-fourth
ingested unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its
own concept, [[microsatellite-instability-detection]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **MSIsensor2 README (niu-lab)** (rank 3, reference implementation) — MSI score = `number of msi
    sites / all valid sites` (as %); recommended **MSI-High cutoff ≥ 20%** (inclusive); tumor-only mode
    keeps the same fraction definition.
  - **Niu et al. (2014), Bioinformatics 30(7):1015–1016** (rank 1) — per-site **chi-square test** of
    tumor vs normal repeat-length distributions, unstable/somatic at default **FDR 0.05**; MSI score =
    "percentage of microsatellite sites with a somatic indel"; original cohort separated MSI vs MSS near
    a **3.5%** score (70/71 MSI > 3.5, 165/168 MSS < 3.5 — dataset-specific, not the production cutoff).
  - **Boland et al. (1998), Cancer Res 58(22):5248–5257** (rank 1, NCI consensus workshop) — over the
    validated **5-marker Bethesda panel** (BAT-25, BAT-26, D2S123, D5S346, D17S250): **MSI-H** = ≥ 2/5
    unstable, **MSI-L** = exactly 1/5, **MSS** = 0/5; MSS vs MSI-L unreliable with only 5 markers.
  - **British Journal of Cancer 2014;111:813** (rank 1, snippet only — full text behind login redirect)
    — revised 2004 Bethesda **fraction** form: MSI-H = ≥ 2/5 or ≥ 40%; MSI-L = 1/5 or ≥ 20% and < 40%;
    MSS = none. Consistent with Boland 1998; used only as a cross-check, not sole authority.

- **Documented corner cases / failure modes:** zero valid loci → score undefined (division by zero,
  throws); `unstable > valid` or negative counts → invalid input (throws); score invariant
  `0 ≤ score ≤ 1` with 20% boundary inclusive; MSS-vs-MSI-L 1-marker call inherently low-confidence
  (larger panel recommended); tumor-only mode leaves the fraction definition unchanged.

- **Datasets (deterministic worked oracles):**
  - **Continuous (MSIsensor2):** 5/25 → 20% → MSI-H (boundary) · 4/25 → 16% → not MSI-H · 0/25 → 0% →
    MSS-range.
  - **Categorical (Bethesda 5-marker, Boland 1998):** 0/5 → MSS · 1/5 → MSI-L · 2/5, 3/5, 5/5 → MSI-H.

- **Coverage recommendations:** MUST test the continuous score `unstable/valid` and MSI-H at the 20%
  boundary (≥ 20% High, < 20% not High); MUST test the Bethesda marker-count classification
  (0→MSS, 1→MSI-L, ≥2→MSI-H); MUST test the score invariant `0 ≤ score ≤ 1` + boundary exactness;
  SHOULD test zero-valid-loci → throws and `unstable > valid` / negative → throws; COULD test the 40%
  Bethesda fraction on larger panels (BJC 2014 cross-check).

## Deviations and assumptions

- **ASSUMPTION — MSI-L band applies only to the discrete marker count, not the continuous score.**
  MSIsensor2 defines only a binary MSI-H cutoff (≥ 20%) on the continuous fraction; it does not define
  an MSI-L band there. The three-way categorical MSS/MSI-L/MSI-H scheme is therefore applied to the
  **discrete marker-count** input (Boland 1998), and the continuous fraction is classified only as the
  source-backed binary MSI-H (≥ 20%) vs not-High. No MSI-L band is invented for the continuous score.

No source contradictions — the continuous MSIsensor fraction definition and the categorical Bethesda
marker-count rule cover disjoint inputs (a percentage of loci vs a count over a 5-marker panel) and
agree where they overlap (an all-stable sample is MSS-range / MSS).
