---
type: source
title: "Evidence: ONCO-TMB-001 (Tumor mutational burden — mutations/Mb + TMB-high)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-TMB-001-Evidence.md
sources:
  - docs/Evidence/ONCO-TMB-001-Evidence.md
source_commit: 701e1721ea2175070a8479c7d18e4e6f38be7076
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-TMB-001

The validation-evidence artifact for test unit **ONCO-TMB-001** — **Tumor Mutational Burden (TMB)**:
the mutations-per-megabase value and the TMB-high classification. The **thirty-fifth ingested unit of
the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its
own concept, [[tumor-mutational-burden]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Chalmers et al. (2017), Genome Medicine 9:34** (rank 1, PMC5395719) — TMB = "number of somatic,
    coding, base substitution, and indel mutations per megabase of genome examined"; FoundationOne
    **315-gene panel = 1.1 Mb** coding denominator; **synonymous alterations are counted** (to reduce
    sampling noise) before filtering; known COSMIC/tumor-suppressor/dbSNP/ExAC germline/driver
    alterations are **removed before counting**; below **~0.5 Mb** sequenced, panel-vs-WES deviation
    rises sharply (estimate unstable).
  - **Marcus et al. (2021), FDA Approval Summary, Clin Cancer Res 27(17):4685** (rank 2, PMC8416776) —
    **TMB-High ≡ TMB ≥ 10 mut/Mb** (inclusive); pembrolizumab tumor-agnostic approval **June 16 2020**;
    companion diagnostic **FoundationOne CDx (F1CDx)**; WES TMB = nonsynonymous SNVs + indels in coding
    regions.
  - **Fancello / Sholl review + Friends of Cancer Research Harmonization Project (Merino et al. 2020,
    J Immunother Cancer 8:e000147)** (rank 1–2, PMC7710563) — harmonized definition = somatic mutations
    per Mb of interrogated sequence; **report in mut/Mb** for comparability; confirms the FDA ≥10 mut/Mb
    cutoff; current F1CDx = **~0.8 Mb over 324 genes**.

- **Documented corner cases / failure modes:** `targetRegionMb = 0` → division by zero (throws);
  TMB-high threshold is inclusive (`≥ 10`, exactly 10.0 is high); small panels (< 0.5 Mb) still compute
  the ratio but are known-unstable (documentation, not an error); negative count/region → invalid input;
  filtering (germline/driver removal) happens upstream — this library counts the caller-supplied somatic
  list per ONCO-SOMATIC-001.

- **Datasets (deterministic worked oracles):**
  - **315-gene panel:** 11 mut / 1.1 Mb = 10.0 → TMB-High.
  - **WES:** 300 mut / 30 Mb = 10.0 → TMB-High.
  - **FDA boundary cases:** 99/10 = 9.9 → not high · 100/10 = 10.0 → high (inclusive) · 150/10 = 15.0 →
    high · 0/10 = 0.0 → not high.

- **Coverage recommendations:** MUST test `CalculateTMB` = count/regionMb (11/1.1=10, 300/30=10,
  150/10=15); MUST test `ClassifyTMB` boundary (9.9→Low, 10.0→High, 15.0→High, 0→Low); MUST test
  `targetRegionMb = 0` throws; SHOULD test small panel (< 0.5 Mb) computes without exception but is
  flagged unstable; SHOULD test negative count/region rejected; COULD test monotonicity invariant.

## Deviations and assumptions

- **ASSUMPTION — two-tier (High vs Not-High) using only the FDA ≥10 cutoff.** The only harmonized,
  authoritative threshold retrieved is the FDA/F1CDx **≥10 mut/Mb** TMB-high cutoff. The three-tier
  "Low <6 / Intermediate 6–20 / High >20" scheme some registries list has **no authoritative source**
  for the 6 and 20 boundaries (they are tumor-type-specific research cut-points). Per the
  evidence-first / no-fabrication policy, `ClassifyTMB` implements only TMB-High = TMB ≥ 10 vs
  TMB-Low = TMB < 10; the unsupported 6/20 constants are not implemented.

No source contradictions among the three references — Chalmers (formula + panel denominators + counting
rules), the FDA summary (≥10 cutoff + companion diagnostic), and the FoCR harmonization review (mut/Mb
reporting unit + cutoff confirmation) agree; the only flagged conflict is the unsupported registry 6/20
scheme, resolved in favour of the single source-backed ≥10 cutoff.
