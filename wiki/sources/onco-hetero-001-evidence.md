---
type: source
title: "Evidence: ONCO-HETERO-001 (intratumor heterogeneity metrics — MATH + Shannon clonal diversity)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-HETERO-001-Evidence.md
sources:
  - docs/Evidence/ONCO-HETERO-001-Evidence.md
source_commit: 44c3fb9d26b8f3d5317135230f8320077e2d1322
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-HETERO-001

The validation-evidence artifact for test unit **ONCO-HETERO-001** — **Tumor Heterogeneity Analysis
(MATH score, Shannon clonal diversity, subclone count, subclonal fraction)**. The **seventeenth
ingested unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in
[[intratumor-heterogeneity-metrics]] — the scalar-summary metric layer over the clonal structure.
[[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (five, mutually consistent — no contradictions):**
  - **Mroz & Rocco (2013)** *Oral Oncology* 49(3):211 (rank 1) — the **MATH** definition, verbatim
    **MATH = 100 · MAD/median**: the ratio of the width (median absolute deviation) to the centre
    (median) of the tumour's mutant-allele-fraction distribution; higher = more ITH.
  - **Mroz et al. (2015)** *PLOS Medicine* 12(2):e1001786 (rank 1) — the **1.4826** MAD scaling
    (makes the expected MAD of a normal variable equal its SD) and the full procedure confirming
    **MATH = 100 · 1.4826 · median(|f − median f|) / median(f)**.
  - **maftools `mathScore.R`** (rank 3) — reference code `abs.med.dev = abs(vaf − median(vaf))`;
    `pat.mad = median(abs.med.dev)·100`; `pat.math = pat.mad·1.4826/median(vaf)`, algebraically
    identical; rescales VAFs by /100 only if `max > 1` (this unit requires VAFs already in [0,1]).
  - **Liu & Zhang (2017)** *BMC Genomics* 18:457, PMC5468233 (rank 1) — **Shannon** ITH
    **S = −Σ pᵢ ln(pᵢ)** with the **natural log**, over clonal frequencies; **richness = number of
    clones** (monoclonal = 1, polyclonal ≥ 2).
  - **Shannon (1948)** (rank 4, corroborating) — H' = −Σ pᵢ ln pᵢ, base *e* usual in ecology; H rises
    with richness and evenness, → 0 as one class dominates.
  - **Landau et al. (2013)** *Cell* 152(4):714 — reused from ONCO-CLONAL-001: subclonal ⟺ CCF < 0.95
    (`ClonalCcfThreshold`), used for the **subclonal fraction**.

- **Documented corner cases / failure modes:** zero median VAF → MATH undefined (division by zero,
  must throw); all-identical VAFs or a single mutation → MAD = 0 ⇒ **MATH = 0**; single clone
  (richness 1, p = 1) → **H = 0**; k equal clones → **H = ln k** (maximum for that richness).

- **Datasets (hand-derived, exact):**
  - **MATH odd** — VAFs {0.10,0.20,0.30,0.40,0.50}, median 0.30, raw MAD 0.10, ×1.4826 = 0.14826 →
    **MATH 49.42**.
  - **MATH even** — VAFs {0.20,0.40,0.60,0.80}, median 0.50, raw MAD 0.20, ×1.4826 = 0.29652 →
    **MATH 59.304**.
  - **Shannon** — {1.0} → **0.0**; {0.5,0.5} → −ln 0.5 = **0.6931471805599453**;
    {0.25,0.25,0.25,0.25} → ln 4 = **1.3862943611198906**.

- **Coverage recommendations (7 items):** MUST — MATH 49.42 / 59.304 on the odd/even oracles; MATH = 0
  when all VAFs identical; Shannon H = ln k for k equal clones and H = 0 for a single clone; subclone
  count = number of occupied CCF clusters (Liu richness); subclonal fraction = #(CCF < 0.95)/n
  (Landau). SHOULD — zero-median VAF throws, null/empty/out-of-range throw. COULD — invariant
  MATH ≥ 0 over varied inputs (registry `ITH_score ≥ 0`).

## Deviations and assumptions

- **ASSUMPTION — Shannon clone fractions pᵢ = per-CCF-cluster mutation proportions.** Liu 2017
  computes Shannon over "clonal frequencies"; this unit takes the pᵢ as the proportion of mutations
  assigned to each CCF cluster (from the ONCO-CCF-001 clustering). Standard operationalisation; using
  cluster sizes (counts) is the natural per-mutation reading — a modelling choice, not
  source-mandated.
- **ASSUMPTION — median for even counts.** R's `median` (maftools) averages the two central order
  statistics; replicated. Sources do not enumerate an even-count rule, so the standard R behaviour is
  adopted (drives the 59.304 even-count oracle).

No source contradictions — the three MATH sources are algebraically identical, Liu 2017 and Shannon
1948 agree on the natural-log Shannon index, and both assumptions are operationalisations the sources
leave open rather than statements the sources contradict.
