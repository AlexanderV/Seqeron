---
type: concept
title: "Tumor mutational burden (TMB) — mutations/Mb + TMB-high classification"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-TMB-001-Evidence.md
source_commit: 701e1721ea2175070a8479c7d18e4e6f38be7076
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-tmb-001-evidence
      evidence: "Test Unit ID: ONCO-TMB-001 ... Algorithm: Tumor Mutational Burden (TMB) — mutations/Mb and TMB-high classification"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:somatic-variant-calling-tumor-normal
      source: onco-tmb-001-evidence
      evidence: "this library counts the caller-supplied somatic-mutation list, leaving germline/driver filtering to the upstream somatic caller, per ONCO-SOMATIC-001"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:microsatellite-instability-detection
      source: onco-tmb-001-evidence
      evidence: "TMB and MSI are the two independent immunotherapy / genomic-biomarker units of the Oncology family; the FDA pembrolizumab TMB-high approval is a distinct tumor-agnostic biomarker alongside MSI-H."
      confidence: medium
      status: current
---

# Tumor mutational burden (TMB)

The Oncology family's **tumor-mutational-burden** unit (**ONCO-TMB-001**): TMB is the number of
somatic coding mutations per megabase of sequenced genome — a tumor-agnostic **immunotherapy
biomarker** (high TMB predicts response to immune-checkpoint inhibitors). This unit computes the
continuous TMB value and classifies **TMB-high** status against the FDA cutoff. The
literature-traced record is [[onco-tmb-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## 1. TMB value (mutations per megabase)

```
TMB = (counted somatic mutations) / (sequenced coding region in Mb)      # units: mut/Mb
```

TMB is "the number of somatic, coding, base substitution, and indel mutations per megabase of
genome examined" (Chalmers et al. 2017, Methods). The denominator is the **coding footprint of the
assay** in Mb — e.g. **1.1 Mb** for the historical FoundationOne 315-gene panel, **~0.8 Mb** for the
324-gene FoundationOne CDx, and **≈30–40 Mb** for whole-exome sequencing. The algorithm takes the
region size as a parameter so any assay can be scored.

**Mutation counting.** All base substitutions and indels in the targeted coding region are counted.
Panel methods (Chalmers 2017) additionally **include synonymous alterations** to reduce sampling
noise, whereas the WES analyses in the FDA/harmonization work counted **nonsynonymous** coding SNVs +
indels. Known driver/germline alterations (COSMIC hotspots, tumor-suppressor truncations, dbSNP/ExAC
germline, somatic-germline-zygosity predictions) are **filtered out before counting** — a raw variant
list overstates TMB. This unit counts the **caller-supplied somatic-mutation list**, leaving that
filtering to the upstream somatic caller ([[somatic-variant-calling-tumor-normal]]).

## 2. TMB-high classification

```
TMB-High  ⇔  TMB ≥ 10 mut/Mb        (boundary inclusive)
TMB-Low   ⇔  TMB < 10 mut/Mb
```

The **≥10 mut/Mb** cutoff is the FDA-approved TMB-high definition from the tumor-agnostic
pembrolizumab approval (June 16 2020), measured by the **FoundationOne CDx** companion diagnostic
(Marcus et al. 2021). The threshold is **inclusive**: exactly 10.0 mut/Mb is TMB-high; 9.9 is not.
The Friends of Cancer Research TMB Harmonization Project (Merino et al. 2020) confirms both the
mut/Mb reporting unit and the ≥10 cutoff.

## 3. Two-tier only — the single-cutoff modelling choice

**ASSUMPTION — Two-tier (High vs Not-High) using only the FDA ≥10 cutoff.** The only authoritative,
harmonized TMB threshold retrieved is the FDA/F1CDx **≥10 mut/Mb** boundary. Some registries list a
three-tier "Low <6 / Intermediate 6–20 / High >20" scheme, but **no authoritative source defines
those 6 and 20 boundaries** — they are tumor-type-specific research cut-points, not a harmonized
standard. Per the evidence-first / no-fabrication policy, `ClassifyTMB` implements only the
source-backed TMB-High = TMB ≥ 10 vs TMB-Low = TMB < 10; the unsupported 6/20 constants are **not**
implemented (they would be fabricated).

## Worked oracles

| Counted mutations | Region (Mb) | TMB (mut/Mb) | Classification |
|-------------------|-------------|--------------|----------------|
| 11 | 1.1 | 10.0 | TMB-High (315-gene panel example) |
| 300 | 30 | 10.0 | TMB-High (WES example) |
| 150 | 10 | 15.0 | TMB-High |
| 100 | 10 | 10.0 | TMB-High (boundary, inclusive) |
| 99 | 10 | 9.9 | Not high (< 10) |
| 0 | 10 | 0.0 | Not high |

## Corner cases and assumptions

- **Region = 0 → throws.** Mb is in the denominator; `CalculateTMB` with `targetRegionMb = 0` is a
  division by zero and rejected.
- **Non-negative inputs.** Mutation count and region size are non-negative by definition; negative
  values are invalid input.
- **Small-denominator instability (documentation, not an error).** Below **~0.5 Mb** of sequenced
  coding region the percentage deviation of panel TMB from whole-exome TMB rises sharply — the panel
  estimate is unreliable (Chalmers 2017). The ratio is still mathematically defined, so the value is
  computed (no exception) but is flagged as known-unstable.
- **Monotonicity invariant.** For fixed region, TMB is non-decreasing in mutation count; for fixed
  count, non-increasing in region size (division property).

## Relation to the oncology family

TMB is an independent **genomic-biomarker** unit: it consumes a somatic-mutation **count** plus the
assay's coding footprint in Mb, not the allele-specific copy-number substrate. It sits alongside the
other tumor-agnostic immunotherapy biomarker, [[microsatellite-instability-detection]] (MSI-H), and
downstream of the somatic caller [[somatic-variant-calling-tumor-normal]] that supplies the filtered
mutation list. High TMB is one of the tumor phenotypes the clinical-interpretation layers
([[cancer-variant-tier-classification-amp-asco-cap]], [[clinical-actionability-oncokb-levels]]) act on,
and it correlates with the neoantigen load handled by [[neoantigen-peptide-generation]].

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for the TMB ratio and TMB-high
classification. The mut/Mb formula and synonymous-counting / pre-count filtering rules are Chalmers
et al. 2017; the ≥10 mut/Mb TMB-high cutoff and F1CDx companion diagnostic are the FDA pembrolizumab
approval (Marcus et al. 2021); the mut/Mb reporting unit and coding-mutation basis are cross-checked
against the Friends of Cancer Research Harmonization Project (Merino et al. 2020). **Not for clinical
or diagnostic use.** No source contradictions — the sole flagged conflict is the unsupported 6/20
three-tier scheme, resolved in favour of the single source-backed ≥10 cutoff (see assumption above).
