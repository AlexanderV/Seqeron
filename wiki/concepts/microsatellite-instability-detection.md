---
type: concept
title: "Microsatellite instability (MSI) detection — unstable-loci fraction + status"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-MSI-001-Evidence.md
source_commit: ea6d7a9858e7a4c0541dacc8d3df2a8b227021d9
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-msi-001-evidence
      evidence: "Test Unit ID: ONCO-MSI-001 ... Algorithm: Microsatellite Instability (MSI) Detection — fraction-of-unstable-loci scoring and status classification"
      confidence: high
      status: current
---

# Microsatellite instability (MSI) detection

The Oncology family's **microsatellite-instability** unit (**ONCO-MSI-001**): microsatellite
instability is the accumulation of insertion/deletion length changes at short tandem-repeat
(microsatellite) loci caused by defective DNA mismatch repair. MSI status is both a Lynch-syndrome
marker and an immunotherapy biomarker. This unit computes an MSI **score** and classifies MSI
**status** two ways — a continuous unstable-loci fraction and the categorical Bethesda marker-count
rule. The literature-traced record is [[onco-msi-001-evidence]]; [[test-unit-registry]] tracks the
unit and [[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## 1. Continuous MSI score (MSIsensor / MSIsensor2)

```
MSI score = (number of unstable/"msi" loci) / (number of valid evaluated loci)   # as a %
```

The score is the **fraction of unstable microsatellite loci among all valid evaluated loci**
(MSIsensor2 README, verbatim: "the msi score (number of msi sites / all valid sites) can be
calculated"). In MSIsensor (Niu et al. 2014) each site is tested by a **chi-square test** comparing
the tumor vs matched-normal repeat-length (allele-frequency) distributions, and a site is called
unstable/somatic when the difference is significant under a default **FDR 0.05**; the final score is
"the percentage of microsatellite sites with a somatic indel". MSIsensor2 keeps the same
fraction-of-unstable-loci definition in **tumor-only** mode (no matched normal).

**MSI-High cutoff:** `MSI score ≥ 20%` (boundary **inclusive** — 20% is MSI-High, 16% is not;
MSIsensor2 README, verbatim "msi high: msi score >= 20%"). Below 20% is classified only as
**not-High** — see the assumption below. MSIsensor's original 2014 cohort separated MSI from MSS
near a much lower **3.5%** score (70/71 MSI > 3.5, 165/168 MSS < 3.5); the 3.5% figure is a
dataset-specific separation point, not the recommended production cutoff (20%).

## 2. Categorical Bethesda classification (Boland et al. 1998)

Over the validated **5-marker Bethesda reference panel** (classic markers BAT-25, BAT-26, D2S123,
D5S346, D17S250), status is decided by the **count** of unstable markers:

| Unstable / 5 markers | Status |
|----------------------|--------|
| 0 / 5 | MSS |
| 1 / 5 | MSI-L (MSI-Low) |
| 2 / 5 | MSI-H |
| 3 / 5 | MSI-H |
| 5 / 5 | MSI-H |

- **MSI-H** — "two or more of the five markers show instability" (≥ 2 of 5).
- **MSI-L** — "only one of the five markers shows instability" (exactly 1 of 5).
- **MSS** — no marker shows instability (0 of 5).

The revised 2004 Bethesda criteria (British Journal of Cancer 2014;111:813, snippet only) restate
this in **fraction** form — MSI-H = ≥ 2/5 markers **or ≥ 40%**; MSI-L = 1/5 **or ≥ 20% and < 40%**;
MSS = no unstable markers — consistent with Boland 1998, used only as a cross-check.

## 3. Two inputs, two classifiers (key modelling choice)

MSIsensor2 defines **only a binary MSI-H cutoff (≥ 20%)** on the continuous score — it does **not**
define an MSI-L band on the continuous fraction. So the two classifiers apply to different inputs:

- **Continuous fraction score** → source-backed **binary** MSI-H (≥ 20%) vs not-High only.
- **Discrete marker count** (5-marker panel) → the three-way Bethesda MSS / MSI-L / MSI-H.

No MSI-L band is invented for the continuous score.

## Worked oracles

| Input | Result |
|-------|--------|
| 5 unstable / 25 valid | 20% → MSI-H (boundary, inclusive) |
| 4 unstable / 25 valid | 16% → not MSI-H |
| 0 unstable / 25 valid | 0% → MSS-range |
| 0 / 5 markers | MSS |
| 1 / 5 markers | MSI-L |
| 2 / 5 markers | MSI-H |

## Corner cases and assumptions

- **Score invariant:** `0 ≤ score ≤ 1`; the 20% boundary is inclusive (≥ 20% High).
- **Zero valid loci → undefined.** The score is `unstable / valid`; with no sufficiently covered
  locus there are zero valid sites and the score is undefined (division by zero → throws).
- **Invalid counts throw.** Counts must satisfy `0 ≤ unstable ≤ valid`; `unstable > valid` or
  negative counts are invalid input.
- **MSS vs MSI-L ambiguity (Boland 1998).** With only the 5-marker panel, distinguishing true MSS
  (0/5) from MSI-L (1/5) is unreliable; a larger panel is recommended. The 1-marker call is
  inherently low-confidence.
- **ASSUMPTION — no MSI-L band on the continuous score.** MSIsensor2 gives only the binary MSI-H
  cutoff; the categorical MSS/MSI-L/MSI-H scheme is therefore applied to the **discrete
  marker-count** input (Boland 1998), and the continuous fraction is classified only as MSI-H
  (≥ 20%) vs not-High. No MSI-L band is fabricated for the fraction score.

## Relation to the oncology family

MSI is an independent **genomic-biomarker** unit: it consumes per-locus stability calls (or an
unstable/valid count), not the allele-specific copy-number substrate. It is a distinct
immunotherapy/mismatch-repair biomarker, orthogonal to the copy-number-scar
[[homologous-recombination-deficiency-score]] (which sums LOH + TAI + LST from allele-specific
segments) and to the clinical-interpretation units [[cancer-variant-tier-classification-amp-asco-cap]]
and [[clinical-actionability-oncokb-levels]]. High MSI (MSI-H) is one of the tumor phenotypes those
interpretation layers act on.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for the MSI unstable-loci fraction and
status classification. The score definition + 20% cutoff are MSIsensor2 (Niu-lab README); the
per-site chi-square/FDR test and the percentage-of-somatic-sites definition are Niu et al. 2014; the
categorical MSS/MSI-L/MSI-H marker-count rule and the 5-marker Bethesda panel are Boland et al. 1998,
cross-checked against the 2004 revised-Bethesda fraction form. **Not for clinical or diagnostic use.**
No source contradictions — the continuous (MSIsensor) and categorical (Bethesda) definitions cover
disjoint inputs and agree where they overlap.
