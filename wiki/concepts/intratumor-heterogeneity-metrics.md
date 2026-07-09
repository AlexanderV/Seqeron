---
type: concept
title: "Intratumor heterogeneity metrics (MATH score + Shannon clonal diversity)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-HETERO-001-Evidence.md
source_commit: 44c3fb9d26b8f3d5317135230f8320077e2d1322
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-hetero-001-evidence
      evidence: "Test Unit ID: ONCO-HETERO-001 ... Algorithm: Tumor Heterogeneity Analysis (MATH score, Shannon clonal diversity, subclone count, subclonal fraction)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:cancer-cell-fraction-clonal-clustering
      source: onco-hetero-001-evidence
      evidence: "Subclone count = number of occupied CCF clusters (Liu 2017 richness) and the Shannon clone fractions pᵢ = proportion of mutations per CCF cluster — both consume the ONCO-CCF-001 CCF clustering output."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:clonal-subclonal-classification-ccf-posterior
      source: onco-hetero-001-evidence
      evidence: "Subclonal fraction = #(CCF < 0.95)/n reuses the Landau 2013 clonal CCF threshold 0.95 (ClonalCcfThreshold), the same 0.95 cut that ONCO-CLONAL-001 uses in its posterior classification rule."
      confidence: high
      status: current
---

# Intratumor heterogeneity metrics (MATH score + Shannon clonal diversity)

The **scalar-summary metric layer** of the Oncology family: reduce a tumor's clonal structure to a
handful of **intratumor heterogeneity (ITH)** numbers — the **MATH score**, the **Shannon clonal
diversity**, the **subclone count** (richness), and the **subclonal fraction**. Distinct from the
clonal-structure *reconstruction* units: it does not estimate per-mutation CCF or classify each
mutation — it **summarizes** the already-reconstructed structure into comparable scalars. Validated
under test unit **ONCO-HETERO-001**; the literature-traced record is [[onco-hetero-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern. Research-grade correctness reference — [[scientific-rigor|research-grade]],
**not for clinical or diagnostic use**.

## 1. MATH — Mutant-Allele Tumor Heterogeneity (Mroz & Rocco 2013)

MATH is the **ratio of the width to the centre of the mutant-allele-fraction (VAF) distribution**,
computed directly from NGS VAFs across mutated loci — **no clustering required**. Three mutually
corroborating sources give the identical formula (Mroz & Rocco 2013 *Oral Oncology*, Mroz 2015 *PLOS
Medicine*, maftools `mathScore.R`):

```
MATH = 100 · MAD / median(VAF)
     = 100 · 1.4826 · median(|VAF − median(VAF)|) / median(VAF)
```

The **1.4826** factor rescales the raw median absolute deviation so that the expected MAD of a
normally distributed variable equals its standard deviation (Mroz 2015 Methods). Higher MATH ⇒ more
intratumour genetic heterogeneity. maftools' reference code is byte-for-byte this:
`abs.med.dev = abs(vaf − median(vaf))`; `pat.mad = median(abs.med.dev) · 100`;
`pat.math = pat.mad · 1.4826 / median(vaf)`. (maftools rescales VAFs by /100 only if `max > 1`; this
unit requires VAFs already in [0, 1], so no rescaling.)

**Worked oracles (hand-derived, exact):**

| VAFs | median | raw MAD | ×1.4826 | MATH |
|------|--------|---------|---------|------|
| 0.10, 0.20, 0.30, 0.40, 0.50 (odd) | 0.30 | 0.10 | 0.14826 | **49.42** |
| 0.20, 0.40, 0.60, 0.80 (even) | 0.50 | 0.20 | 0.29652 | **59.304** |

## 2. Shannon clonal diversity (Liu 2017 / Shannon 1948)

The **Shannon index over clone fractions**, using the **natural logarithm** (base *e*, the ecology
convention):

```
H = −Σᵢ pᵢ · ln(pᵢ)          # pᵢ = fraction of the i-th clone
```

H rises with both **richness** (number of clones) and **evenness**, and → 0 as one clone dominates
or a single clone remains. Liu et al. (2017) compute it over "clonal frequencies"; here the pᵢ are
the **proportion of mutations assigned to each CCF cluster** (see Assumptions).

**Worked oracles:** `{1.0}` → **0.0**; `{0.5, 0.5}` → −ln 0.5 = **0.6931471805599453**;
`{0.25, 0.25, 0.25, 0.25}` → ln 4 = **1.3862943611198906**. In general k equal clones give the
maximum **H = ln k** for that richness.

This is the **same Shannon index math** as the metagenomics [[alpha-diversity]] unit (both natural
log, in nats) but a different domain — clone fractions here vs taxon abundances there.

## 3. Subclone count (richness) and subclonal fraction

- **Subclone count** = the **number of occupied CCF clusters** — Liu 2017's *richness* (monoclonal =
  1 clone; polyclonal = ≥ 2). This consumes the [[cancer-cell-fraction-clonal-clustering]]
  (ONCO-CCF-001) clustering directly.
- **Subclonal fraction** = **#(CCF < 0.95) / n** — the fraction of mutations whose cancer cell
  fraction falls below the **Landau 2013 clonal CCF threshold 0.95** (`ClonalCcfThreshold`, the same
  0.95 cut used by the posterior classifier [[clonal-subclonal-classification-ccf-posterior]]).

## Corner cases and failure modes

- **Zero median VAF** — MATH divides by median(VAF); a median of 0 makes MATH **undefined** (division
  by zero) and must be **rejected / throw**.
- **All identical VAFs / single mutation** — MAD = 0 ⇒ **MATH = 0**, a valid minimal value (no
  heterogeneity). A single mutation likewise gives median = that value, MAD = 0 → MATH = 0.
- **Single clone** — richness 1, p = 1, **H = −1·ln 1 = 0** (minimum diversity).
- **Invariant** — `ITH_score ≥ 0` (MATH ≥ 0 over all inputs); null / empty / out-of-range inputs
  throw.

## Assumptions (both source-consistent, neither source-mandated)

1. **Clone fractions for Shannon = per-cluster mutation proportions.** Liu 2017 computes Shannon over
   "clonal frequencies"; this unit takes pᵢ = proportion of mutations assigned to each CCF cluster
   (from ONCO-CCF-001 clustering). A standard operationalisation — using cluster **sizes** (counts) is
   the natural per-mutation diversity reading, but which pᵢ source (cluster CCF value vs cluster size)
   is a modelling choice.
2. **Median for even counts.** R's `median` (used by maftools) averages the two central order
   statistics; replicated here. Sources do not enumerate an even-count rule, so standard R behaviour
   is adopted (drives the 59.304 even-count oracle).

## Relationship to the clonal-structure units

This is the **metric / summary** layer sitting **downstream** of the reconstruction units. It does not
overlap them:

- [[cancer-cell-fraction-clonal-clustering]] (ONCO-CCF-001) **produces** the CCF vector and the
  clusters that MATH's subclone count and Shannon's pᵢ consume.
- [[clonal-subclonal-classification-ccf-posterior]] (ONCO-CLONAL-001) supplies the **0.95 clonal
  threshold** that the subclonal fraction reuses.
- MATH itself needs neither — it is a pure VAF-distribution statistic, computable straight from the
  per-locus mutant-allele fractions.

Both reconstruction units in turn consume the purity / copy-number substrate from the upstream
[[allele-specific-copy-number-ascat]].

## Scope and limitations

Five mutually consistent sources (Mroz & Rocco 2013; Mroz 2015; maftools `mathScore.R`; Liu 2017;
Shannon 1948) plus the reused Landau 2013 threshold — **no contradictions**. MATH and Shannon are
orthogonal ITH facets (VAF spread vs clone-fraction entropy). **Not for clinical or diagnostic use.**
