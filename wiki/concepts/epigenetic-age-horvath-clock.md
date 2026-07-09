---
type: concept
title: "Epigenetic age (Horvath DNAm methylation clock)"
tags: [epigenetics, algorithm]
sources:
  - docs/Evidence/EPIGEN-AGE-001-Evidence.md
  - docs/algorithms/Epigenetics/Epigenetic_Age_Estimation.md
source_commit: e90a75989c52785a566e034477454379fd0535e0
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: epigen-age-001-evidence
      evidence: "Test Unit ID: EPIGEN-AGE-001 ... Algorithm: Epigenetic Age Estimation (Horvath DNA methylation clock) ... Algorithm Group: Epigenetics"
      confidence: high
      status: current
---

# Epigenetic age (Horvath DNAm methylation clock)

Estimating **DNA-methylation ("epigenetic") age** from methylation **β-values** measured at clock CpG
sites, following the **Horvath (2013) multi-tissue clock**. This is the **first ingested unit of the
Epigenetics family**; two siblings now exist — [[bisulfite-methylation-calling]] *produces* the per-CpG
β-values this clock consumes, and [[chromatin-state-prediction]] annotates chromatin state from
histone marks (CpG islands and DMRs remain queued). Validated under test
unit **EPIGEN-AGE-001**;
the validation record is [[epigen-age-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

DNA methylation at **CpG dinucleotides** changes systematically with age. Horvath (2013) selected
**353 CpGs** by **elastic-net regression** (α = 0.5, λ = 0.0226; 193 positively, 160 negatively
correlated with age) of a *transformed* version of chronological age, yielding a predictor that
generalises across tissues. Each site's methylation is a β-value in [0, 1] (fraction methylated).

## Two-stage model

DNAm age is computed in two stages (`EpigeneticsAnalyzer.CalculateEpigeneticAge`):

**1. Linear predictor** in *transformed-age* units — a weighted sum over the clock CpGs plus an
intercept:

```
Y = intercept + Σ_i (coef_i · β_i)
```

**2. Inverse calibration** `F⁻¹` (`anti.trafo`) mapping transformed age back to years, with the
calibration-break constant **`adult.age = 20`**:

```
F⁻¹(Y) = (1 + 20)·exp(Y) − 1   if Y < 0      (exponential branch)
F⁻¹(Y) = (1 + 20)·Y + 20        if Y ≥ 0      (linear branch)
```

The forward calibration `F` is **logarithmic below adult age and linear above** — reflecting the
paper's observation that the methylation–age dependence is "logarithmic until adulthood … linear
later in life". The two branches meet continuously at `(0, 20)`.

The embedded multi-tissue clock uses **intercept `0.695507258`** (`CoefficientTraining[1]` of
Additional file 3) and the **353-CpG `CoefficientTraining` weights**, retrieved verbatim from the
Springer supplement and cross-verified **byte-identical** against an independent GitHub mirror (all
353 `(CpGmarker, coefficient)` pairs). Exposed as `HorvathMultiTissueCoefficients` /
`HorvathMultiTissueIntercept`, the two-branch transform as `HorvathAntiTransform`.

## Worked oracles

**Transform in isolation** (`anti.trafo`, adult.age = 20):

| Linear predictor Y | Branch | DNAm age |
|--------------------|--------|----------|
| 0.0 | linear (Y ≥ 0) | 20.0 |
| 1.0 | linear | 41.0 |
| −1.0 | exponential | 21·e⁻¹ − 1 = 6.7254682646002895 |
| −2.5 | exponential | 21·e⁻²·⁵ − 1 = 0.7237849711018749 |

**Linear-predictor assembly** (intercept + Σ coef·β): intercept 0.695507258 with `(coef, β)` pairs
(0.0127, 0.5), (−0.0312, 0.8), (0.0245, 0.3) → `Y = 0.684247258` → linear branch `21·Y + 20` =
**34.369192418** years.

## Invariants

- **Y ≥ 0** → age = `21·Y + 20` (linear branch).
- **Y < 0** → age = `21·exp(Y) − 1` (exponential branch).
- **Y = 0** → age = **20.0** exactly (`Y < 0` false → linear branch → `21·0 + 20`).
- Age is **strictly increasing and continuous** in Y; the branches meet at `(0, 20)`.
- **CpGs absent from the coefficient table do not change the result** — only clock CpGs enter the
  weighted sum; an empty methylation map gives `F⁻¹(intercept)` (intercept-only, no CpG
  contributions).

## Contract and edge cases

- `methylationAtClockCpGs == null` → `ArgumentNullException`; `coefficients == null` →
  `ArgumentNullException`; **empty** `coefficients` → `ArgumentException` (an empty clock has no
  defined output). An empty methylation map is **valid**.
- CpG ids are matched by exact dictionary key (case sensitivity follows the supplied comparer).
- For extreme negative Y the exponential branch approaches −1 — a **mathematical artefact, not a
  biological age**; realistic predictors stay positive.

## Scope and limitations

The EPIGEN-AGE-001 Evidence unit covers only the **multi-tissue 353-CpG clock**. The algorithm doc
(`Epigenetic_Age_Estimation.md`, reviewed one day later) has since extended the implementation with
two more embedded clocks that share this concept's linear-predictor structure:

- **Horvath (2018) skin-&-blood clock** — 391 CpGs, intercept −0.447119319, **same** `anti.trafo`
  (adult.age = 20) path; `CalculateSkinBloodAge`.
- **Levine (2018) DNAm PhenoAge** — 513 CpGs, intercept 60.664, **no transform** (the linear
  predictor is already in years); `CalculatePhenoAge`.

The implementation does **not** perform the full Horvath pipeline's upstream **BMIQ normalisation** of
input β-values or array-platform handling — callers must pre-normalise. A
[[research-grade-limitations|research-grade]] simplification: a weighted-sum-plus-transform scorer over
a dictionary, not the end-to-end normalisation/QC pipeline. No occurrence search is performed, so the
repository suffix tree is not applicable.
