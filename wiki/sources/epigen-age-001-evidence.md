---
type: source
title: "Evidence: EPIGEN-AGE-001 (Epigenetic age, Horvath DNAm clock)"
tags: [validation, epigenetics]
doc_path: docs/Evidence/EPIGEN-AGE-001-Evidence.md
sources:
  - docs/Evidence/EPIGEN-AGE-001-Evidence.md
source_commit: e90a75989c52785a566e034477454379fd0535e0
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: EPIGEN-AGE-001

The validation-evidence artifact for test unit **EPIGEN-AGE-001** — **epigenetic ("DNAm") age
estimation** by the **Horvath (2013) multi-tissue DNA-methylation clock**. This is the **first
ingested unit of the Epigenetics family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The algorithm itself is synthesized in
its own concept, [[epigenetic-age-horvath-clock]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources:**
  - **Horvath S. (2013)** "DNA methylation age of human tissues and cell types", *Genome Biology*
    14:R115 (authority rank 1, primary paper) — the **353-CpG** multi-tissue predictor selected by
    **elastic-net regression** (α = 0.5, λ = 0.0226; 193 positively / 160 negatively age-correlated),
    a linear combination of CpG β-values **in transformed-age units**, and the calibration shape
    "logarithmic dependence until adulthood that slows to a linear dependence later in life".
  - **aldringsvitenskap/epigeneticclock** reference R (`horvath2013.R`, rank 3) — the verbatim forward
    `trafo` and inverse `anti.trafo` transforms with **`adult.age = 20`**; the inverse is the
    two-branch `DNAm age = (1+20)·exp(x) − 1` when `x < 0`, else `(1+20)·x + 20`.
  - **aldringsvitenskap/epigeneticclock** reference R (`StepwiseAnalysis.R`, rank 3) — the verbatim
    `predictedAge = anti.trafo(CoefficientTraining[1] + meth %*% CoefficientTraining[-1])`: the
    **intercept** is `CoefficientTraining[1]`, per-CpG weights are `CoefficientTraining[-1]`, and the
    linear predictor `intercept + Σ coef·β` is passed to `anti.trafo`.
  - **Horvath (2013) Additional file 3** (Springer supplement, rank 1) — the **353-CpG coefficient
    table** (`CoefficientTraining` column relating CpGs to *transformed* age); **intercept =
    0.695507258**; example weights cg00075967 = 0.12933661, cg00864867 = 1.59976405, cg09809672
    (EDARADD) = −0.391318905, cg27544190 = −0.869124446. **Memory-trap note:** cg22454769 is **not**
    in the 353-CpG table (confirmed absent).
  - **aldringsvitenskap/epigeneticclock `AdditionalFile3.csv`** (GitHub mirror, rank 3) — a
    **field-by-field diff of all 353 `(CpGmarker, CoefficientTraining)` pairs against the Springer
    supplement showed zero differences (byte-identical)**. This cross-check is what authorised
    embedding the table.
  - **perishky/meffonym** (rank 3, second independent implementation) — confirms the Horvath DNAmAge
    model applies the `anti.trafo` inverse log-transform, corroborating the two-branch transform.

- **Documented corner cases / failure modes:** the two-branch boundary at `x = 0` takes the **linear**
  branch → exactly **20** (adult age); negative linear predictors (young / hypomethylated) take the
  **exponential** branch, approaching −1 as `x → −∞` (mathematical limit, not biological); CpGs
  **without a coefficient are ignored** (only clock CpGs enter the matrix product).

- **Datasets (documented oracles):**
  - *`anti.trafo` transform table* — `Y = 0.0` → **20.0** (linear); `Y = 1.0` → **41.0** (linear);
    `Y = −1.0` → `21·e⁻¹ − 1` = **6.7254682646002895** (exp); `Y = −2.5` → `21·e⁻²·⁵ − 1` =
    **0.7237849711018749** (exp).
  - *Linear-predictor assembly* — intercept 0.695507258 with `(coef, β)` pairs
    (0.0127, 0.5), (−0.0312, 0.8), (0.0245, 0.3) → `Y = 0.684247258`; linear branch
    `21·Y + 20` = **34.369192418**.

- **Test-coverage recommendations:** MUST — linear branch `21·Y + 20` (Y = 0.684247258 → 34.369192418);
  boundary `Y = 0` → 20.0; negative branch `21·exp(Y) − 1` (Y = −1 → 6.7254682646002895); CpGs absent
  from the table ignored; null/empty coefficient table and null methylation map raise the documented
  exceptions. SHOULD — `HorvathAntiTransform` directly at `x = 0` and a strong-negative value. COULD —
  empty methylation map → age = `anti.trafo(intercept)` (intercept-only).

## Deviations and assumptions

- **RESOLVED (2026-06-22): the 353-CpG table is embedded.** Previously the coefficients were
  caller-supplied because they had not been retrieved; they were then retrieved verbatim from
  Additional file 3 and cross-verified byte-identical against the independent GitHub mirror (all 353
  pairs + intercept). Exposed as `HorvathMultiTissueCoefficients` / `HorvathMultiTissueIntercept` plus
  a parameterless `CalculateEpigeneticAge(methylation)` overload; the caller-supplied overload is
  retained. **No correctness-affecting assumption remains for the multi-tissue clock.**
- **Scope (not an assumption):** only the **multi-tissue** clock (353 CpGs) is in scope for this unit.
  The skin-&-blood clock (Horvath 2018) and PhenoAge (Levine 2018) are different models with different
  coefficient tables; callers needing them use the caller-supplied overload.

No source contradictions — the primary paper, the two reference R implementations, the Springer
supplement, its byte-identical mirror, and meffonym are mutually consistent on the two-stage model
(linear predictor → two-branch `anti.trafo`, `adult.age = 20`, intercept 0.695507258).

> Reconciliation note: the algorithm doc `docs/algorithms/Epigenetics/Epigenetic_Age_Estimation.md`
> (Last Reviewed 2026-06-23, one day after this Evidence file) has since **extended** the
> implementation to embed the Horvath-2018 skin-&-blood (391 CpG) and Levine-2018 PhenoAge (513 CpG,
> no transform) clocks with their own parameterless overloads — beyond this Evidence unit's
> multi-tissue scope. See [[epigenetic-age-horvath-clock]].
