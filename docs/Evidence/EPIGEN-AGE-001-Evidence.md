# Evidence Artifact: EPIGEN-AGE-001

**Test Unit ID:** EPIGEN-AGE-001
**Algorithm:** Epigenetic Age Estimation (Horvath DNA methylation clock)
**Date Collected:** 2026-06-13

---

## Online Sources

### Horvath S (2013) — DNA methylation age of human tissues and cell types (primary paper)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4015143/
**Accessed:** 2026-06-13 (fetched via WebFetch; original NCBI PMC URL https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4015143/ 301-redirected to pmc.ncbi.nlm.nih.gov, then fetched)
**Authority rank:** 1 (peer-reviewed paper, Genome Biology)

**Key Extracted Points:**

1. **Clock size:** The multi-tissue age predictor uses **353 CpGs** automatically selected by elastic-net regression (193 positively correlated, 160 negatively correlated with age). Extracted from the fetched Methods/Results text.
2. **Penalized regression:** A *transformed* version of chronological age is regressed on the CpGs using elastic net (alpha = 0.5, lambda = 0.0226). The predictor is therefore a linear combination of CpG methylation values plus an intercept, in transformed-age units.
3. **Calibration shape:** The paper states the calibration curve shows "a logarithmic dependence until adulthood that slows to a linear dependence later in life" — i.e. a piecewise log/linear transform with a break at adult age. The explicit formula for transform F and its inverse is in Additional file 2 (not reproduced in the main HTML), so the exact two-branch formula was taken from reference implementations below.

### Horvath 2013 reference R implementation (aldringsvitenskap/epigeneticclock)

**URL:** https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/horvath2013.R
**Accessed:** 2026-06-13 (fetched verbatim via `curl -sL` and via WebFetch)
**Authority rank:** 3 (reference implementation of the Horvath 2013 method)

**Key Extracted Points (verbatim R code):**

1. **adult.age constant = 20** (default parameter of both transform functions).
2. **Forward transform `trafo`:**
   ```r
   trafo = function(x, adult.age = 20) {
     x = (x + 1) / (1 + adult.age)
     y = ifelse(x <= 1, log(x), x - 1)
     y
   }
   ```
3. **Inverse transform `anti.trafo`:**
   ```r
   anti.trafo = function(x, adult.age = 20) {
     ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
   }
   ```
   i.e. DNAm age = `(1+20)·exp(x) − 1` when `x < 0`, else `(1+20)·x + 20`.

### Horvath 2013 reference R implementation — age computation (StepwiseAnalysis.R)

**URL:** https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/StepwiseAnalysis.R
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Predicted age formula (verbatim):**
   ```r
   predictedAge = as.numeric(anti.trafo(
     datClock$CoefficientTraining[1] + as.matrix(datMethClock) %*% as.numeric(datClock$CoefficientTraining[-1])))
   ```
   The **intercept** is `CoefficientTraining[1]`; the **per-CpG coefficients** are `CoefficientTraining[-1]`; the linear predictor `intercept + Σ coef_i · β_i` is passed to `anti.trafo` to yield DNAm age in years.

### Horvath 2013 Additional file 3 — the 353-CpG coefficient table (PRIMARY supplement)

**URL:** https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2013-14-10-r115/MediaObjects/13059_2013_3156_MOESM3_ESM.csv
**Accessed:** 2026-06-22 (fetched via WebFetch and downloaded via `curl -sL`)
**Authority rank:** 1 (official supplement of the peer-reviewed paper, Genome Biology 14:R115)

**Key Extracted Points:**

1. **Columns:** `CpGmarker, CoefficientTraining, CoefficientTrainingShrunk, varByCpG, ...` (23 columns). The age-relating weights are in `CoefficientTraining`; the file header states "These coefficient values relate CpGs to a transformed version of age".
2. **Intercept:** the `(Intercept)` row has `CoefficientTraining = 0.695507258`.
3. **353 CpG rows:** rows 5..357 of the CSV are the 353 probe coefficients (e.g. cg00075967 = 0.12933661, cg00374717 = 0.005017857, cg00864867 = 1.59976405, cg09809672 (EDARADD) = -0.391318905, last row cg27544190 = -0.869124446).
4. **Note (memory-trap probe):** cg22454769 is **not** in the multi-tissue 353-CpG table — confirmed absent in the retrieved CSV. (Recorded because the task warned not to trust that probe from memory.)

### Cross-verification source — aldringsvitenskap/epigeneticclock AdditionalFile3.csv (GitHub mirror)

**URL:** https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/AdditionalFile3.csv
**Accessed:** 2026-06-22 (fetched via WebFetch and `curl -sL`)
**Authority rank:** 3 (faithful mirror of the supplement bundled with a reference R implementation)

**Key Extracted Points:**

1. Same header columns and same `(Intercept) = 0.695507258`.
2. A complete field-by-field diff of all 353 `(CpGmarker, CoefficientTraining)` pairs against the Springer supplement showed **zero differences** — the two independent sources are byte-identical. This is the cross-check that authorised embedding the table.

### Second independent reference implementation (meffonym, perishky/meffonym)

**URL:** https://github.com/perishky/meffonym (searched + fetched R/horvath.r)
**Accessed:** 2026-06-13 (WebSearch + WebFetch)
**Authority rank:** 3 (reference implementation, Bioconductor-style R package)

**Key Extracted Points:**

1. Confirms the Horvath DNAmAge model applies an `anti.trafo` inverse log-transform (default behaviour of `meffonym.score()` for the Horvath model) to the linear model estimate to convert child/young ages, matching the two-branch transform above. Corroborates source #2/#3.

---

## Documented Corner Cases and Failure Modes

### From the reference R implementations (sources #2, #3)

1. **Two-branch transform at x = 0:** When the linear predictor `x = 0`, the `x < 0` branch is false, so the linear branch applies: `(1+20)·0 + 20 = 20`. The boundary at x = 0 evaluates to exactly the adult age (20).
2. **Negative linear predictor (young / hypomethylated profiles):** When `x < 0`, the exponential branch applies, producing ages below 20 that approach −1 as `x → −∞` (mathematical limit; biologically ages stay positive for realistic inputs).
3. **CpGs without a coefficient:** Only CpGs present in the coefficient table contribute (the R code multiplies the methylation matrix by the clock coefficient vector; CpGs not in the clock are not part of `datMethClock`). Methylation values for non-clock CpGs are ignored.

---

## Test Datasets

### Dataset: Synthetic worked examples derived from the published transform

**Source:** Horvath (2013) reference R code `anti.trafo` (source #2) — exact arithmetic computed from the published formula.

| Linear predictor Y | Branch | DNAm age = anti.trafo(Y) |
|--------------------|--------|--------------------------|
| 0.0 | linear (Y ≥ 0) | 20.0 |
| 1.0 | linear | 41.0 |
| −1.0 | exponential | 21·e⁻¹ − 1 = 6.7254682646002895 |
| −2.5 | exponential | 21·e⁻²·⁵ − 1 = 0.7237849711018749 |

### Dataset: Linear-predictor assembly (intercept + Σ coef·β)

**Source:** Horvath (2013) reference R code `StepwiseAnalysis.R` (source #3).

| Component | Value |
|-----------|-------|
| intercept | 0.695507258 |
| (coef, β) pairs | (0.0127, 0.5), (−0.0312, 0.8), (0.0245, 0.3) |
| Y = intercept + Σ coef·β | 0.684247258 |
| DNAm age = 21·Y + 20 (linear branch) | 34.369192418 |

---

## Assumptions

1. **RESOLVED (2026-06-22): the 353-CpG coefficient table is now embedded.** Previously the table was externalised to the caller because it had not been retrieved. It has since been retrieved verbatim from Horvath (2013) Additional file 3 (Springer supplement, source #5) and cross-verified byte-identical against an independent GitHub mirror (source #6) for all 353 `(CpGmarker, CoefficientTraining)` pairs plus the intercept (0.695507258). The built-in clock is exposed via `HorvathMultiTissueCoefficients` / `HorvathMultiTissueIntercept` and a parameterless `CalculateEpigeneticAge(methylation)` overload; the caller-supplied overload is retained unchanged. No correctness-affecting assumption remains for the multi-tissue clock.
2. **Scope note (not an assumption):** only the **multi-tissue** clock (353 CpGs) is embedded. The skin-&-blood clock (Horvath 2018) and PhenoAge (Levine 2018) are different models with different coefficient tables and are out of scope for this unit; callers needing them use the caller-supplied overload.

---

## Recommendations for Test Coverage

1. **MUST Test:** Linear branch with intercept + multiple CpGs gives `21·Y + 20` — Evidence: sources #2, #3 (worked example Y = 0.684247258 → 34.369192418).
2. **MUST Test:** Boundary Y = 0 → exactly 20.0 (adult age) — Evidence: source #2 (`x < 0` false → linear branch).
3. **MUST Test:** Negative branch Y < 0 → `21·exp(Y) − 1` (Y = −1 → 6.7254682646002895) — Evidence: source #2.
4. **MUST Test:** CpGs absent from the coefficient table are ignored — Evidence: source #3 (only clock CpGs enter the matrix product).
5. **MUST Test:** Null/empty coefficient table and null methylation map raise the documented exceptions — Evidence: contract/precondition (empty clock has no defined output).
6. **SHOULD Test:** `HorvathAntiTransform` directly at the x = 0 boundary and a strong-negative value — Rationale: isolates the published transform from the predictor assembly.
7. **COULD Test:** Empty methylation map with non-empty coefficients → age = anti.trafo(intercept) (no CpG contributions) — Rationale: confirms intercept-only behaviour.

---

## References

1. Horvath S. (2013). DNA methylation age of human tissues and cell types. Genome Biology 14:R115. https://doi.org/10.1186/gb-2013-14-10-r115 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC4015143/)
2. aldringsvitenskap/epigeneticclock. Horvath 2013 reference R implementation, `horvath2013.R` (trafo / anti.trafo, adult.age = 20). https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/horvath2013.R (accessed 2026-06-13)
3. aldringsvitenskap/epigeneticclock. Horvath 2013 reference R implementation, `StepwiseAnalysis.R` (predictedAge = anti.trafo(intercept + meth %*% coef)). https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/StepwiseAnalysis.R (accessed 2026-06-13)
4. perishky/meffonym. R package implementing DNA methylation clocks incl. Horvath anti.trafo. https://github.com/perishky/meffonym (accessed 2026-06-13)
5. Horvath S. (2013). Additional file 3 (353-CpG `CoefficientTraining` table), Genome Biology 14:R115. https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2013-14-10-r115/MediaObjects/13059_2013_3156_MOESM3_ESM.csv (accessed 2026-06-22)
6. aldringsvitenskap/epigeneticclock. AdditionalFile3.csv (GitHub mirror of the supplement; cross-verification source). https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/AdditionalFile3.csv (accessed 2026-06-22)

---

## Change History

- **2026-06-13**: Initial documentation.
- **2026-06-22**: Embedded the published Horvath (2013) 353-CpG multi-tissue clock coefficients (intercept 0.695507258 + 353 weights from Additional file 3), cross-verified byte-identical against an independent GitHub mirror; resolved the "coefficients are caller-supplied" assumption (sources #5, #6).
