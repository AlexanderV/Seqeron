# Epigenetic Age Estimation (Horvath DNA Methylation Clock)

| Field | Value |
|-------|-------|
| Algorithm Group | Epigenetics |
| Test Unit ID | EPIGEN-AGE-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Production |
| Last Reviewed | 2026-06-23 |

## 1. Overview

Estimates DNA methylation ("epigenetic") age from methylation β-values measured at clock CpG sites,
following three published epigenetic clocks: the Horvath (2013) multi-tissue clock [1], the Horvath
(2018) skin &amp; blood clock [5], and the Levine (2018) DNAm PhenoAge clock [6]. The Horvath clocks are a
two-stage computation — a linear predictor `Y = intercept + Σ coef_i · β_i` over the clock CpGs, followed
by the Horvath inverse calibration `F⁻¹` (`anti.trafo`, adult.age = 20) that maps transformed age to years
[2][3]. PhenoAge is a one-stage linear predictor in years with **no** transform [6]. All three coefficient
tables (353, 391, and 513 CpGs) and their intercepts are **embedded** from the papers' supplements and
cross-verified against an independent reimplementation, so a caller can compute each age directly via a
parameterless overload. Caller-supplied-coefficients overloads are also retained for other clocks.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

DNA methylation at CpG dinucleotides changes systematically with age. Horvath (2013) selected 353 CpGs
by elastic-net regression of a *transformed* version of chronological age, yielding a predictor that
generalises across tissues [1]. Methylation is measured as a β-value in [0, 1] (fraction methylated).

### 2.2 Core Model

Let `β_i` be the methylation value at clock CpG `i` with coefficient `coef_i`, and `intercept` the model
intercept. The linear predictor in transformed-age units is [3]:

> `Y = intercept + Σ_i (coef_i · β_i)`

DNAm age in years is the inverse calibration `F⁻¹(Y)`, with adult-age constant `adult.age = 20` [2]:

> `F⁻¹(Y) = (1 + 20)·exp(Y) − 1`   if `Y < 0`
> `F⁻¹(Y) = (1 + 20)·Y + 20`        if `Y ≥ 0`

This is the `anti.trafo` function from the Horvath 2013 reference R code [2]; the forward calibration `F`
is logarithmic below adult age and linear above, reflecting the paper's observation that methylation–age
dependence is "logarithmic until adulthood … linear later in life" [1]. The **Horvath (2018) skin &amp;
blood clock** uses the SAME `anti.trafo` (adult.age = 20) with its own 391-CpG table and intercept
−0.447119319 [5][7]. The **Levine (2018) PhenoAge clock** is different: DNAm PhenoAge = `intercept + Σ
weight_i · β_i` is returned directly with **no** transform — the linear predictor is already in years [6].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Coefficients and intercept come from a validly trained clock on comparable normalised data | Output is not interpretable as age; arbitrary numbers in, arbitrary numbers out |
| ASM-02 | Methylation β-values are on [0, 1] and measured at the same CpGs the clock was trained on | Mis-scaled predictor; biased age |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | For `Y ≥ 0`, age = `21·Y + 20` | Linear branch of `anti.trafo` [2] |
| INV-02 | For `Y < 0`, age = `21·exp(Y) − 1` | Exponential branch of `anti.trafo` [2] |
| INV-03 | At `Y = 0`, age = 20.0 exactly | `Y < 0` false → linear branch → `21·0 + 20` [2] |
| INV-04 | Age is strictly increasing and continuous in `Y` | Both branches strictly increasing; meet at `(0, 20)` [2] |
| INV-05 | CpGs absent from the coefficient table do not change the result | Only clock CpGs enter the weighted sum [3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| methylationAtClockCpGs | `IReadOnlyDictionary<string,double>` | required | CpG id → methylation β-value | values typically [0,1]; non-null |
| coefficients | `IReadOnlyDictionary<string,double>` | required (explicit overload) | CpG id → clock coefficient | non-null, non-empty |
| intercept | `double` | 0.0 (explicit overload) | Model intercept added before the inverse transform | finite |

The parameterless overloads use embedded tables: `CalculateEpigeneticAge(methylation)` →
`HorvathMultiTissueCoefficients` / intercept 0.695507258; `CalculateSkinBloodAge(methylation)` →
`HorvathSkinBloodCoefficients` / intercept −0.447119319 (anti.trafo path); `CalculatePhenoAge(methylation)`
→ `PhenoAgeCoefficients` / intercept 60.664 (no transform).

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `double` | Estimated DNAm age in years (`F⁻¹` of the linear predictor) |

### 3.3 Preconditions and Validation

- `methylationAtClockCpGs == null` → `ArgumentNullException`.
- `coefficients == null` → `ArgumentNullException`.
- `coefficients.Count == 0` → `ArgumentException` (an empty clock has no defined output).
- An empty methylation map is valid: the result is `F⁻¹(intercept)` (no CpG contributions).
- CpG ids are matched by exact dictionary key; case sensitivity follows the supplied dictionary's comparer.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs (null map, null/empty coefficients).
2. Initialise the linear predictor `Y` to `intercept`.
3. For each `(cpg, β)` in the methylation map, if `cpg` is in the coefficient table add `coef·β` to `Y`.
4. Return `HorvathAntiTransform(Y)` (the two-branch `F⁻¹`).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **adult.age = 20**: the calibration break between the logarithmic and linear regimes [2].
- **Intercept = 0.695507258**: `CoefficientTraining[1]` of the embedded multi-tissue clock [4].
- The 353-CpG Horvath multi-tissue coefficient table is **embedded** (`HorvathMultiTissueCoefficients`),
  retrieved from Additional file 3 [4] and cross-verified byte-identical against an independent mirror
  (Evidence sources #5, #6). Other clocks' tables are still supplied via the caller overload.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateEpigeneticAge | O(n) | O(1) | n = number of methylation entries; one dictionary lookup each |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [EpigeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs) and the embedded table [EpigeneticsAnalyzer.HorvathClock.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.HorvathClock.cs)

- `EpigeneticsAnalyzer.CalculateEpigeneticAge(methylationAtClockCpGs)`: computes DNAm age using the built-in Horvath multi-tissue clock.
- `EpigeneticsAnalyzer.CalculateEpigeneticAge(methylationAtClockCpGs, coefficients, intercept)`: computes the linear predictor and applies the inverse transform with caller-supplied coefficients.
- `EpigeneticsAnalyzer.HorvathMultiTissueCoefficients` / `HorvathMultiTissueIntercept`: the embedded 353-CpG table and intercept.
- `EpigeneticsAnalyzer.HorvathAntiTransform(transformedAge)`: the published two-branch `anti.trafo`.

### 5.2 Current Behavior

- No search/matching is performed, so the repository suffix tree is **not applicable** (this is a weighted
  sum over a dictionary, not occurrence finding).
- The Horvath multi-tissue coefficient table is built in; the parameterless overload computes age directly.
  The explicit overload still accepts caller coefficients for other clocks, and a null/empty coefficient
  table there is rejected rather than silently substituted.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Linear predictor `Y = intercept + Σ coef_i · β_i` over clock CpGs [3].
- Inverse calibration `F⁻¹`: `21·exp(Y)−1` for `Y < 0`, `21·Y + 20` otherwise, with adult.age = 20 [2].
- The published 353-CpG multi-tissue coefficient table and intercept (0.695507258), embedded verbatim from
  Additional file 3 [4] (cross-verified byte-identical against an independent mirror).

**Intentionally simplified:**

- (none)

**Not implemented:**

- BMIQ "gold standard" normalisation of input β-values that the full Horvath pipeline performs upstream;
  **users should rely on:** normalising inputs before calling this method.
- Other Horvath clocks (skin-&-blood 2018, PhenoAge 2018) have different coefficient tables and are not
  embedded; **users should rely on:** the caller-supplied-coefficients overload for those models.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | 353-CpG coefficient table embedded | Resolved | Multi-tissue clock computes age directly; no caller table needed | resolved 2026-06-22 | Retrieved from Additional file 3 [4], cross-verified byte-identical (Evidence #5,#6) |
| 2 | Prior code used `exp(x)−1` single-branch transform and hard-coded example coefficients | Deviation (fixed) | Produced incorrect ages | fixed | Replaced with two-branch `anti.trafo` + published coefficients |
| 3 | Only the multi-tissue clock is embedded | Scope | Skin-&-blood / PhenoAge need the caller overload | accepted | Out of scope for EPIGEN-AGE-001 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `Y = 0` | age = 20.0 | Linear branch at boundary [2] |
| `Y < 0` | age = `21·exp(Y)−1` (< 20) | Exponential branch [2] |
| Empty methylation map | age = `F⁻¹(intercept)` | No CpG contributions; intercept still applies [3] |
| CpG not in coefficient table | ignored | Only clock CpGs enter the sum [3] |
| null map / null coefficients | `ArgumentNullException` | Contract |
| empty coefficients | `ArgumentException` | Empty clock undefined |

### 6.2 Limitations

- Only the Horvath multi-tissue clock is built in; other clocks require caller-supplied coefficients.
- No input normalisation (BMIQ) or array-platform handling; callers must pre-normalise.
- The transform is biologically meaningful for realistic predictors; for extreme negative `Y` the formula
  approaches −1 (a mathematical artefact, not a biological age).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var coefficients = new Dictionary<string, double>
{
    ["cg00000029"] = 0.0127,
    ["cg00000165"] = -0.0312,
    ["cg00000363"] = 0.0245,
};
var methylation = new Dictionary<string, double>
{
    ["cg00000029"] = 0.5,
    ["cg00000165"] = 0.8,
    ["cg00000363"] = 0.3,
};
// intercept 0.695507258 → Y = 0.684247258 → age = 21*Y + 20 = 34.369192418
double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, coefficients, 0.695507258);
```

**Numerical walk-through:** `Y = 0.695507258 + 0.0127·0.5 + (−0.0312)·0.8 + 0.0245·0.3 = 0.684247258`;
`Y ≥ 0` → age `= 21·0.684247258 + 20 = 34.369192418` years.

**Built-in clock example:**

```csharp
// Uses the embedded Horvath 2013 353-CpG multi-tissue table + intercept 0.695507258.
var betas = new Dictionary<string, double> { ["cg00864867"] = 1.0 };
// Y = 0.695507258 + 1.59976405·1.0 = 2.295271308 → 21·Y + 20 = 68.200697468
double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(betas);
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [EPIGEN-AGE-001-Evidence.md](../../../docs/Evidence/EPIGEN-AGE-001-Evidence.md)
- Related algorithms: [Methylation_Analysis](./Methylation_Analysis.md)

## 8. References

1. Horvath S. 2013. DNA methylation age of human tissues and cell types. Genome Biology 14:R115. https://doi.org/10.1186/gb-2013-14-10-r115
2. aldringsvitenskap/epigeneticclock. 2013-method reference R implementation, `horvath2013.R` (trafo / anti.trafo, adult.age = 20). https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/horvath2013.R
3. aldringsvitenskap/epigeneticclock. 2013-method reference R implementation, `StepwiseAnalysis.R` (predictedAge = anti.trafo(intercept + meth·coef)). https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/StepwiseAnalysis.R
4. Horvath S. 2013. Additional file 3 — 353-CpG `CoefficientTraining` table (intercept 0.695507258 + 353 weights), Genome Biology 14:R115. https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2013-14-10-r115/MediaObjects/13059_2013_3156_MOESM3_ESM.csv (cross-verified against https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/AdditionalFile3.csv)
5. perishky/meffonym. R package implementing DNA methylation clocks (independent anti.trafo confirmation). https://github.com/perishky/meffonym
