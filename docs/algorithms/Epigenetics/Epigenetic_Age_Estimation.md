# Epigenetic Age Estimation (Horvath DNA Methylation Clock)

| Field | Value |
|-------|-------|
| Algorithm Group | Epigenetics |
| Test Unit ID | EPIGEN-AGE-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Estimates DNA methylation ("epigenetic") age from methylation β-values measured at clock CpG sites,
following the Horvath (2013) multi-tissue epigenetic clock [1]. The estimate is a two-stage,
specification-driven computation: a linear predictor `Y = intercept + Σ coef_i · β_i` over the clock
CpGs, followed by the Horvath inverse calibration `F⁻¹` that maps transformed age to years [2][3].
Because the clock's coefficient table (353 CpGs) is a large published table [1], the implementation is a
**framework**: it performs the source-defined math but the caller supplies the coefficient table and
intercept for the clock they intend to use.

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
dependence is "logarithmic until adulthood … linear later in life" [1].

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
| coefficients | `IReadOnlyDictionary<string,double>` | required | CpG id → clock coefficient | non-null, non-empty |
| intercept | `double` | 0.0 | Model intercept added before the inverse transform | finite |

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
- The clock coefficient table itself (353 CpGs for Horvath) is **not bundled**; it is supplied by the
  caller. See §5.3 "Not implemented".

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateEpigeneticAge | O(n) | O(1) | n = number of methylation entries; one dictionary lookup each |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [EpigeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs)

- `EpigeneticsAnalyzer.CalculateEpigeneticAge(methylationAtClockCpGs, coefficients, intercept)`: computes the linear predictor and applies the inverse transform.
- `EpigeneticsAnalyzer.HorvathAntiTransform(transformedAge)`: the published two-branch `anti.trafo`.

### 5.2 Current Behavior

- No search/matching is performed, so the repository suffix tree is **not applicable** (this is a weighted
  sum over a dictionary, not occurrence finding).
- Coefficients are caller-supplied; there is no built-in default table. A null/empty coefficient table is
  rejected rather than silently substituted.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Linear predictor `Y = intercept + Σ coef_i · β_i` over clock CpGs [3].
- Inverse calibration `F⁻¹`: `21·exp(Y)−1` for `Y < 0`, `21·Y + 20` otherwise, with adult.age = 20 [2].

**Intentionally simplified:**

- (none)

**Not implemented:**

- The 353-CpG Horvath coefficient table and intercept; **users should rely on:** supplying the published
  coefficients (Additional file 3 of [1]) or another clock's table via the `coefficients`/`intercept`
  parameters. Bundling fabricated coefficients is a defect and is deliberately avoided.
- BMIQ "gold standard" normalisation of input β-values that the full Horvath pipeline performs upstream;
  **users should rely on:** normalising inputs before calling this method.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Coefficient table externalised | Assumption | Output meaningfulness depends on caller-supplied weights | accepted | ASM-01; framework design per task policy on large published tables |
| 2 | Prior code used `exp(x)−1` single-branch transform and hard-coded example coefficients | Deviation (fixed) | Produced incorrect ages | fixed | Replaced with two-branch `anti.trafo` + caller coefficients |

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

- Not a self-contained clock: without a valid coefficient table the output is not an age.
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

### 7.3 Related Tests, Evidence, or Documents

- Tests: [EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [EPIGEN-AGE-001-Evidence.md](../../../docs/Evidence/EPIGEN-AGE-001-Evidence.md)
- Related algorithms: [Methylation_Analysis](./Methylation_Analysis.md)

## 8. References

1. Horvath S. 2013. DNA methylation age of human tissues and cell types. Genome Biology 14:R115. https://doi.org/10.1186/gb-2013-14-10-r115
2. aldringsvitenskap/epigeneticclock. 2013-method reference R implementation, `horvath2013.R` (trafo / anti.trafo, adult.age = 20). https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/horvath2013.R
3. aldringsvitenskap/epigeneticclock. 2013-method reference R implementation, `StepwiseAnalysis.R` (predictedAge = anti.trafo(intercept + meth·coef)). https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/StepwiseAnalysis.R
4. perishky/meffonym. R package implementing DNA methylation clocks (independent anti.trafo confirmation). https://github.com/perishky/meffonym
