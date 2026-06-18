# Validation Report: EPIGEN-AGE-001 — Epigenetic Age Estimation (Horvath DNA methylation clock)

- **Validated:** 2026-06-15   **Area:** Epigenetics
- **Canonical method(s):** `EpigeneticsAnalyzer.CalculateEpigeneticAge(methylationAtClockCpGs, coefficients, intercept)`, `EpigeneticsAnalyzer.HorvathAntiTransform(transformedAge)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (one test weakness fixed in-session; framework externalises the coefficient table)

## Stage A — Description

### Sources opened this session (independent of repo artifacts)
1. **aldringsvitenskap/epigeneticclock `horvath2013.R`** (fetched verbatim) —
   `https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/horvath2013.R`
   Confirms verbatim:
   - `anti.trafo = function(x, adult.age = 20) { ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age) }`
   - `trafo` forward transform; default `adult.age = 20`.
2. **aldringsvitenskap/epigeneticclock `StepwiseAnalysis.R`** (fetched verbatim) —
   `https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/StepwiseAnalysis.R`
   Confirms: `predictedAge = anti.trafo(CoefficientTraining[1] + as.matrix(datMethClock) %*% CoefficientTraining[-1])`
   i.e. age = `anti.trafo(intercept + Σ coef_i·β_i)`.
3. **Bioconductor wateRmelon `R/horv.R`** (independent 2nd implementation) —
   `https://rdrr.io/bioc/wateRmelon/src/R/horv.R`
   Same `trafo`/`anti.trafo` verbatim, `adult.age = 20`; age = `anti.trafo(data %*% coef2 + coeff[1])`,
   coefficients selected by intersection with available CpGs (→ non-clock CpGs ignored).
4. **perishky/meffonym** (independent 3rd implementation, via rdrr.io search) — same `trafo`/`anti.trafo`
   verbatim, `adult.age = 20`.
5. **Horvath S (2013) Genome Biology 14:R115** (PMC4015143) — 353 CpGs by elastic net; chronological age
   is regressed in a *transformed* form (log up to adulthood, linear after); the explicit transform lives
   in Additional file 2, reproduced by the reference implementations above.

### Formula check
- Inverse calibration `F⁻¹(Y)`: `21·exp(Y) − 1` for `Y < 0`, else `21·Y + 20`. Matches all THREE independent
  reference implementations byte-for-byte (sources #1, #3, #4). The `<` (strict) at the boundary places
  `Y = 0` on the linear branch → exactly 20 (adult age). ✓
- Linear predictor `Y = intercept + Σ coef_i·β_i`. Matches source #2 (and #3). ✓
- Non-clock CpG handling: only CpGs in the coefficient table contribute (source #2/#3 matrix product;
  source #3 explicit intersection). ✓

### Edge-case semantics check
- `Y = 0` → 20.0 (linear branch, `x < 0` false). Sourced (#1/#3). ✓
- `Y < 0` → exponential branch, < 20, → −1 as Y→−∞ (mathematical limit). Sourced (#1). ✓
- Empty methylation map → `F⁻¹(intercept)`; intercept always applies. Consistent with source #2. ✓
- null map / null coefficients / empty coefficients → exceptions: a defensible programming contract
  (an empty clock has no defined age); not a biological claim, no external source needed beyond contract. ✓

### Independent cross-check (hand-computed from the sourced formula)
| Input | Branch | Sourced expected |
|-------|--------|------------------|
| Y = 0.684247258 (=0.695507258+0.0127·0.5−0.0312·0.8+0.0245·0.3) | linear | 21·Y+20 = **34.369192418** |
| Y = 0 | linear | **20.0** |
| Y = −1.0 | exp | 21·e⁻¹−1 = **6.7254682646002895** |
| anti.trafo(−2.5) | exp | 21·e⁻²·⁵−1 = **0.7237849711018749** |
| anti.trafo(1.0) | linear | **41.0** |
| Y = 0.3 (M4: 0.1+0.5·0.4) | linear | 21·0.3+20 = **26.3** |
All recomputed in Python from the published formula; every value matches the spec/tests.

**Stage A findings:** Description (TestSpec, Evidence, algorithm doc) is correct and faithfully traces the
two-branch `anti.trafo` and the linear predictor to authoritative reference implementations. The
framework decision (caller supplies the 353-CpG table, not bundled/fabricated) is sound and documented.
No divergence. **PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs`
- `HorvathAntiTransform` (lines 1177–1182): `transformedAge < 0 ? (1+20)·exp(x)−1 : (1+20)·x+20`,
  `HorvathAdultAge = 20.0` (line 1123). Byte-for-byte the sourced `anti.trafo`. ✓
- `CalculateEpigeneticAge` (lines 1142–1168): null/null/empty guards (ArgumentNullException ×2,
  ArgumentException), `Y = intercept + Σ coef·β` over CpGs present in the coefficient table, then
  `HorvathAntiTransform(Y)`. Matches source #2/#3. ✓

### Cross-verification table recomputed vs code (full suite run)
All 12 unit tests pass; expected values are the externally hand-computed ones above (34.369192418, 20.0,
6.7254682646002895, 0.7237849711018749, 41.0, 26.3, and the ArgumentNull/ArgumentException contracts).

### Variant/delegate consistency
`CalculateEpigeneticAge` delegates the transform to `HorvathAntiTransform`; both tested directly and agree.

### Test quality audit (HARD gate)
- **Sourced, not code echoes:** M1/M2/M3/M5/S1/S2/C1 assert exact externally-derived constants. ✓
- **Defect found & fixed:** **M4** (non-clock CpG ignored) originally asserted only
  `ageWithExtra == ageWithoutExtra` — a relative check a *consistent-but-wrong* predictor (e.g. dropped
  intercept) would still pass. Strengthened to also assert the exact sourced value **26.3** (Y=0.3 →
  21·0.3+20), wrapped in `Assert.Multiple`. No weakening of any other assertion.
- **Coverage:** every public method exercised; both `anti.trafo` branches; boundary Y=0; intercept-only;
  monotonicity (INV-04 property, legitimately `GreaterThan` — a property, not a known single value); all
  three error contracts (null map, null coef, empty coef). ✓
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6542`; `dotnet build` 0 errors; the changed
  test file introduces no new warnings (the 4 build warnings are pre-existing, in unrelated files). ✓

### Findings / defects
- Stage-B defect (test weakness): M4 was a relative-equality assertion that would survive a deliberately
  wrong implementation. **Fixed in-session** (exact sourced value 26.3 added). No code defect found — the
  implementation matches the sourced formula exactly.
- Note (not a defect): `CalculateEpigeneticAge` is a framework — the 353-CpG Horvath coefficient table is
  caller-supplied, not bundled (per task policy forbidding fabricated coefficients). Every constant the
  code itself uses (adult.age=20, two-branch transform, predictor assembly) is source-backed.
- Out of scope: `PredictImprintedGenes` (listed under the unit in the registry) is unrelated to the age
  clock and has no retrieved authoritative basis; not modified, not tested here (own future unit).

## Verdict & follow-ups
- **Stage A: PASS.** Description matches three independent reference implementations verbatim.
- **Stage B: PASS-WITH-NOTES.** Code is correct; one test weakness (M4) found and completely fixed;
  framework coefficient externalisation documented and sound.
- **End-state: CLEAN.** Defect fully fixed in-session; build + full suite green (6542/0).
