# Validation Report: ONCO-EXPR-001 — Tumor Gene Expression Outlier (z-score) and Signature Score

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CalculateExpressionZScore(double, IReadOnlyList<double>)`, `OncologyAnalyzer.IdentifyOutlierGenes(IReadOnlyDictionary<string,double>, IReadOnlyDictionary<string,IReadOnlyList<double>>, double)`, `OncologyAnalyzer.CalculateSignatureScore(IReadOnlyList<double>)` (+ `ExpressionDirection` enum, `ExpressionOutlier` record)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** CLEAN

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **cBioPortal z-score normalization spec** — https://docs.cbioportal.org/z-score-normalization-script/ (WebFetch). Confirms verbatim the z-score formula `(r - mu)/sigma` where `r` is the raw expression value and `mu`/`sigma` are the mean and standard deviation of the reference base population. Confirms zero-SD handling: "Z-Score <- NA when standard deviation = 0". The prose does **not** specify n vs n−1 for sigma (resolved by the reference implementation, below).
2. **cBioPortal-core `NormalizeExpressionLevels.java`** — fetched live via `gh api repos/cBioPortal/cbioportal-core/contents/.../NormalizeExpressionLevels.java` (base64-decoded). Inspected the actual source:
   - `avg`: `avg = avg/(double)v.length` → arithmetic mean (lines 596–603).
   - `std`: `std=std/(double)(v.length-1); std=Math.sqrt(std)` → **sample standard deviation, divisor (n−1)** (lines 605–613). This authoritatively resolves the n-vs-(n−1) question the prose spec leaves open.
   - `getZ`: `if (0.0d == std) fatalError("cannot normalize relative to distribution with standard deviation of 0.0.")` (lines 426–432) → zero-SD is a hard error, not a silent value.
3. **cBioPortal FAQ** (corroborated via WebSearch) — default outlier threshold ±2: samples with z-scores ">2 or <-2 … are considered altered"; z>2 ⇒ over, z<−2 ⇒ under; ±2 is a *commonly recommended default*, strict inequality.
4. **Lee E. et al. (2008), PLoS Comput Biol 4(11):e1000217** — https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000217 (WebFetch). Confirms verbatim: per-gene z-scores `z_ij` are standardized to "mean μi = 0 and standard deviation σi = 1 over all samples j"; and the combined score: "The individual z_ij of each member gene in the gene set are averaged into a combined z-score … (the square root of the number of member genes is used in the denominator to stabilize the variance of the mean)" ⇒ **a = (Σ z_i)/√k**.
5. **GSVA vignette** (WebSearch corroboration) — lists the "combined z-score" method attributed to Lee et al. 2008: per-gene z then Σz/√k.

### Formula check
- z = (r − μ)/σ — matches source 1 verbatim. ✓
- σ = sqrt(Σ(rᵢ−μ)²/(n−1)) — matches reference impl `std()` (source 2, divisor n−1). ✓
- Outlier: z > +t over / z < −t under, strict, default t=2 — matches sources 3. ✓
- Signature a = Σz/√k — matches source 4 verbatim. ✓

### Edge-case semantics
- σ = 0 (constant cohort) → error: matches `fatalError` in reference impl (source 2). ✓
- n ≤ 1 → sample SD undefined (divisor n−1=0): defined as error. Consistent with the (n−1) definition; reasonable, sourced inference. ✓
- z = ±t exactly → not an outlier (strict >/<): matches FAQ wording ">2 or <-2". ✓
- k = 1 → a = z₁/√1 = z₁. ✓; k = 0 → undefined → error. ✓

### Independent cross-check (hand computation, traced to sources)
Cohort {2,2,4,6,6}: μ=(2+2+4+6+6)/5=4; Σ(rᵢ−μ)²=4+4+0+4+4=16; var=16/(5−1)=4; σ=2.
- z(10)=(10−4)/2=**3.0**; z(4)=**0.0**; z(−1)=**−2.5**; z(6)=(6−4)/2=**1.0** (population-SD would give 1.118 — distinguishing test); z(8)=(8−4)/2=**2.0** (boundary, not outlier); z(−2)=**−3.0** (reflection of z(10)).
- Signature {3,1,−1,1}: Σ=4, √4=2, a=**2.0**; {2.5}: a=**2.5**; {1.5,1.5,1.5,1.5}: Σ=6, a=6/2=**3.0**.

All values trace to the external formulas (sources 1,2,4), not to the implementation.

### Findings / divergences
None. The cBioPortal prose alone does not fix the SD divisor, but the cited reference implementation does (n−1); the Evidence/TestSpec correctly rely on it. The ±2 threshold is documented as a *default* (caller-overridable), which the API honours via the `threshold` parameter.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:7436–7613`.

### Formula realised correctly?
- `CalculateExpressionZScore` (7478–7518): mean over n, then `sd = sqrt(sumSquaredDeviations/(n-1))` — sample SD (n−1), matching source 2 exactly; throws `ArgumentException` for n<2 and for sd==0.0; null guard. ✓
- `IdentifyOutlierGenes` (7537–7575): per-gene z via the canonical method; strict `z > threshold` ⇒ Over, `z < -threshold` ⇒ Under; null guards on both dictionaries; `threshold <= 0` ⇒ `ArgumentOutOfRangeException`; missing reference cohort ⇒ `ArgumentException`. Iteration order = sample dictionary order. ✓
- `CalculateSignatureScore` (7593–7611): null guard; k==0 ⇒ `ArgumentException`; returns `sum / Math.Sqrt(k)`. ✓

### Cross-verification table recomputed vs code (via passing tests)
| Input | Expected (source) | Code result |
|-------|-------------------|-------------|
| z(10), {2,2,4,6,6} | 3.0 | 3.0 ✓ |
| z(4) | 0.0 | 0.0 ✓ |
| z(−1) | −2.5 | −2.5 ✓ |
| z(6) (n−1 distinguishing) | 1.0 (not 1.118) | 1.0 ✓ |
| z(8) boundary | not outlier | excluded ✓ |
| signature {3,1,−1,1} | 2.0 | 2.0 ✓ |
| signature {2.5} | 2.5 | 2.5 ✓ |
| signature {1.5×4} | 3.0 | 3.0 ✓ |
| constant cohort {5,5,5} | throw | ArgumentException ✓ |

### Variant/delegate consistency
`IdentifyOutlierGenes` delegates z-scoring to `CalculateExpressionZScore` (single source of truth). No `*Fast`/instance duplicates exist. Consistent. ✓

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** every expected value derives from z=(x−μ)/σ and a=Σz/√k computed by hand from the cited sources. M4 deliberately picks z(6)=1.0 so a population-SD (n) implementation would fail (1.118) — guards the documented pitfall.
- **No green-washing:** exact-value asserts with 1e-10 floating tolerance only; no Greater/AtLeast/Contains/ranges where an exact value is known; no skipped/ignored tests; no widened tolerances.
- **Coverage:** all three public methods + both result types exercised. Branches: happy paths (M1–M3, M9, M10), n−1 vs n distinction (M4), z=0 (M2/INV-01), reflection (C1/INV-03), equal-z variance stabilization (C2/INV-05), single-gene k=1 (M10/INV-06); error/edge: null cohort (S1), n<2 (S2), σ=0 (M11), strict boundary (M7), over/under/normal classification (M5/M6/M8), null sample/cohorts (S5), missing cohort (S6), non-positive threshold, empty signature (S3), null signature (S4). Every Stage-A edge case is covered.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6694`; changed files: none (no code/test change required); build 0 errors.
- Note: TestSpec §5.4 labels "19" tests; the file has 18 methods (M5/M6/M8 are combined into one multi-assert test, and a non-positive-threshold guard test was added). Coverage is complete; the count label is cosmetic, not a defect.

**Test-quality gate: PASS.**

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS.** Description matches cBioPortal spec + reference implementation (sample SD n−1, σ=0 error, strict ±2 default) and Lee et al. 2008 (a=Σz/√k), all independently retrieved this session.
- **Stage B: PASS.** Code realises the validated formulas exactly; all edge cases handled; tests assert exact sourced values and cover every branch.
- **End-state: CLEAN.** No defect found; no code or test change needed. Full suite green (6694 passed, 0 failed); build clean.
- Follow-ups: none.
