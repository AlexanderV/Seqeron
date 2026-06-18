# Test Specification: EPIGEN-AGE-001

**Test Unit ID:** EPIGEN-AGE-001
**Area:** Epigenetics
**Algorithm:** Epigenetic Age Estimation (Horvath DNA methylation clock)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Horvath S. (2013). DNA methylation age of human tissues and cell types. Genome Biology 14:R115. | 1 | https://doi.org/10.1186/gb-2013-14-10-r115 / https://pmc.ncbi.nlm.nih.gov/articles/PMC4015143/ | 2026-06-13 |
| 2 | Horvath 2013 reference R, `horvath2013.R` (trafo / anti.trafo, adult.age=20) | 3 | https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/horvath2013.R | 2026-06-13 |
| 3 | Horvath 2013 reference R, `StepwiseAnalysis.R` (predictedAge = anti.trafo(intercept + meth·coef)) | 3 | https://raw.githubusercontent.com/aldringsvitenskap/epigeneticclock/master/StepwiseAnalysis.R | 2026-06-13 |
| 4 | perishky/meffonym (independent anti.trafo confirmation) | 3 | https://github.com/perishky/meffonym | 2026-06-13 |

### 1.2 Key Evidence Points

1. DNAm age = `anti.trafo(intercept + Σ coef_i · β_i)` over clock CpGs — source #3.
2. `anti.trafo(x) = (1+20)·exp(x) − 1` if `x < 0`, else `(1+20)·x + 20`; adult.age = 20 — source #2.
3. The clock uses 353 elastic-net-selected CpGs; coefficient table is large/published, not bundled here — source #1.
4. Only CpGs present in the coefficient table contribute to the sum — source #3.

### 1.3 Documented Corner Cases

- Boundary at linear predictor x = 0: `x < 0` is false → linear branch → `21·0 + 20 = 20` (source #2).
- Negative predictor → exponential branch, ages below 20 (source #2).
- Non-clock CpGs ignored (source #3).

### 1.4 Known Failure Modes / Pitfalls

1. Using an `exp(x) − 1` single-branch transform (the prior repo code) instead of the two-branch `anti.trafo` is incorrect — source #2.
2. Omitting the intercept term mis-scales every result — source #3.
3. Fabricated coefficient tables produce meaningless ages — source #1 (coefficients are the model).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateEpigeneticAge(methylationAtClockCpGs, coefficients, intercept)` | EpigeneticsAnalyzer | Canonical | Linear predictor + Horvath inverse transform |
| `HorvathAntiTransform(transformedAge)` | EpigeneticsAnalyzer | Canonical | Published two-branch `anti.trafo` (also exercised directly) |
| `PredictImprintedGenes(...)` | EpigeneticsAnalyzer | Out of scope | Listed under unit in Registry but unrelated to the age clock; its ad-hoc imprinting score has no retrieved authoritative basis — see §7. Not modified, not tested here. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | For linear predictor Y ≥ 0, age = `(1+20)·Y + 20` (linear branch) | Yes | Source #2 |
| INV-2 | For linear predictor Y < 0, age = `(1+20)·exp(Y) − 1` (exponential branch) | Yes | Source #2 |
| INV-3 | At Y = 0 exactly, age = 20.0 (adult age) | Yes | Source #2 |
| INV-4 | Age is monotonically non-decreasing in the linear predictor Y | Yes | Source #2 (both branches strictly increasing; continuous at 0) |
| INV-5 | CpGs absent from the coefficient table do not affect the result | Yes | Source #3 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Linear branch, intercept + 3 CpGs | intercept 0.695507258; (0.0127,0.5),(−0.0312,0.8),(0.0245,0.3) → Y=0.684247258 | age = 34.369192418 | Sources #2, #3 |
| M2 | Boundary Y = 0 | intercept 0, single coef 0.0 → Y=0 | age = 20.0 | Source #2 |
| M3 | Negative branch | intercept 0, coef −2.0, β 0.5 → Y=−1.0 | age = 21·e⁻¹ − 1 = 6.7254682646002895 | Source #2 |
| M4 | Non-clock CpG ignored | extra CpG not in coefficient table present in methylation map | result equals same call without that CpG | Source #3 |
| M5 | `HorvathAntiTransform` direct, x=−2.5 | published transform in isolation | 0.7237849711018749 | Source #2 |
| M6 | Null methylation map | null first arg | `ArgumentNullException` | Contract |
| M7 | Null coefficients | null coefficient table | `ArgumentNullException` | Contract |
| M8 | Empty coefficients | empty coefficient table | `ArgumentException` | Contract (empty clock undefined) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty methylation map, intercept only | no CpG contributions, intercept 1.0 | age = 21·1.0 + 20 = 41.0 | Intercept-only behaviour |
| S2 | `HorvathAntiTransform` boundary x=0 | isolates transform boundary | 20.0 | INV-3 |
| S3 | Monotonicity | age(Y2) > age(Y1) for Y2 > Y1 across the boundary | strictly increasing | INV-4 property test |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Positive `HorvathAntiTransform`, x=1.0 | linear branch isolated | 41.0 | Cross-check with M-series |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for `CalculateEpigeneticAge`. Searched `tests/Seqeron/Seqeron.Genomics.Tests/` — only sibling epigenetics fixtures exist (Bisulfite, ChromatinState, CpGDetection). No prior coverage of the age clock.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New unit |
| M2 | ❌ Missing | New unit |
| M3 | ❌ Missing | New unit |
| M4 | ❌ Missing | New unit |
| M5 | ❌ Missing | New unit |
| M6 | ❌ Missing | New unit |
| M7 | ❌ Missing | New unit |
| M8 | ❌ Missing | New unit |
| S1 | ❌ Missing | New unit |
| S2 | ❌ Missing | New unit |
| S3 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests.cs` — all cases for this unit.
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests.cs | Canonical fixture | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact value 34.369192418 |
| M2 | ✅ Covered | Boundary 20.0 |
| M3 | ✅ Covered | 6.7254682646002895 |
| M4 | ✅ Covered | Non-clock CpG ignored |
| M5 | ✅ Covered | 0.7237849711018749 |
| M6 | ✅ Covered | ArgumentNullException |
| M7 | ✅ Covered | ArgumentNullException |
| M8 | ✅ Covered | ArgumentException |
| S1 | ✅ Covered | 41.0 intercept-only |
| S2 | ✅ Covered | 20.0 transform boundary |
| S3 | ✅ Covered | Monotonicity property |
| C1 | ✅ Covered | 41.0 transform positive |

**Total in-scope cases:** 12 | **✅:** 12

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Clock coefficient set is caller-supplied (353-CpG Horvath table not bundled; framework design) | Implementation contract; §2 Notes |

---

## 7. Open Questions / Decisions

1. **Decision:** `CalculateEpigeneticAge` is implemented as a generic linear-predictor framework taking caller-supplied coefficients + intercept, per the task policy forbidding fabricated coefficient tables. The prior repo code's hard-coded 5-CpG "example coefficients" and single-branch `exp(x)−1` transform were defects and were removed.
2. **Decision:** `PredictImprintedGenes` is listed under this unit in the Registry but is unrelated to the epigenetic-age clock and its `ImprintingScore = diff/(maternal+paternal+0.01)` / `HasDMR = diff > 0.5` formulas have no retrieved authoritative basis. It is left out of scope for EPIGEN-AGE-001 (its own evidence-based validation would be a separate unit) and is neither modified nor tested here.
