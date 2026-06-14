# Test Specification: ONCO-TMB-001

**Test Unit ID:** ONCO-TMB-001
**Area:** Oncology
**Algorithm:** Tumor Mutational Burden (mutations/Mb) and TMB-high classification
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Chalmers et al. (2017), Genome Medicine 9:34 — landscape of TMB | 1 | https://doi.org/10.1186/s13073-017-0424-2 (full: https://pmc.ncbi.nlm.nih.gov/articles/PMC5395719) | 2026-06-14 |
| 2 | Marcus et al. (2021) FDA Approval Summary: Pembrolizumab for TMB-High Solid Tumors, Clin Cancer Res 27(17):4685 | 2 | https://doi.org/10.1158/1078-0432.CCR-21-0327 (full: https://pmc.ncbi.nlm.nih.gov/articles/PMC8416776/) | 2026-06-14 |
| 3 | Fancello et al. / FoCR TMB Harmonization (Merino et al. 2020) review | 1–2 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7710563/ | 2026-06-14 |

### 1.2 Key Evidence Points

1. TMB = number of somatic, coding, base-substitution and indel mutations **per megabase** of genome examined; the value is `mutations / regionMb` in mut/Mb — Chalmers 2017 (src 1).
2. The FoundationOne 315-gene panel denominator is **1.1 Mb** of coding genome — Chalmers 2017 (src 1).
3. FDA TMB-High cutoff for pembrolizumab is **TMB ≥ 10 mut/Mb** (inclusive), companion diagnostic FoundationOne CDx — Marcus 2021 (src 2).
4. Harmonized reporting unit is **mut/Mb**; TMB = somatic mutations per Mb of interrogated sequence — FoCR/Merino (src 3).

### 1.3 Documented Corner Cases

- Below ~0.5 Mb of sequenced region, panel-TMB deviation from WES TMB rises sharply (value still defined; instability is a documented limitation, not an error) — Chalmers 2017.
- The ≥10 threshold is inclusive: exactly 10.0 mut/Mb is TMB-High — Marcus 2021.

### 1.4 Known Failure Modes / Pitfalls

1. Counting unfiltered variants (germline/known drivers) overstates TMB — Chalmers 2017 filters them before counting; this unit counts the caller-supplied somatic mutations.
2. Dividing by a zero region size is undefined (division by zero) — formula has Mb in the denominator.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateTMB(int mutationCount, double targetRegionMb)` | OncologyAnalyzer | Canonical | TMB = mutationCount / targetRegionMb (Chalmers 2017) |
| `CalculateTMB(IEnumerable<SomaticCall>, double targetRegionMb)` | OncologyAnalyzer | Delegate | counts Somatic-status calls, then delegates to the scalar overload |
| `ClassifyTMB(double tmb)` | OncologyAnalyzer | Canonical | High when tmb ≥ 10 mut/Mb (FDA), else Low |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | `CalculateTMB(n, r) = n / r` for r > 0 | Yes | Chalmers 2017 (src 1) |
| INV-02 | TMB ≥ 0 for n ≥ 0, r > 0 | Yes | division of non-negatives |
| INV-03 | TMB is non-decreasing in mutation count for fixed region; non-increasing in region for fixed count | Yes | division property |
| INV-04 | `ClassifyTMB(tmb) = High ⇔ tmb ≥ 10` (inclusive boundary) | Yes | Marcus 2021 (src 2) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | TMB 315-panel worked example | 11 mutations / 1.1 Mb | 10.0 mut/Mb | Chalmers 2017 (1.1 Mb panel; TMB = mut/Mb) |
| M2 | TMB exome example | 300 mutations / 30 Mb | 10.0 mut/Mb | TMB = mut/Mb (Chalmers/FoCR) |
| M3 | TMB high count | 150 / 10 Mb | 15.0 mut/Mb | TMB = mut/Mb |
| M4 | TMB zero mutations | 0 / 10 Mb | 0.0 mut/Mb | division of 0 |
| M5 | TMB region = 0 | 5 / 0 Mb | throws ArgumentOutOfRangeException | division by zero undefined |
| M6 | ClassifyTMB below cutoff | tmb = 9.9 | Low | FDA cutoff ≥10 (9.9 < 10) |
| M7 | ClassifyTMB at cutoff | tmb = 10.0 | High | FDA cutoff inclusive (= 10) |
| M8 | ClassifyTMB above cutoff | tmb = 15.0 | High | FDA cutoff ≥10 |
| M9 | ClassifyTMB zero | tmb = 0.0 | Low | 0 < 10 |
| M10 | CalculateTMB from somatic calls | 3 Somatic + 1 Germline + 1 NotDetected / 1.0 Mb | 3.0 mut/Mb (only Somatic counted) | counts somatic mutations / Mb |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Small panel < 0.5 Mb | 2 / 0.3 Mb | ≈ 6.6667 mut/Mb (defined, no throw) | Chalmers 2017: value defined; instability documented, not an error |
| S2 | Negative mutation count | -1 / 1.0 | throws ArgumentOutOfRangeException | counts are non-negative |
| S3 | Negative / NaN region | 5 / -1.0; 5 / NaN | throws ArgumentOutOfRangeException | region must be > 0 and finite |
| S4 | Null somatic-call collection | null / 1.0 | throws ArgumentNullException | input contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV-03 monotonicity | sweep counts/regions | TMB non-decreasing in count, non-increasing in region | division property |
| C2 | INV-04 boundary sweep | tmb just below / at / above 10 | classification flips only at 10 | exact boundary |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for `CalculateTMB` / `ClassifyTMB` (methods do not yet exist). Sibling Oncology tests under `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_*_Tests.cs` provide the convention.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M10, S1–S4, C1–C2 | ❌ Missing | new unit; no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateTMB_Tests.cs` — all cases for this unit.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_CalculateTMB_Tests.cs` | Canonical | 16 |

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
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | S1 | ❌ Missing | Implemented | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented | ✅ Done |
| 13 | S3 | ❌ Missing | Implemented | ✅ Done |
| 14 | S4 | ❌ Missing | Implemented | ✅ Done |
| 15 | C1 | ❌ Missing | Implemented | ✅ Done |
| 16 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | CalculateTMB_11over1_1Mb_Returns10 |
| M2 | ✅ | CalculateTMB_300over30Mb_Returns10 |
| M3 | ✅ | CalculateTMB_150over10Mb_Returns15 |
| M4 | ✅ | CalculateTMB_ZeroMutations_ReturnsZero |
| M5 | ✅ | CalculateTMB_ZeroRegion_Throws |
| M6 | ✅ | ClassifyTMB_BelowCutoff_IsLow |
| M7 | ✅ | ClassifyTMB_AtCutoff_IsHigh |
| M8 | ✅ | ClassifyTMB_AboveCutoff_IsHigh |
| M9 | ✅ | ClassifyTMB_Zero_IsLow |
| M10 | ✅ | CalculateTMB_FromSomaticCalls_CountsOnlySomatic |
| S1 | ✅ | CalculateTMB_SmallPanel_ComputesRatioWithoutThrowing |
| S2 | ✅ | CalculateTMB_NegativeCount_Throws |
| S3 | ✅ | CalculateTMB_InvalidRegion_Throws |
| S4 | ✅ | CalculateTMB_NullCalls_Throws |
| C1 | ✅ | CalculateTMB_Monotonicity_HoldsOverSweep |
| C2 | ✅ | ClassifyTMB_BoundarySweep_FlipsOnlyAtTen |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Two-tier classification (High ⇔ tmb ≥ 10) using the only source-backed cutoff; the Registry's unsourced 6/20 Low/Intermediate/High boundaries are NOT implemented (would be fabricated). | INV-04, M6–M9, C2 |

---

## 7. Open Questions / Decisions

1. **Threshold conflict (resolved by evidence).** The Registry by-area note lists "Low (<6/Mb), Intermediate (6–20/Mb), High (>20/Mb)". No authoritative source was retrieved defining the 6 or 20 boundaries; the only harmonized, FDA-approved cutoff retrieved is **TMB-High = TMB ≥ 10 mut/Mb** (Marcus 2021; F1CDx companion diagnostic). Per the evidence-first policy, `ClassifyTMB` implements the source-backed ≥10 cutoff (Low/High), and the unsourced 6/20 boundaries are not implemented. Checklist entry updated to reflect the FDA cutoff.
