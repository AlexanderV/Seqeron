# Test Specification: ONCO-VAF-001

**Test Unit ID:** ONCO-VAF-001
**Area:** Oncology
**Algorithm:** Variant Allele Frequency Analysis (empirical VAF, Wilson 95% score interval, purity/ploidy correction)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wilson E.B. (1927), via Wikipedia Binomial proportion CI | 4 (cites primary 1) | https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval | 2026-06-14 |
| 2 | GATK Mutect2 FAQ / AlleleFraction | 3 | https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2 | 2026-06-14 |
| 3 | Tarabichi et al. (2017), Subclonal Architecture | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5538405/ | 2026-06-14 |
| 4 | CNAqc, Genome Biology (2024) | 1 | https://doi.org/10.1186/s13059-024-03170-5 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Empirical VAF = altReads / totalReads (alt AD / Σ AD) — source 2.
2. Wilson 95% score interval: center = (p̂ + z²/2n)/(1+z²/n); margin = (z/(1+z²/n))·√(p̂(1−p̂)/n + z²/4n²); z = 1.96 for 95% — source 1.
3. Wilson interval stays in [0,1] (no overshoot, non-zero width at p̂=0,1) — source 1.
4. Expected clonal VAF v = m·π / (2(1−π) + π·n_tot); diploid het ⇒ π/2 (π=0.8 ⇒ 0.4) — sources 3, 4.
5. Purity/ploidy correction (inversion): adjusted m·CCF = VAF·(2(1−π) + π·n_tot)/π — source 4.

### 1.3 Documented Corner Cases

- totalReads = 0 ⇒ VAF undefined (0/0); treated as 0 (no coverage) — source 2.
- altReads > totalReads (alignment artifact, VAF > 1) ⇒ invalid input — source 2.
- p̂ = 0 ⇒ Wilson lower = 0; p̂ = 1 ⇒ Wilson upper = 1 — source 1.
- purity = 0 ⇒ correction divides by 0, undefined — sources 3, 4.

### 1.4 Known Failure Modes / Pitfalls

1. Using the Wald interval gives overshoot (<0 or >1) at extreme p̂; Wilson must be used — source 1.
2. Mutect2 `AF` is a Bayesian estimate, not alt/total; this unit computes the empirical ratio — source 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateVAF(int altReads, int totalReads)` | OncologyAnalyzer | Canonical | Empirical VAF = alt/total |
| `CalculateVAFConfidenceInterval(int altReads, int totalReads, double confidence)` | OncologyAnalyzer | Canonical | Wilson score interval |
| `AdjustVAFForPurity(double vaf, double purity, double ploidy)` | OncologyAnalyzer | Canonical (Correction) | CNAqc inversion |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 0 ≤ VAF ≤ 1 | Yes | Definition (alt ≤ total) — source 2 |
| INV-2 | lower ≤ center ≤ upper for the Wilson interval | Yes | Source 1 |
| INV-3 | 0 ≤ lower and upper ≤ 1 (no overshoot) | Yes | Source 1 |
| INV-4 | AdjustVAFForPurity(π/2, π, 2) = 1 (diploid het round-trip) | Yes | Sources 3, 4 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | VAF 25/100 | Empirical VAF | 0.25 | Source 2 |
| M2 | VAF 50/100 | Empirical VAF | 0.50 | Source 2 |
| M3 | VAF 10/10 | Full support | 1.00 | Source 2 |
| M4 | VAF 0/10 | No alt support | 0.00 | Source 2 |
| M5 | VAF totalReads=0 | No coverage | 0.0 (defined) | Source 2 |
| M6 | VAF altReads>totalReads | 11/10 artifact | ArgumentOutOfRangeException | Source 2 |
| M7 | VAF negative count | alt=-1 | ArgumentOutOfRangeException | Source 2 |
| M8 | Wilson 25/100 | 95% CI center/lo/hi | center 0.2592487019, lo 0.1754509400, hi 0.3430464637 | Source 1 |
| M9 | Wilson 50/100 | 95% CI symmetric center | center 0.5, lo 0.4038298286, hi 0.5961701714 | Source 1 |
| M10 | Wilson 5/20 | small n | center 0.2902825314, lo 0.1118600528, hi 0.4687050100 | Source 1 |
| M11 | Wilson 0/10 no overshoot | p̂=0 | lo 0.0, hi 0.2775401688 | Source 1 |
| M12 | Wilson 10/10 no overshoot | p̂=1 | lo 0.7224598312, hi 1.0 | Source 1 |
| M13 | Adjust 0.4,0.8,2 | diploid het | 1.0 | Sources 3, 4 |
| M14 | Adjust 0.2,0.5,2 | diploid | 0.8 | Source 4 |
| M15 | Adjust 0.3,0.5,4 | tetraploid | 1.8 | Source 4 |
| M16 | Adjust purity=0 | division by zero | ArgumentOutOfRangeException | Sources 3, 4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Wilson confidence out of (0,1) | confidence=0 / 1 / 1.5 | ArgumentOutOfRangeException | Validation |
| S2 | Adjust negative vaf | vaf=-0.1 | ArgumentOutOfRangeException | Validation |
| S3 | Adjust ploidy ≤ 0 | ploidy=0 | ArgumentOutOfRangeException | Validation |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV sweep | 0≤VAF≤1 and lo≤center≤hi over many inputs | all hold | Property test (INV-1..3) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for `CalculateVAF`, `CalculateVAFConfidenceInterval`, or `AdjustVAFForPurity` (these methods did not exist before this unit). Sibling unit ONCO-SOMATIC-001 has `OncologyAnalyzer_CallSomaticMutations_Tests.cs`, unrelated.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M16 | ❌ Missing | New methods, no prior tests |
| S1–S3 | ❌ Missing | New methods |
| C1 | ❌ Missing | New methods |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateVAF_Tests.cs` — all cases for the three methods.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_CalculateVAF_Tests.cs | Canonical | 20 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1–M7 | ❌ Missing | Implemented CalculateVAF tests | ✅ Done |
| 2 | M8–M12 | ❌ Missing | Implemented Wilson CI tests | ✅ Done |
| 3 | M13–M16 | ❌ Missing | Implemented AdjustVAFForPurity tests | ✅ Done |
| 4 | S1–S3 | ❌ Missing | Implemented validation tests | ✅ Done |
| 5 | C1 | ❌ Missing | Implemented property/invariant test | ✅ Done |

**Total items:** 5
**✅ Done:** 5 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M16 | ✅ | Implemented, evidence-based exact values |
| S1–S3 | ✅ | Implemented |
| C1 | ✅ | Implemented (invariant sweep) |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Wilson z = 1.96 for default 95% (source-cited verbatim) | M8–M12 |
| 2 | AdjustVAFForPurity fixes normal copy number = 2 (autosomal diploid) | M13–M16 |

---

## 7. Open Questions / Decisions

1. Mutect2 reports a Bayesian `AF`; this unit deliberately implements the empirical alt/total VAF (the standard, model-free definition). Recorded in Evidence and algorithm doc §5.3.
