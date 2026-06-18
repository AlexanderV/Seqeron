# Test Specification: ONCO-SOMATIC-001

**Test Unit ID:** ONCO-SOMATIC-001
**Area:** Oncology
**Algorithm:** Somatic Mutation Calling (tumor vs matched normal classification)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Saunders et al. (2012), Strelka, Bioinformatics 28(14):1811–1817 | 1 | https://doi.org/10.1093/bioinformatics/bts271 | 2026-06-14 |
| 2 | Kim et al. (2018), Strelka2, Nature Methods 15:591–594 | 1 | https://doi.org/10.1038/s41592-018-0051-x | 2026-06-14 |
| 3 | Benjamin et al. / Broad (2019), GATK Mutect2 mutect.tex | 3 | https://raw.githubusercontent.com/broadinstitute/gatk/master/docs/mutect/mutect.tex | 2026-06-14 |
| 4 | Yan et al. (2021), Sci. Rep. 11:11640 | 1 | https://doi.org/10.1038/s41598-021-91142-1 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Somatic state S = {(f_t, f_n): f_t ≠ f_n}, restricted to homozygous-reference normal genotype P(S, G_n=ref/ref|D) — Saunders et al. (2012).
2. A variant is somatic when present in the tumor and absent from the matched normal; given a matched normal, Mutect2 skips variants clearly present in the germline; "If we have no matched normal, ℓ_n = 1" — Benjamin et al. (2019).
3. Tumor limit of detection: WES VAF LoD = 5%; calls ≤ 5% VAF are frequently sequencing errors — Yan et al. (2021).
4. VAF = altReads / totalReads, continuous allele frequency compared between tumor and normal — Kim et al. (2018).

### 1.3 Documented Corner Cases

- Tumor-only mode (no matched normal) — Benjamin et al. (2019): ℓ_n = 1.
- Low tumor purity / impure samples — Benjamin et al. (2019): allele fraction deviates from diploid 1/2.
- Sub-5% VAF subclonal calls — Yan et al. (2021): frequently errors → not detected at standard LoD.
- LOH / copy-number regions over-called by raw somatic probability — Saunders et al. (2012).

### 1.4 Known Failure Modes / Pitfalls

1. Clonal-hematopoiesis (CHIP) contamination of the normal raises normal VAF → classified germline — Benjamin et al. (2019).
2. Reporting sub-5% VAF variants → false positives — Yan et al. (2021).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CallSomaticMutations(variants, tumorVafThreshold, normalVafThreshold)` | OncologyAnalyzer | Canonical | Core tumor/normal classification |
| `Classify(variant, ...)` | OncologyAnalyzer | Canonical | Single-variant form of the above |
| `FilterGermlineVariants(variants, ...)` | OncologyAnalyzer | Canonical | Returns somatic subset only |
| `CalculateSomaticScore(variant)` | OncologyAnalyzer | Canonical | Monotone separation score in [0,1] |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Somatic ⇔ f_t ≥ tumorThreshold ∧ f_n ≤ normalThreshold | Yes | Saunders 2012; Yan 2021 |
| INV-2 | f_t < tumorThreshold ⇒ NotDetected (regardless of normal) | Yes | Yan 2021 (5% LoD) |
| INV-3 | 0 ≤ SomaticScore ≤ 1; score = max(0, f_t − f_n) | Yes | derived (ASSUMPTION) |
| INV-4 | FilterGermlineVariants output = exactly the Somatic-status subset | Yes | Mutect2 (Benjamin 2019) |
| INV-5 | totalReads = 0 ⇒ VAF = 0 (allele absent at uncovered site) | Yes | VAF def. (Kim 2018) |
| INV-6 | Output count = input count; order preserved | Yes | implementation contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Clear somatic | tumor 25/100, normal 0/100 | Somatic, score 0.25 | Saunders 2012 (f_t≠f_n, ref/ref) |
| M2 | Germline het | tumor 48/100, normal 50/100 | Germline, score 0 | Benjamin 2019 (germline filter) |
| M3 | Sub-LoD tumor | tumor 2/100, normal 0/100 | NotDetected, score 0 | Yan 2021 (5% LoD) |
| M4 | Tumor-only mode | tumor 20/100, normal 0/0 | Somatic, score 0.20 | Benjamin 2019 (ℓ_n=1) |
| M5 | Tumor threshold boundary | tumor 5/100 (=0.05), normal 0/100 | Somatic | DefaultTumorVafThreshold=0.05 |
| M6 | Normal threshold boundary | tumor 30/100, normal 1/100 (=0.01) | Somatic, score 0.29 | DefaultNormalVafThreshold=0.01 |
| M7 | CHIP-like normal | tumor 30/100, normal 3/100 (=0.03) | Germline | Benjamin 2019 (present in normal) |
| M8 | FilterGermlineVariants | mixed panel | returns only the Somatic calls | Benjamin 2019 |
| M9 | CalculateSomaticScore | tumor 25/100, normal 5/100 | 0.20 (=0.25−0.05) | INV-3 (separation) |
| M10 | Order & count preserved | 3-variant panel | 3 calls in input order | INV-6 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Just above normal threshold | tumor 30/100, normal 2/100 (=0.02) | Germline | f_n > 0.01 |
| S2 | Custom thresholds reclassify | tumor 3/100, normal 0, tumorThreshold=0.02 | Somatic | parameter effect |
| S3 | Score 0 when normal ≥ tumor | tumor 10/100, normal 20/100 | score 0 | max(0, f_t−f_n) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty input | no variants | empty result | guard |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior tests existed for OncologyAnalyzer somatic calling. The class `OncologyAnalyzer` and the canonical test file are created new in this unit. (The Oncology project previously contained only `ImmuneAnalyzer`.)

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M10 | ❌ Missing | brand-new unit, no prior tests |
| S1–S3 | ❌ Missing | brand-new unit |
| C1, null/empty guards | ❌ Missing | brand-new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CallSomaticMutations_Tests.cs` — all cases for this unit.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_CallSomaticMutations_Tests.cs | Canonical unit fixture | 18 |

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
| 14 | C1 | ❌ Missing | Implemented | ✅ Done |
| 15 | Guards (null, neg reads, alt>total, bad threshold) | ❌ Missing | Implemented | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M10 | ✅ Covered | exact evidence-based values |
| S1–S3 | ✅ Covered | edge cases |
| C1 + guards | ✅ Covered | null/empty/invalid input |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Normal "absent" ceiling = 1% VAF (configurable parameter) | INV-1, M6, M7, S1 |
| 2 | Somatic score = max(0, f_t − f_n) | INV-3, M9, S3 |

---

## 7. Open Questions / Decisions

1. Full Bayesian caller LOD/somatic-quality models (Strelka somaticLOD, Mutect2 TLOD) are out of scope; the rule-based tumor/normal classification per Saunders 2012 / Yan 2021 is implemented and tested. The score is a documented surrogate, not a caller LOD.
