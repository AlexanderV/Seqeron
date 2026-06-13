# Test Specification: META-TAXA-001

**Test Unit ID:** META-TAXA-001
**Area:** Metagenomics
**Algorithm:** Significant Taxa Detection (Wilcoxon rank-sum / Mann–Whitney U, normal approximation)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Mann & Whitney (1947), via Wikipedia "Mann–Whitney U test" | 1 (primary) / 4 (Wikipedia) | https://doi.org/10.1214/aoms/1177730491 ; https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test | 2026-06-13 |
| 2 | SciPy `scipy.stats.mannwhitneyu` documentation | 3 | https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html | 2026-06-13 |
| 3 | Xia & Sun (2017), *Genes & Diseases* | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC6128532/ | 2026-06-13 |
| 4 | Abramowitz & Stegun 7.1.26 (erf), via John D. Cook | 2 | https://www.johndcook.com/blog/python_erf/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. `U1 = R1 − n1(n1+1)/2`; `U1 + U2 = n1·n2` — Source 1.
2. Normal approximation: `m_U = n1·n2/2`, `σ_U = sqrt(n1·n2·(n1+n2+1)/12)`, `z = (U − m_U)/σ_U` — Source 1.
3. Tie correction: `σ_ties = sqrt(n1·n2·(n1+n2+1)/12 − n1·n2·Σ(t_k³−t_k)/(12·n·(n−1)))`; tied values get midranks — Source 1.
4. Complement statistic `U2 = n1·n2 − U1`; continuity correction reduces `|U − μ|` by 0.5; default on for asymptotic method — Source 2.
5. Reference output: x=[19,22,16,29,24], y=[20,11,17,12] → U1=17, U2=3, p(cc)=0.11134688653314041, p(no-cc)=0.0864107329737 — Source 2.
6. Wilcoxon rank-sum identifies statistically significant differences in microbial taxa/OTUs between two sample groups — Source 3.
7. Repository `NormalCDF`/`Erf` is A&S 7.1.26 with the listed constants; |ε|≤1.5×10⁻⁷, so p-values match exact normal CDF to ≈1×10⁻⁶ — Source 4.

### 1.3 Documented Corner Cases

- Ties → midranks + tie-corrected σ (Source 1).
- All-tied / zero-variance groups → σ → 0, z undefined → report p = 1 (degenerate case from σ formula).
- Small samples → asymptotic approximation only (Source 2).
- Multiple taxa tested independently; no multiplicity correction applied here (Source 3).

### 1.4 Known Failure Modes / Pitfalls

1. Using a normal/t approximation without tie correction overstates σ and biases p-values conservative — Source 1.
2. Omitting the continuity correction shifts p toward smaller values vs SciPy default — Source 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindSignificantTaxa(profiles, groups, pThreshold, useContinuityCorrection)` | MetagenomicsAnalyzer | Canonical | Per-taxon two-group Mann–Whitney; returns (Taxon, U, Z, PValue, Significant). |
| `MannWhitneyU(group1, group2, useContinuityCorrection)` | MetagenomicsAnalyzer | Canonical | Core statistic + asymptotic two-tailed p-value; tie-corrected σ. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `U1 + U2 = n1·n2` | Yes | Source 1 |
| INV-2 | `0 ≤ U ≤ n1·n2` | Yes | Source 1 (consequence of INV-1, both ≥0) |
| INV-3 | p-value ∈ [0, 1] | Yes | Two-tailed `2·SF(|z|)` clamped |
| INV-4 | Significant ⇔ `PValue < pThreshold` | Yes | Method contract / Source 3 |
| INV-5 | Identical groups (all tied) → p = 1, not significant | Yes | σ→0 degenerate case |
| INV-6 | Test is symmetric: swapping group1/group2 leaves p unchanged | Yes | Source 1 (z uses max-side U / |U−μ|) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | SciPy U values | MannWhitneyU on x=[19,22,16,29,24], y=[20,11,17,12] | U1=17, U2=3 | Source 2 |
| M2 | SciPy σ and z | same data, no continuity | σ=sqrt(200/12)=4.08248290463863; z=1.7146428199482247 | Source 1+2 |
| M3 | SciPy p (no cc) | same data, continuity off | p≈0.0864107329737 (within 1e-6) | Source 2 |
| M4 | SciPy p (cc) | same data, continuity on (default) | z=1.5921683328090657; p≈0.11134688653314041 (within 1e-6) | Source 2 |
| M5 | Tortoise/hare U | tortoise ranks {12,6,5,4,3,2}-style data vs hare | U_T=11, U_H=25 | Source 1 |
| M6 | Complement invariant | U1+U2 = n1·n2 for arbitrary data | equality exact | Source 1 (INV-1) |
| M7 | Ties / midranks | data with tied values uses tie-corrected σ | σ smaller than no-tie formula; matches midrank derivation | Source 1 |
| M8 | All-tied groups | both groups identical constant | p = 1.0, not significant | INV-5 |
| M9 | FindSignificantTaxa flags | two groups, one taxon clearly separated, one overlapping | separated taxon Significant=true, overlapping false; p<threshold ⇔ significant | Source 3 / INV-4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Symmetry | swap group1/group2 in MannWhitneyU | p identical | INV-6 |
| S2 | p in [0,1] | random-ish separated data | 0 ≤ p ≤ 1 | INV-3 |
| S3 | FindSignificantTaxa zero abundance fill | taxon absent in some profiles | absent → 0 abundance, still tested | ASSUMPTION-3 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Continuity toggle direction | cc p vs no-cc p, same data | p(cc) > p(no-cc) | Source 2 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for `FindSignificantTaxa` or `MannWhitneyU`. No existing `FindSignificantTaxa` method in `MetagenomicsAnalyzer` (a `DifferentialAbundance` method exists using a Welch-t normal approximation — out of scope for this unit, untouched).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing → ✅ | New |
| M2 | ❌ Missing → ✅ | New |
| M3 | ❌ Missing → ✅ | New |
| M4 | ❌ Missing → ✅ | New |
| M5 | ❌ Missing → ✅ | New |
| M6 | ❌ Missing → ✅ | New |
| M7 | ❌ Missing → ✅ | New |
| M8 | ❌ Missing → ✅ | New |
| M9 | ❌ Missing → ✅ | New |
| S1 | ❌ Missing → ✅ | New |
| S2 | ❌ Missing → ✅ | New |
| S3 | ❌ Missing → ✅ | New |
| C1 | ❌ Missing → ✅ | New |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_FindSignificantTaxa_Tests.cs` — all cases for `MannWhitneyU` and `FindSignificantTaxa`.
- **Remove:** nothing (no pre-existing tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| MetagenomicsAnalyzer_FindSignificantTaxa_Tests.cs | Canonical unit tests | 13 + null/empty edge tests |

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
| 10 | S1 | ❌ Missing | Implemented | ✅ Done |
| 11 | S2 | ❌ Missing | Implemented | ✅ Done |
| 12 | S3 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | MannWhitneyU_ScipyExample_ReturnsU17And3 |
| M2 | ✅ Covered | MannWhitneyU_ScipyExample_NoContinuity_MatchesSigmaAndZ |
| M3 | ✅ Covered | MannWhitneyU_ScipyExample_NoContinuity_MatchesPValue |
| M4 | ✅ Covered | MannWhitneyU_ScipyExample_WithContinuity_MatchesZAndPValue |
| M5 | ✅ Covered | MannWhitneyU_TortoiseHare_ReturnsU11And25 |
| M6 | ✅ Covered | MannWhitneyU_AnyInput_USumEqualsProductOfSizes |
| M7 | ✅ Covered | MannWhitneyU_WithTies_UsesTieCorrectedSigma |
| M8 | ✅ Covered | MannWhitneyU_IdenticalGroups_PValueIsOne |
| M9 | ✅ Covered | FindSignificantTaxa_SeparatedAndOverlapping_FlagsCorrectly |
| S1 | ✅ Covered | MannWhitneyU_SwapGroups_PValueUnchanged |
| S2 | ✅ Covered | MannWhitneyU_AnyInput_PValueInUnitInterval |
| S3 | ✅ Covered | FindSignificantTaxa_AbsentTaxon_TreatedAsZeroAbundance |
| C1 | ✅ Covered | MannWhitneyU_ContinuityCorrection_IncreasesPValue |

In-scope cases: 13. ✅ = 13. (Plus null/empty/single-group validation tests.)

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Continuity correction default-on (SciPy default `use_continuity=True`) | M4, C1 |
| 2 | Two-tailed alternative (`p = 2·SF(|z|)`) | M3, M4 |
| 3 | Taxon absent in a profile → abundance 0 in that group's vector | M9, S3 |

All three are reference-implementation/standard-practice defaults (not invented constants); each is exposed/observable and source-anchored.

---

## 7. Open Questions / Decisions

1. Exact (permutation) p-values for very small n are not implemented; the asymptotic normal approximation is used per Source 1/2. Documented as a limitation, not a defect.
2. Multiplicity correction across taxa is the caller's responsibility (Source 3); not part of this unit.
