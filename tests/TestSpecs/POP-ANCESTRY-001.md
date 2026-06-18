# Test Specification: POP-ANCESTRY-001

**Test Unit ID:** POP-ANCESTRY-001
**Area:** PopGen
**Algorithm:** Ancestry Estimation (supervised / projection ADMIXTURE — EM with fixed reference allele frequencies)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Alexander, Novembre & Lange (2009), *Genome Research* 19(9):1655–1664 (ADMIXTURE) | 1 | https://doi.org/10.1101/gr.094052.109 (open: https://faculty.eeb.ucla.edu/Novembre/AlexanderEtAl_GR_2009.pdf) | 2026-06-13 |
| 2 | Alexander & Lange (2011), *BMC Bioinformatics* 12:246 (supervised mode) | 1 | https://doi.org/10.1186/1471-2105-12-246 | 2026-06-13 |
| 3 | ADMIXTURE 1.4 Software Manual, §2.10 / §2.14 | 3 | https://dalexander.github.io/admixture/admixture-manual.pdf | 2026-06-13 |

### 1.2 Key Evidence Points

1. Genotype g_ij ∈ {0,1,2} = number of copies of allele 1 — Alexander et al. (2009), Methods.
2. Log-likelihood (Eq. 2): `L = Σ_i Σ_j { g_ij ln(Σ_k q_ik f_kj) + (2−g_ij) ln(Σ_k q_ik (1−f_kj)) }` — Alexander et al. (2009).
3. EM update for ancestry (Eq. 4): `q_ik^{n+1} = 1/(2J) Σ_j [ g_ij a_ijk + (2−g_ij) b_ijk ]`, `a_ijk = q_ik f_kj / Σ_m q_im f_mj`, `b_ijk = q_ik(1−f_kj)/Σ_m q_im(1−f_mj)` — Alexander et al. (2009).
4. Constraints: `Σ_k q_ik = 1`, `q_ik ≥ 0`, `0 ≤ f_kj ≤ 1` — Alexander et al. (2009).
5. Convergence rule (Eq. 5): stop when ΔL < ε; ADMIXTURE default ε = 10⁻⁴ — Alexander et al. (2009).
6. Supervised / projection mode fixes the reference allele frequencies F and estimates only Q — ADMIXTURE manual §2.10, §2.14; Alexander & Lange (2011).

### 1.3 Documented Corner Cases

- Label permutation: Eq. 2 has K! equivalent maxima (label-invariant); pinned by labelled reference panels in this unit — Alexander et al. (2009).
- EM is a monotone-ascent algorithm (slow convergence noted) — Alexander et al. (2009).

### 1.4 Known Failure Modes / Pitfalls

1. Using a per-population genotype likelihood instead of the admixed mixture `Σ_m q_im f_mj` is NOT the ADMIXTURE model and gives wrong proportions — derived from Eq. 4 vs. Eq. 2 (Evidence §Test Datasets).
2. Reference panels must be on the same SNP set / order as the genotype vector (alignment by index) — ADMIXTURE manual §2.14 (same `.bim`).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `EstimateAncestry(individuals, referencePops, maxIterations)` | PopulationGeneticsAnalyzer | Canonical | Supervised/projection EM (Eq. 4), F fixed |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Each returned individual's proportions sum to 1 (Σ_k q_ik = 1) | Yes | Alexander et al. (2009), constraint Σ_k q_ik = 1 |
| INV-2 | All proportions are in [0,1] (q_ik ≥ 0 and ≤ 1) | Yes | Alexander et al. (2009), constraint q_ik ≥ 0 with Σ=1 |
| INV-3 | Log-likelihood (Eq. 2) is non-decreasing across EM iterations | Yes | Alexander et al. (2009), EM ascent / Eq. 5 |
| INV-4 | One result per valid individual; result keyed by reference population id | Yes | API contract / supervised mode (manual §2.10) |
| INV-5 | A uniform individual under identical reference panels stays uniform (uninformative fixed point) | Yes | Eq. 4 (Evidence Test Dataset 3) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | One iteration, symmetric 2-pop panel | f=[[0.8,0.2],[0.2,0.8]], g=[2,0], q⁰ uniform, maxIterations=1 | q_A=0.8, q_B=0.2 (exact) | Alexander et al. (2009) Eq. 4 (Test Dataset 1) |
| M2 | Single SNP g=2 | f=[[0.9],[0.1]], g=[2], 1 iteration | q_A=0.9, q_B=0.1 | Eq. 4 closed form (Test Dataset 2) |
| M3 | Single SNP g=0 | f=[[0.9],[0.1]], g=[0], 1 iteration | q_A=0.1, q_B=0.9 | Eq. 4 closed form (Test Dataset 2) |
| M4 | Proportions sum to 1 | any valid input | Σ_k q_ik = 1 within 1e-10 | INV-1 (constraint Σ_k q_ik = 1) |
| M5 | Convergence to source | f=[[0.8,0.2],[0.2,0.8]], g=[2,0], maxIterations=100 | q_A→1.0, q_B→0.0 | Eq. 4 ascent (Test Dataset 1) |
| M6 | Identical panels stay uniform | f_A=f_B, uniform individual | q=(0.5,0.5) | Eq. 4 (Test Dataset 3), INV-5 |
| M7 | Heterozygote symmetric stays uniform | f=[[0.9],[0.1]], g=[1] | q=(0.5,0.5) | Eq. 4 (Test Dataset 2) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Log-likelihood non-decreasing | record L at each iteration on a 3-pop case | L monotone non-decreasing | INV-3, property-based (O(n·k) but documents the EM ascent invariant) |
| S2 | Empty individuals | no individuals supplied | empty result | validation behavior |
| S3 | Empty reference panels | no reference pops | empty result | validation behavior |
| S4 | Mismatched genotype length | individual genotype length ≠ J | that individual skipped | validation behavior (Evidence Assumption 2 baseline) |
| S5 | Three-population panel result keys | refs A,B,C | result has exactly keys {A,B,C}, sums to 1 | INV-4 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Reference-panel order independence | permute ref panels | proportions permute consistently | Eq. 2 label invariance |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Three pre-existing tests for `EstimateAncestry` lived in `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzerTests.cs` (`EstimateAncestry_SinglePopulation_Returns100Percent`, `EstimateAncestry_TwoPopulations_SumsToOne`, `EstimateAncestry_EmptyInput_ReturnsEmpty`). All were ⚠ Weak: permissive `.Within(0.01)`, no exact evidence values, and they validated the previous (non-conforming) implementation.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 one-iteration Eq. 4 | ❌ Missing | no exact-value test of the EM update existed |
| M2 single SNP g=2 | ❌ Missing | |
| M3 single SNP g=0 | ❌ Missing | |
| M4 sum to 1 | ⚠ Weak | old test used `.Within(0.01)` on wrong impl |
| M5 convergence to source | ⚠ Weak | old "100%" test used `.Within(0.01)` on f=1.0 panel |
| M6 identical panels uniform | ❌ Missing | |
| M7 symmetric heterozygote | ❌ Missing | |
| S1 log-likelihood ascent | ❌ Missing | |
| S2 empty individuals | ⚠ Weak | old empty test conflated empty-inds and empty-refs |
| S3 empty reference panels | ❌ Missing | |
| S4 mismatched genotype length | ❌ Missing | |
| S5 three-panel keys | ❌ Missing | |
| S6 missing genotype skipped | ❌ Missing | |
| S7 all-missing uniform prior | ❌ Missing | |
| C1 panel-order invariance | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_EstimateAncestry_Tests.cs` — all evidence-based ancestry tests.
- **Remove:** the three ⚠ Weak tests and their `#region Ancestry Analysis Tests` from `PopulationGeneticsAnalyzerTests.cs` (replaced by a pointer comment).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PopulationGeneticsAnalyzer_EstimateAncestry_Tests.cs` | Canonical unit file | 15 |
| `PopulationGeneticsAnalyzerTests.cs` | Ancestry region removed (pointer comment) | 0 (for ancestry) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ | Implemented exact Eq. 4 one-iteration test | ✅ Done |
| 2 | M2 | ❌ | Implemented single-SNP g=2 exact test | ✅ Done |
| 3 | M3 | ❌ | Implemented single-SNP g=0 exact test | ✅ Done |
| 4 | M4 | ⚠ | Rewrote: exact `Sigma=1` with `.Within(1e-10)` | ✅ Done |
| 5 | M5 | ⚠ | Rewrote: convergence to (1,0) | ✅ Done |
| 6 | M6 | ❌ | Implemented identical-panels fixed point | ✅ Done |
| 7 | M7 | ❌ | Implemented symmetric heterozygote | ✅ Done |
| 8 | S1 | ❌ | Implemented log-likelihood ascent property test | ✅ Done |
| 9 | S2 | ⚠ | Rewrote: empty individuals only | ✅ Done |
| 10 | S3 | ❌ | Implemented empty reference panels | ✅ Done |
| 11 | S4 | ❌ | Implemented length-mismatch skip | ✅ Done |
| 12 | S5 | ❌ | Implemented three-panel key/sum test | ✅ Done |
| 13 | S6 | ❌ | Implemented missing-genotype skip | ✅ Done |
| 14 | S7 | ❌ | Implemented all-missing uniform prior | ✅ Done |
| 15 | C1 | ❌ | Implemented panel-order invariance | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | `EstimateAncestry_OneIterationSymmetricPanel_MatchesEq4` |
| M2 | ✅ | `EstimateAncestry_SingleSnpHomozygousAllele1_MatchesEq4` |
| M3 | ✅ | `EstimateAncestry_SingleSnpHomozygousAllele2_MatchesEq4` |
| M4 | ✅ | `EstimateAncestry_ArbitraryInput_ProportionsSumToOne` |
| M5 | ✅ | `EstimateAncestry_DiagnosticIndividual_ConvergesToSource` |
| M6 | ✅ | `EstimateAncestry_IdenticalPanels_StaysUniform` |
| M7 | ✅ | `EstimateAncestry_SymmetricHeterozygote_StaysUniform` |
| S1 | ✅ | `EstimateAncestry_LogLikelihood_IsNonDecreasingAcrossIterations` |
| S2 | ✅ | `EstimateAncestry_NoIndividuals_ReturnsEmpty` |
| S3 | ✅ | `EstimateAncestry_NoReferencePanels_ReturnsEmpty` |
| S4 | ✅ | `EstimateAncestry_MismatchedGenotypeLength_SkipsIndividual` |
| S5 | ✅ | `EstimateAncestry_ThreePanels_KeysMatchAndSumToOne` |
| S6 | ✅ | `EstimateAncestry_MissingGenotype_SkipsThatSnp` |
| S7 | ✅ | `EstimateAncestry_AllGenotypesMissing_ReturnsUniformPrior` |
| C1 | ✅ | `EstimateAncestry_PermutedPanels_PermutesProportions` |

In-scope cases: 15. ✅ count: 15. No ❌/⚠ remaining.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | EM runs up to `maxIterations` with early stop at ΔL < 10⁻⁴ (Eq. 5); budget is API shape, not a correctness constant | M5, S1 |
| 2 | Genotype value outside {0,1,2} treated as missing; individual whose genotype length ≠ J is skipped | S4 |

---

## 7. Open Questions / Decisions

1. The pre-existing implementation used a per-population genotype-likelihood responsibility update, which is NOT ADMIXTURE Eq. 4 (it produces 0.941 instead of 0.8 after one iteration on M1). Decision: correct the implementation to Eq. 4 (Present-but-nonconforming → fixed).
