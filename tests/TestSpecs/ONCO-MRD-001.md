# Test Specification: ONCO-MRD-001

**Test Unit ID:** ONCO-MRD-001
**Area:** Oncology
**Algorithm:** Minimal (Molecular) Residual Disease Detection — tumor-informed panel-level ctDNA MRD call
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Reinert et al. (2019), JAMA Oncol 5(8):1124–1131 | 1 | https://pubmed.ncbi.nlm.nih.gov/31070691/ (DOI:10.1001/jamaoncol.2019.0528) | 2026-06-15 |
| 2 | Natera Signatera analytical-validation white paper (2020) | 3 | https://www.natera.com/wp-content/uploads/2020/11/Oncology-Clinical-A-personalized-tumor-informed-approach-to-detect-molecular-residual-disease-SGN_SR_WP.pdf | 2026-06-15 |
| 3 | Wan et al. (2020), Sci Transl Med 12(548):eaaz8084 (INVAR/IMAF) | 1 | https://www.science.org/doi/10.1126/scitranslmed.aaz8084 | 2026-06-15 |
| 4 | Tumor-informed ctDNA MRD review, Table 1 (quotes Reinert/Signatera rule) | 4 | https://pmc.ncbi.nlm.nih.gov/articles/PMC9265001/ | 2026-06-15 |
| 5 | Avanzini et al. (2020), Sci Adv 6(50):eabc4308 (Poisson p=1−e^(−λ)) | 1 | DOI:10.1126/sciadv.abc4308 | 2026-06-15 |

### 1.2 Key Evidence Points

1. Tumor-informed MRD tracks up to **16** patient-specific somatic SNVs selected from tumor WES — Signatera white paper.
2. A plasma sample is **MRD/ctDNA-positive when ≥ 2 of the tracked variants are detected**; < 2 ⇒ negative — PMC9265001 Table 1 (verbatim, quoting Reinert 2019).
3. Panel-level Poisson detection probability `p = 1 − e^(−n·f·m)` (n=genome equivalents, f=VAF, m=tracked mutations) — Signatera white paper Figure 2; reuses ONCO-CTDNA-001 `CtDnaDetectionProbability`.
4. ctDNA burden summarized as the **integrated (depth-weighted) mutant allele fraction (IMAF)** across tracked loci — Wan 2020.

### 1.3 Documented Corner Cases

- Exactly 1 variant detected ⇒ MRD-negative (below the ≥2 threshold) — PMC9265001/Reinert 2019.
- 0 variants detected ⇒ MRD-negative.
- Empty marker panel ⇒ undefined (no loci to interrogate) — invalid input.

### 1.4 Known Failure Modes / Pitfalls

1. Tracking ≤ 8 markers compromises sensitivity at ≤ 0.1% VAF (calling rule unaffected) — Signatera white paper.
2. Per-locus signal must exceed background; loci with no alt reads do not count as detected — Wan 2020.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DetectMRD(plasmaVariants, positivityThreshold, minSupportingReads)` | OncologyAnalyzer | Canonical | Panel-level ≥2 rule; reports detected count, IMAF, Poisson p |
| `TrackVariantsOverTime(timepoints, ...)` | OncologyAnalyzer | Canonical | Longitudinal per-timepoint MRD status + first-positive timepoint |
| `IsVariantDetected(variant, minSupportingReads)` | OncologyAnalyzer | Internal | Per-locus presence (alt reads ≥ min) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | MRD-positive ⟺ DetectedVariantCount ≥ positivityThreshold (default 2) | Yes | PMC9265001 Table 1 / Reinert 2019 |
| INV-2 | 0 ≤ DetectedVariantCount ≤ TrackedVariantCount | Yes | counting |
| INV-3 | IMAF ∈ [0, 1] and equals Σalt / Σtotal over tracked loci | Yes | Wan 2020 |
| INV-4 | Panel Poisson p = 1 − e^(−n·f·m) ∈ [0, 1], non-decreasing in m | Yes | white paper Fig 2 |
| INV-5 | TrackVariantsOverTime preserves timepoint order; FirstPositiveIndex = earliest positive (or −1) | Yes | longitudinal monitoring (Reinert 2019) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | 2 of 3 detected | panel of 3 markers, 2 with alt reads | MRD-positive, DetectedVariantCount=2 | PMC9265001 Table 1 |
| M2 | 1 of 3 detected | only 1 marker has alt reads | MRD-negative, DetectedVariantCount=1 | PMC9265001 Table 1 |
| M3 | 0 of 3 detected | no marker has alt reads | MRD-negative, DetectedVariantCount=0 | PMC9265001 Table 1 |
| M4 | 3 of 3 detected | all markers have alt reads | MRD-positive, DetectedVariantCount=3 | PMC9265001 Table 1 |
| M5 | TrackedVariantCount | panel of 16 markers reported | TrackedVariantCount=16 | Signatera white paper |
| M6 | IMAF worked example | loci (3,200),(1,150),(0,180) | IMAF = 4/530 = 0.0075471698…; detected=2; positive | Wan 2020 (IMAF) |
| M7 | Panel Poisson p, m=16 | n=1000, f=0.001, m=16 ⇒ λ=16 | p = 1−e^(−16) = 0.9999998874648253 | white paper Fig 2 |
| M8 | Longitudinal status | 4 timepoints det=[0,1,2,3] | status=[neg,neg,pos,pos]; FirstPositiveIndex=2 | Reinert 2019 (serial) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | threshold=1 | single detected variant with threshold 1 | MRD-positive | parameterized threshold |
| S2 | threshold=3, 2 detected | stricter threshold | MRD-negative | parameterized threshold |
| S3 | minSupportingReads=3 | locus with 2 alt reads not counted | detected excludes it | per-variant presence cutoff |
| S4 | all-negative timeline | det=[0,0,0] | all neg; FirstPositiveIndex=−1 | INV-5 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | null panel | DetectMRD(null) | ArgumentNullException | validation |
| C2 | empty panel | DetectMRD(empty) | ArgumentException | validation |
| C3 | invalid threshold | positivityThreshold < 1 | ArgumentOutOfRangeException | validation |
| C4 | null timepoints | TrackVariantsOverTime(null) | ArgumentNullException | validation |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for `DetectMRD` / `TrackVariantsOverTime`. ctDNA primitives tested in `OncologyAnalyzer_CtDnaAnalysis_Tests.cs` (ONCO-CTDNA-001), reused here for the Poisson model.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| M8 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |
| C3 | ❌ Missing | new unit |
| C4 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectMRD_Tests.cs` — all cases for this unit.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_DetectMRD_Tests.cs` | canonical | 16 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | ❌ Missing | implemented | ✅ Done |
| 9 | S1 | ❌ Missing | implemented | ✅ Done |
| 10 | S2 | ❌ Missing | implemented | ✅ Done |
| 11 | S3 | ❌ Missing | implemented | ✅ Done |
| 12 | S4 | ❌ Missing | implemented | ✅ Done |
| 13 | C1 | ❌ Missing | implemented | ✅ Done |
| 14 | C2 | ❌ Missing | implemented | ✅ Done |
| 15 | C3 | ❌ Missing | implemented | ✅ Done |
| 16 | C4 | ❌ Missing | implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `DetectMRD_TwoOfThreeDetected_PositiveCall` |
| M2 | ✅ Covered | `DetectMRD_OneOfThreeDetected_NegativeCall` |
| M3 | ✅ Covered | `DetectMRD_ZeroDetected_NegativeCall` |
| M4 | ✅ Covered | `DetectMRD_ThreeOfThreeDetected_PositiveCall` |
| M5 | ✅ Covered | `DetectMRD_SixteenMarkers_ReportsTrackedCount` |
| M6 | ✅ Covered | `DetectMRD_ImafWorkedExample_DepthWeightedMeanVaf` |
| M7 | ✅ Covered | `DetectMRD_PanelDetectionProbability_PoissonM16` |
| M8 | ✅ Covered | `TrackVariantsOverTime_RisingSignal_FirstPositiveAtIndexTwo` |
| S1 | ✅ Covered | `DetectMRD_ThresholdOne_SingleDetectedIsPositive` |
| S2 | ✅ Covered | `DetectMRD_ThresholdThree_TwoDetectedIsNegative` |
| S3 | ✅ Covered | `DetectMRD_MinSupportingReadsThree_LowSupportNotDetected` |
| S4 | ✅ Covered | `TrackVariantsOverTime_AllNegative_FirstPositiveMinusOne` |
| C1 | ✅ Covered | `DetectMRD_NullPanel_Throws` |
| C2 | ✅ Covered | `DetectMRD_EmptyPanel_Throws` |
| C3 | ✅ Covered | `DetectMRD_InvalidThreshold_Throws` |
| C4 | ✅ Covered | `TrackVariantsOverTime_NullTimepoints_Throws` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Per-variant "detected" = alt reads ≥ minSupportingReads (default 1); panel-level ≥2 rule is source-exact. Threshold is a tunable parameter. | DetectMRD per-locus detection |

---

## 7. Open Questions / Decisions

1. INVAR's exact per-locus GLRT/trinucleotide background model is not publicly specified verbatim; the per-variant presence cutoff is therefore exposed as `minSupportingReads` (default 1, minimal source-consistent presence). The panel-level positivity rule (≥2 detected) and the Poisson/IMAF formulas are source-exact.
