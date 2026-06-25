# Test Specification: ONCO-MRD-001

**Test Unit ID:** ONCO-MRD-001
**Area:** Oncology
**Algorithm:** Minimal (Molecular) Residual Disease Detection — tumor-informed panel-level ctDNA MRD call
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

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
| 6 | INVAR2 reference implementation (nrlab-CRUK), detectionFunctions.R / generalisedLikelihoodRatioTest.R | 3 | https://github.com/nrlab-CRUK/INVAR2 | 2026-06-23 |
| 7 | Lanczos (1964), J. SIAM Numer. Anal. 1(1):86–96 (log-gamma) | 1 | https://doi.org/10.1137/0701008 | 2026-06-23 |
| 8 | Silverman (1986), Density Estimation for Statistics and Data Analysis (Gaussian kernel eq. 2.2a, bandwidth eq. 3.31) | 1 | https://ned.ipac.caltech.edu/level5/March02/Silverman/paper.pdf | 2026-06-24 |
| 9 | R `bw.nrd0` / `density` (Silverman's rule; Gaussian kernel; `adjust`) | 3 | https://stat.ethz.ch/R-manual/R-devel/library/stats/html/bandwidth.html ; .../density.html ; source https://raw.githubusercontent.com/wch/r-source/trunk/src/library/stats/R/bandwidths.R | 2026-06-24 |

### 1.2 Key Evidence Points

1. Tumor-informed MRD tracks up to **16** patient-specific somatic SNVs selected from tumor WES — Signatera white paper.
2. A plasma sample is **MRD/ctDNA-positive when ≥ 2 of the tracked variants are detected**; < 2 ⇒ negative — PMC9265001 Table 1 (verbatim, quoting Reinert 2019).
3. Panel-level Poisson detection probability `p = 1 − e^(−n·f·m)` (n=genome equivalents, f=VAF, m=tracked mutations) — Signatera white paper Figure 2; reuses ONCO-CTDNA-001 `CtDnaDetectionProbability`.
4. ctDNA burden summarized as the **integrated (depth-weighted) mutant allele fraction (IMAF)** across tracked loci — Wan 2020.
5. **INVAR background subtraction + AF-weighted GLRT** (caller supplies per-locus background `e`): per-read mixture `q = AF·(1−e)·p + (1−AF)·e·p + e·(1−p)`; log-likelihood `Σ[lchoose(R,M)+M·log(q)+(R−M)·log(1−q)]/n`; EM for `p̂` (init 0.01, 200 iters); detection statistic `LR = logL(p̂)−logL(0)` — INVAR2 `calc_log_likelihood` / `estimate_p_EM` / `calc_likelihood_ratio`.
6. **IMAFv2** = depth-weighted mean of `max(0, locusVAF − background)` (per-context background subtraction then aggregation) — INVAR2 `calculateIMAFv2`.
7. Zero background floored to `1/depth`; only loci with tumour AF > 0 are informative — INVAR2 `doMain`.
8. **KDE-smoothed fragment-size weight (opt-in)**: tumour cfDNA is shorter, so the size profiles `P0`/`P1` up-weight short tumour-like fragments in the with-RL likelihood. INVAR2's `estimate_real_length_probability` smooths the per-length counts with R `density()` (Gaussian kernel, Silverman `bw.nrd0`, `adjust = 0.03`) and integrates the smoothed density over each integer bin `[L−0.5, L+0.5]`. `FragmentSizeProfile.FromKernelDensity` implements this: Gaussian KDE `f̂(t)=Σ wᵢ·φ((t−xᵢ)/h)/h` (Silverman eq. 2.2a), bandwidth `h = adjust·0.9·min(σ̂,IQR/1.34)·n^(−1/5)` or explicit, integrated analytically via `Φ(z)=½[1+erf(z/√2)]` and renormalised. Default `FragmentSizeProfile` constructor (discrete `COUNT/TOTAL`) unchanged.

### 1.3 Documented Corner Cases

- Exactly 1 variant detected ⇒ MRD-negative (below the ≥2 threshold) — PMC9265001/Reinert 2019.
- 0 variants detected ⇒ MRD-negative.
- Empty marker panel ⇒ undefined (no loci to interrogate) — invalid input.

### 1.4 Known Failure Modes / Pitfalls

1. Tracking ≤ 8 markers compromises sensitivity at ≤ 0.1% VAF (calling rule unaffected) — Signatera white paper.
2. Per-locus signal must exceed background; loci with no alt reads do not count as detected — Wan 2020.
3. Without per-locus background subtraction, a flat read pool counts background errors as ctDNA signal; without AF weighting, low-SNR loci dilute the statistic — INVAR2 motivation.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DetectMRD(plasmaVariants, positivityThreshold, minSupportingReads)` | OncologyAnalyzer | Canonical | Panel-level ≥2 rule; reports detected count, IMAF, Poisson p |
| `TrackVariantsOverTime(timepoints, ...)` | OncologyAnalyzer | Canonical | Longitudinal per-timepoint MRD status + first-positive timepoint |
| `EstimateInvarSignal(loci, detectionThreshold)` | OncologyAnalyzer | Canonical | INVAR background-subtracted, AF-weighted GLRT: IMAFv2, ML `p̂`, LR, detection call |
| `IntegratedMutantAlleleFractionV2(loci)` | OncologyAnalyzer | Canonical | Background-subtracted, depth-weighted aggregate (IMAFv2) |
| `EstimateInvarSignalWithSize(molecules, sizeProfile, detectionThreshold)` | OncologyAnalyzer | Canonical | INVAR fragment-size-weighted GLRT (with-RL): ML `p̂`, LR, detection call |
| `FragmentSizeProfile(mutantCounts, normalCounts, …)` (ctor) | OncologyAnalyzer.FragmentSizeProfile | Canonical (default) | Discrete empirical size profile `PROPORTION = COUNT/TOTAL` (default; unchanged) |
| `FragmentSizeProfile.FromKernelDensity(mutantCounts, normalCounts, bandwidth, bandwidthAdjust, …)` | OncologyAnalyzer.FragmentSizeProfile | Canonical (opt-in) | KDE-smoothed size profile: Gaussian KDE (Silverman bandwidth / `adjust`) integrated over each integer bin |
| `SuppressOutlierLoci(loci, α, afThreshold, maxMutantReads)` | OncologyAnalyzer | Canonical | INVAR patient-specific outlier suppression (Bonferroni binomial test) |
| `EstimateLocusBackground(controlObservations, controlProportion, maxBgAf)` | OncologyAnalyzer | Canonical | Control-derived per-locus background error + locus-noise verdict |
| `PassesBothStrandsFilter(altF, altR)` | OncologyAnalyzer | Canonical | INVAR both-strands filter |
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
| INV-6 | IMAFv2 ≥ 0; a locus with VAF ≤ background contributes 0 | Yes | INVAR2 calculateIMAFv2 (pmax(0,·)) |
| INV-7 | Pure-background sample ⇒ p̂ ≈ 0 and LR ≈ 0 (not detected) | Yes | INVAR2 GLRT (synthetic dataset inj=0) |
| INV-8 | LR monotone non-decreasing in injected ctDNA signal | Yes | INVAR2 GLRT (synthetic dataset) |
| INV-9 | AF-weighted LR ≥ flat-pooled (mean-AF) LR on same data | Yes | INVAR2 GLRT (AF-weighting dataset) |
| INV-10 | Short-tumour-fragment size-weighted LR ≥ no-size LR on same molecules | Yes | INVAR2 calc_likelihood_ratio_with_RL (size dataset) |
| INV-11 | Flat size profile (P1==P0) ⇒ size-weighted LR == no-size LR | Yes | INVAR2 with-RL (size factor cancels) |
| INV-12 | A locus whose mutant-read count exceeds the null binomial tail (≤ Bonferroni α/n) is flagged outlier | Yes | INVAR2 repolish |
| INV-13 | Estimated background = Σ(ALT_F+ALT_R)/Σ DP over controls; locus-noise pass ⟺ signalFrac < proportion AND bg < maxBg | Yes | INVAR2 createLociErrorRateTable |
| INV-14 | BothStrands pass ⟺ (ALT_F>0 AND ALT_R>0) OR no alt reads | Yes | INVAR2 BOTH_STRANDS.PASS |
| INV-15 | KDE-smoothed size masses ≥ 0, integrate to 1 over the support, and are unimodal around a single observed mode | Yes | Silverman 1986 (Gaussian kernel ∫K=1, smooth single bump) |

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
| M9 | IMAFv2 two loci above bg | (50,1000,bg0.01),(20,1000,bg0.01) | IMAFv2 = (0.04·1000+0.01·1000)/2000 = 0.025 | INVAR2 calculateIMAFv2 |
| M10 | IMAFv2 locus below bg | second locus VAF 0.005 < bg 0.01 | contributes 0 ⇒ IMAFv2 = 0.02 | INVAR2 pmax(0,·) |
| M11 | IMAFv2 pure background | every VAF == bg | IMAFv2 = 0 | INVAR2 (subtraction removes noise) |
| M12 | GLRT pure background | 50 loci, M=1/1000, e=0.001, AF=0.4 | p̂ ≈ 0 (≈3.3e-5), LR ≈ 0; not detected | INVAR2 GLRT (inj=0) |
| M13 | GLRT recovers 1% | 50 loci, M=5/1000, AF=0.4, e=0.001 | p̂ ≈ 0.01002, LR ≈ 4.06; detected | INVAR2 GLRT (inj=0.01) |
| M14 | GLRT recovers 5% | 50 loci, M=21/1000, AF=0.4 | p̂ ≈ 0.0501, LR ≈ 44.14; detected | INVAR2 GLRT (inj=0.05) |
| M15 | GLRT monotonicity | M/locus=[1,5,5,9,21] | LR non-decreasing; endpoints ≈0 and ≈44.14 | INVAR2 GLRT (dataset) |
| M16 | AF weighting boosts LR | mixture AF{0.5,0.05}, inj 0.008 | weighted LR ≈ 2.66 > flat LR ≈ 1.91 | INVAR2 GLRT (AF-weighting dataset) |
| M17 | Size weighting boosts LR | 200 mols, short tumour frags, AF 0.4, e 0.002 | size LR = 0.19691792427890276, p̂ = 0.12042621132507245; size LR > no-size LR | INVAR2 calc_likelihood_ratio_with_RL (size dataset) |
| M18 | Flat profile = no-size | same mols, P1==P0 | size-weighted LR == no-size LR | INVAR2 with-RL (cancellation) |
| M19 | Outlier flagged & removed | 9×(1/1000) + 1×(50/1000), AF 0.4, e 0.001 | clean tail 0.6323…, planted tail 3.726e-66 ≤ 0.005 ⇒ planted is outlier, clean are not | INVAR2 repolish |
| M20 | Outlier removal recovers signal | drop planted locus | kept = 9 loci; IMAFv2 = 0 (clean = background); inflated before removal | INVAR2 repolish + calculateIMAFv2 |
| M21 | Background = pooled control AF | 20 controls, 2 alt/1000 each | BACKGROUND_AF = 40/20000 = 0.002 (= injected) | INVAR2 createLociErrorRateTable |
| M22 | Locus-noise clean passes | 20 controls, signal in 1 (1 alt) | bg 5e-5, signalFrac 0.05 ⇒ LOCUS_NOISE.PASS true | INVAR2 createLociErrorRateTable |
| M23 | Recurrent signal fails | signal in 5/20 | signalFrac 0.25 ≥ 0.1 ⇒ LOCUS_NOISE.PASS false | INVAR2 createLociErrorRateTable |
| M24 | High background fails | 1 control, 250 alt/1000 | bg 0.0125 ≥ 0.01 ⇒ LOCUS_NOISE.PASS false | INVAR2 createLociErrorRateTable |
| M25 | Both-strands filter | (3,2),(3,0),(0,4),(0,0) | true, false, false, true | INVAR2 BOTH_STRANDS.PASS |
| K1 | KDE exact bin integral | 1 obs at len 100, h=0.5, support {99,100,101} | P(100)=0.684537604065696; P(101)=P(99)=0.15773119796715201 | Silverman eq.2.2a + Φ(1),Φ(3); analytic Gaussian-bin integral |
| K2 | KDE integrates to 1 | counts {120:80,170:20} & {120:20,170:80}, h=10, support 60..300 | Σ masses = 1 (mutant and normal) | Silverman eq.2.2 (∫K=1) + renormalisation |
| K3 | KDE unimodal | one dominant length 130, h=8 | strictly increasing to 130, strictly decreasing after | Silverman (smooth single Gaussian bump) |
| K4 | KDE smooths sparse bins | counts {120,170}, h=6 | unobserved 122 gets positive KDE weight > discrete uniform fall-back | INVAR2 estimate_real_length_probability (smoothing) |
| K6 | KDE drives with-RL GLRT | KDE short-tumour vs long-normal profiles | size-weighted LR > no-size LR | INVAR2 calc_likelihood_ratio_with_RL (with KDE P0/P1) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | threshold=1 | single detected variant with threshold 1 | MRD-positive | parameterized threshold |
| S2 | threshold=3, 2 detected | stricter threshold | MRD-negative | parameterized threshold |
| S3 | minSupportingReads=3 | locus with 2 alt reads not counted | detected excludes it | per-variant presence cutoff |
| S4 | all-negative timeline | det=[0,0,0] | all neg; FirstPositiveIndex=−1 | INV-5 |
| S5 | detectionThreshold gates | LR≈4.06 vs thresholds 0 and 5 | detected at 0, not detected at 5 | specificity knob |
| S6 | zero background | e=0, clear signal | LR finite, p̂>0, detected | INVAR2 1/depth floor |
| K5 | KDE auto bandwidth / adjust | two clusters, adjust 0.3 vs 1.5 | both normalised; larger adjust ⇒ more valley mass | R bw.nrd0 Silverman rule + adjust multiplier |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | null panel | DetectMRD(null) | ArgumentNullException | validation |
| C2 | empty panel | DetectMRD(empty) | ArgumentException | validation |
| C3 | invalid threshold | positivityThreshold < 1 | ArgumentOutOfRangeException | validation |
| C4 | null timepoints | TrackVariantsOverTime(null) | ArgumentNullException | validation |
| C7 | IMAFv2 no coverage | locus depth 0 | IMAFv2 = 0 | validation |
| C8 | IMAFv2 null | IntegratedMutantAlleleFractionV2(null) | ArgumentNullException | validation |
| C9 | EstimateInvarSignal null | EstimateInvarSignal(null) | ArgumentNullException | validation |
| C10 | no informative locus | all tumour AF = 0 | ArgumentException | validation |
| C11 | negative threshold | detectionThreshold < 0 | ArgumentOutOfRangeException | validation |
| C12 | tumour AF > 1 | AF = 1.5 | ArgumentOutOfRangeException | validation |
| C13 | background ≥ 1 | e = 1.0 | ArgumentOutOfRangeException | validation |
| C14 | size: null molecules | EstimateInvarSignalWithSize(null, profile) | ArgumentNullException | validation |
| C15 | size: null profile | EstimateInvarSignalWithSize(mols, null) | ArgumentNullException | validation |
| C16 | size: no informative molecule | all AF = 0 | ArgumentException | validation |
| C17 | outlier: null loci | SuppressOutlierLoci(null) | ArgumentNullException | validation |
| C18 | outlier: empty loci | SuppressOutlierLoci(empty) | ArgumentException | validation |
| C19 | outlier: invalid α | α = 0 | ArgumentOutOfRangeException | validation |
| C20 | background: null controls | EstimateLocusBackground(null) | ArgumentNullException | validation |
| C21 | background: empty controls | EstimateLocusBackground(empty) | ArgumentException | validation |
| C22 | KDE: null counts | FromKernelDensity(null, …) / (…, null) | ArgumentNullException | validation |
| C23 | KDE: invalid params | bandwidth 0 / adjust 0 / maxLength<minLength | ArgumentOutOfRangeException | validation |
| C24 | KDE: single observed length | one length, count 100 | uniform fall-back (INVAR2 length(counts)>1 guard) | validation |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Initial unit: 16 tests for `DetectMRD` / `TrackVariantsOverTime` in `OncologyAnalyzer_DetectMRD_Tests.cs`. ctDNA primitives tested in `OncologyAnalyzer_CtDnaAnalysis_Tests.cs` (ONCO-CTDNA-001), reused for the Poisson model.
- 2026-06-23 INVAR extension: no prior tests for `EstimateInvarSignal` / `IntegratedMutantAlleleFractionV2` — all new (M9–M16, S5–S6, C7–C13).
- 2026-06-23 residual closure: no prior tests for `EstimateInvarSignalWithSize` / `SuppressOutlierLoci` / `EstimateLocusBackground` / `PassesBothStrandsFilter` — all new (M17–M25, C14–C21).
- 2026-06-24 KDE smoothing: no prior tests for `FragmentSizeProfile.FromKernelDensity` — all new (K1–K6, C22–C24).

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
| M9–M16 | ❌ Missing | INVAR extension (GLRT / IMAFv2) |
| S5–S6 | ❌ Missing | INVAR extension |
| C7–C13 | ❌ Missing | INVAR extension (validation) |
| M17–M25 | ❌ Missing | residual closure (size / outlier / background / both-strands) |
| C14–C21 | ❌ Missing | residual closure (validation) |
| K1–K6 | ❌ Missing | KDE smoothing (FromKernelDensity) |
| C22–C24 | ❌ Missing | KDE smoothing (validation) |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectMRD_Tests.cs` — all cases for this unit.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_DetectMRD_Tests.cs` | canonical | 66 |

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
| 17 | M9 | ❌ Missing | implemented | ✅ Done |
| 18 | M10 | ❌ Missing | implemented | ✅ Done |
| 19 | M11 | ❌ Missing | implemented | ✅ Done |
| 20 | M12 | ❌ Missing | implemented | ✅ Done |
| 21 | M13 | ❌ Missing | implemented | ✅ Done |
| 22 | M14 | ❌ Missing | implemented | ✅ Done |
| 23 | M15 | ❌ Missing | implemented | ✅ Done |
| 24 | M16 | ❌ Missing | implemented | ✅ Done |
| 25 | S5 | ❌ Missing | implemented | ✅ Done |
| 26 | S6 | ❌ Missing | implemented | ✅ Done |
| 27 | C7 | ❌ Missing | implemented | ✅ Done |
| 28 | C8 | ❌ Missing | implemented | ✅ Done |
| 29 | C9 | ❌ Missing | implemented | ✅ Done |
| 30 | C10 | ❌ Missing | implemented | ✅ Done |
| 31 | C11 | ❌ Missing | implemented | ✅ Done |
| 32 | C12 | ❌ Missing | implemented | ✅ Done |
| 33 | C13 | ❌ Missing | implemented | ✅ Done |
| 34 | M17 | ❌ Missing | implemented | ✅ Done |
| 35 | M18 | ❌ Missing | implemented | ✅ Done |
| 36 | M19 | ❌ Missing | implemented | ✅ Done |
| 37 | M20 | ❌ Missing | implemented | ✅ Done |
| 38 | M21 | ❌ Missing | implemented | ✅ Done |
| 39 | M22 | ❌ Missing | implemented | ✅ Done |
| 40 | M23 | ❌ Missing | implemented | ✅ Done |
| 41 | M24 | ❌ Missing | implemented | ✅ Done |
| 42 | M25 | ❌ Missing | implemented | ✅ Done |
| 43 | C14 | ❌ Missing | implemented | ✅ Done |
| 44 | C15 | ❌ Missing | implemented | ✅ Done |
| 45 | C16 | ❌ Missing | implemented | ✅ Done |
| 46 | C17 | ❌ Missing | implemented | ✅ Done |
| 47 | C18 | ❌ Missing | implemented | ✅ Done |
| 48 | C19 | ❌ Missing | implemented | ✅ Done |
| 49 | C20 | ❌ Missing | implemented | ✅ Done |
| 50 | C21 | ❌ Missing | implemented | ✅ Done |
| 51 | K1 | ❌ Missing | implemented (exact Φ-derived bin integral) | ✅ Done |
| 52 | K2 | ❌ Missing | implemented | ✅ Done |
| 53 | K3 | ❌ Missing | implemented | ✅ Done |
| 54 | K4 | ❌ Missing | implemented | ✅ Done |
| 55 | K5 | ❌ Missing | implemented | ✅ Done |
| 56 | K6 | ❌ Missing | implemented | ✅ Done |
| 57 | C22 | ❌ Missing | implemented | ✅ Done |
| 58 | C23 | ❌ Missing | implemented | ✅ Done |
| 59 | C24 | ❌ Missing | implemented | ✅ Done |

**Total items:** 59
**✅ Done:** 59 | **⛔ Blocked:** 0 | **Remaining:** 0

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
| M9 | ✅ Covered | `IntegratedMutantAlleleFractionV2_TwoLociAboveBackground_DepthWeightedSubtractedMean` |
| M10 | ✅ Covered | `IntegratedMutantAlleleFractionV2_LocusBelowBackground_ContributesZero` |
| M11 | ✅ Covered | `IntegratedMutantAlleleFractionV2_PureBackground_Zero` |
| M12 | ✅ Covered | `EstimateInvarSignal_PureBackground_NotDetectedAndZeroFraction` |
| M13 | ✅ Covered | `EstimateInvarSignal_InjectedOnePercent_RecoversFractionAndDetects` |
| M14 | ✅ Covered | `EstimateInvarSignal_InjectedFivePercent_RecoversFraction` |
| M15 | ✅ Covered | `EstimateInvarSignal_RisingSignal_LikelihoodRatioMonotoneIncreasing` |
| M16 | ✅ Covered | `EstimateInvarSignal_AfWeighting_HigherLikelihoodRatioThanFlatPooling` |
| S5 | ✅ Covered | `EstimateInvarSignal_HighDetectionThreshold_WeakSignalNotDetected` |
| S6 | ✅ Covered | `EstimateInvarSignal_ZeroBackground_FiniteAndDetects` |
| C7 | ✅ Covered | `IntegratedMutantAlleleFractionV2_NoCoverage_Zero` |
| C8 | ✅ Covered | `IntegratedMutantAlleleFractionV2_Null_Throws` |
| C9 | ✅ Covered | `EstimateInvarSignal_NullLoci_Throws` |
| C10 | ✅ Covered | `EstimateInvarSignal_NoInformativeLocus_Throws` |
| C11 | ✅ Covered | `EstimateInvarSignal_NegativeThreshold_Throws` |
| C12 | ✅ Covered | `EstimateInvarSignal_TumourAfAboveOne_Throws` |
| C13 | ✅ Covered | `EstimateInvarSignal_BackgroundRateAtOne_Throws` |
| M17 | ✅ Covered | `EstimateInvarSignalWithSize_ShortTumourFragments_HigherLikelihoodRatioThanNoSize` |
| M18 | ✅ Covered | `EstimateInvarSignalWithSize_FlatProfile_NoSizeAdvantage` |
| M19 | ✅ Covered | `SuppressOutlierLoci_PlantedHighSignalLocus_FlaggedAndRemoved` |
| M20 | ✅ Covered | `SuppressOutlierLoci_RemovingOutlier_RecoversCleanSignal` |
| M21 | ✅ Covered | `EstimateLocusBackground_CleanControls_RecoversInjectedErrorRate` |
| M22 | ✅ Covered | `EstimateLocusBackground_CleanLocus_PassesNoiseFilter` |
| M23 | ✅ Covered | `EstimateLocusBackground_RecurrentSignal_FailsNoiseFilter` |
| M24 | ✅ Covered | `EstimateLocusBackground_HighBackgroundAf_FailsNoiseFilter` |
| M25 | ✅ Covered | `PassesBothStrandsFilter_StrandSupport_MatchesInvarRule` |
| C14 | ✅ Covered | `EstimateInvarSignalWithSize_NullMolecules_Throws` |
| C15 | ✅ Covered | `EstimateInvarSignalWithSize_NullProfile_Throws` |
| C16 | ✅ Covered | `EstimateInvarSignalWithSize_NoInformativeMolecule_Throws` |
| C17 | ✅ Covered | `SuppressOutlierLoci_Null_Throws` |
| C18 | ✅ Covered | `SuppressOutlierLoci_Empty_Throws` |
| C19 | ✅ Covered | `SuppressOutlierLoci_InvalidAlpha_Throws` |
| C20 | ✅ Covered | `EstimateLocusBackground_Null_Throws` |
| C21 | ✅ Covered | `EstimateLocusBackground_Empty_Throws` |
| K1 | ✅ Covered | `FromKernelDensity_SingleObservationExplicitBandwidth_MatchesAnalyticGaussianIntegral` |
| K2 | ✅ Covered | `FromKernelDensity_OverSupport_IntegratesToOne` |
| K3 | ✅ Covered | `FromKernelDensity_SingleMode_IsUnimodalAroundTheMode` |
| K4 | ✅ Covered | `FromKernelDensity_SparseHistogram_SmoothsMassOntoEmptyBins` |
| K5 | ✅ Covered | `FromKernelDensity_AutoBandwidth_AdjustControlsSmoothness` |
| K6 | ✅ Covered | `FromKernelDensity_UsedInWithSizeGlrt_ShortTumourFragmentsRaiseLikelihood` |
| C22 | ✅ Covered | `FromKernelDensity_NullCounts_Throws` |
| C23 | ✅ Covered | `FromKernelDensity_InvalidParameters_Throws` |
| C24 | ✅ Covered | `FromKernelDensity_SinglePoint_FallsBackToUniform` |

---

## 6. Assumption Register

**Total assumptions:** 1 active (1 resolved)

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Per-variant "detected" = alt reads ≥ minSupportingReads (default 1); panel-level ≥2 rule is source-exact. Threshold is a tunable parameter. | DetectMRD per-locus detection |
| 2 | RESOLVED (2026-06-24): the **default** `FragmentSizeProfile` uses the discrete empirical proportion `COUNT/TOTAL` per length bin (unchanged); the INVAR2 KDE-smoothed estimate is now **opt-in** via `FragmentSizeProfile.FromKernelDensity` — a Gaussian KDE (Silverman 1986 eq. 2.2a) with Silverman's-rule bandwidth (R `bw.nrd0`, `adjust` multiplier) integrated analytically over each integer bin. No correctness-affecting assumption remains: callers select discrete (default) or KDE (opt-in). | FragmentSizeProfile.FromKernelDensity (opt-in) |

---

## 7. Open Questions / Decisions

1. INVAR's exact per-locus GLRT/trinucleotide background model is not publicly specified verbatim in the paper; the EXACT formulas are taken from the INVAR2 reference implementation (nrlab-CRUK/INVAR2). `EstimateInvarSignal` implements the no-size GLRT (mixture model, EM `p̂`, `LR`) and IMAFv2 verbatim from that source; the per-variant presence cutoff of the legacy `DetectMRD` remains `minSupportingReads` (default 1). The panel-level ≥2 rule and Poisson/IMAF formulas are unchanged and source-exact.
2. RESIDUAL CLOSED (2026-06-23): fragment-size weighting (`EstimateInvarSignalWithSize` from `calc_likelihood_ratio_with_RL`), patient-specific outlier suppression (`SuppressOutlierLoci` from `repolish`), locus-noise filtering and control-derived background-error estimation (`EstimateLocusBackground` / `PassesBothStrandsFilter` from `createLociErrorRateTable`) are now implemented, ported verbatim from INVAR2.
3. KDE SMOOTHING ADDED (2026-06-24): the size-histogram KDE smoothing (formerly assumption #2) is now implemented opt-in as `FragmentSizeProfile.FromKernelDensity` — a Gaussian KDE (Silverman 1986 eq. 2.2a) with Silverman's-rule bandwidth (R `bw.nrd0`, `adjust` multiplier) integrated analytically over each integer bin via `Φ` and renormalised over the support, matching INVAR2's `estimate_real_length_probability` (`density()` + integrate, `adjust = 0.03`). The default discrete `COUNT/TOTAL` profile is unchanged. The remaining honest residual is that the in-pipeline per-locus background re-estimation stays an optional separate step (`EstimateLocusBackground`) rather than being wired into the GLRT automatically — matching INVAR2's own separation of the parse stage from the detection GLRT.
