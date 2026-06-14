# Test Specification: ONCO-CTDNA-001

**Test Unit ID:** ONCO-CTDNA-001
**Area:** Oncology
**Algorithm:** ctDNA Analysis (Poisson limit-of-detection, tumor-fraction estimation, mean-VAF summarization)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Newman et al. 2014, CAPP-Seq, Nat Med 20(5):548–554 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4016134/ | 2026-06-15 |
| 2 | US Patent 11,085,084 B2 (Poisson detection model; restates Avanzini 2020) | 2 | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11085084 | 2026-06-15 |
| 3 | Devonshire et al. 2014, cfDNA standardisation | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4182654/ | 2026-06-15 |
| 4 | Alcaide et al. 2020, ddPCR cfDNA, Sci Rep 10:12564 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7387491/ | 2026-06-15 |
| 5 | Pessoa et al. 2023, ctDNA/MRD review | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC10314661/ | 2026-06-15 |
| 6 | Antonello et al. 2024, CNAqc, Genome Biol 25:38 | 1 | https://doi.org/10.1186/s13059-024-03170-5 | 2026-06-15 |

### 1.2 Key Evidence Points

1. ctDNA detection follows a Poisson model with mean λ = n·d (n = sequenced genome equivalents, d = mutant allele fraction); single-reporter detection probability x = 1 − e^(−nd), k-reporter probability p = 1 − e^(−ndk) — Patent US11085084 (verbatim).
2. Worked example: n = 15,000 GE, d = 0.001 (0.1% VAF) ⇒ λ = 15 expected mutant molecules — Pessoa et al. 2023.
3. One haploid genome equivalent = 3.3 pg DNA; ≈ 303 haploid GE per ng cfDNA — Devonshire 2014; Alcaide 2020.
4. CAPP-Seq validated detection range 0.025%–10% mutant allele fraction; median pre-treatment ctDNA fraction ~0.1%; ctDNA level summarized as a fraction across SNV/indel reporters — Newman 2014.
5. For a clonal heterozygous SNV at a copy-neutral diploid locus, observed VAF v = π/2, so tumor fraction = 2·v — Antonello 2024 (CNAqc; m=1, n_tot=2 special case).

### 1.3 Documented Corner Cases

- **Poisson-limited low-input regime (λ < 3):** detection is stochastic; a true positive can be missed by sampling alone (Patent US11085084).
- **λ = 0** (n = 0 or d = 0): P(detect) = 1 − e⁰ = 0 (Poisson model).
- **Below assay LoD / background floor:** allele fractions below ~0.02–0.025% are not reliably distinguishable from error (Newman 2014).

### 1.4 Known Failure Modes / Pitfalls

1. Treating detection as deterministic in the low-burden regime hides Poisson false negatives — Patent US11085084.
2. Conflating allele fraction (VAF) with tumor fraction; for diploid heterozygous SNVs TF = 2·VAF, not VAF — Antonello 2024.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CtDnaDetectionProbability(int genomeEquivalents, double mutantAlleleFraction, int reporterCount)` | OncologyAnalyzer | Canonical | p = 1 − e^(−n·d·k) |
| `IsCtDnaDetected(int genomeEquivalents, double mutantAlleleFraction, int reporterCount, double minDetectionProbability)` | OncologyAnalyzer | Canonical | detect ⇔ λ≥1 AND p≥threshold |
| `CalculateTumorFraction(IEnumerable<VariantObservation>)` | OncologyAnalyzer | Canonical | 2 × mean clonal het VAF |
| `CalculateMeanVaf(IEnumerable<VariantObservation>)` | OncologyAnalyzer | Canonical | mean of altReads/totalReads |
| `HaploidGenomeEquivalents(double cfDnaNanograms)` | OncologyAnalyzer | Canonical | ng → GE via 3.3 pg/GE |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Detection probability p = 1 − e^(−n·d·k) ∈ [0, 1] | Yes | Patent US11085084 |
| INV-2 | p = 0 ⇔ λ = n·d·k = 0 (n=0 or d=0) | Yes | Poisson model, 1−e⁰=0 |
| INV-3 | p is non-decreasing in n, d, and k (strictly increasing while d>0,n>0) | Yes | Patent US11085084 (λ monotone) |
| INV-4 | CalculateTumorFraction = 2 × mean VAF; result clamped to [0,1] | Yes | Antonello 2024 (TF=2·v) |
| INV-5 | HaploidGenomeEquivalents(x ng) = x·1000/3.3; linear, GE(0)=0 | Yes | Devonshire 2014 (3.3 pg/GE) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | DetectionProbability worked example | n=15000, d=0.001, k=1 ⇒ λ=15 | p = 1 − e^(−15) = 0.99999969409… | Patent US11085084; Pessoa 2023 (n·d=15) |
| M2 | DetectionProbability k-reporter | n=1000, d=0.001, k=1 ⇒ λ=1 | p = 1 − e^(−1) = 0.6321205588… | Patent US11085084 (p=1−e^(−ndk)) |
| M3 | DetectionProbability multi-reporter raises p | n=1000, d=0.001, k=10 ⇒ λ=10 | p = 1 − e^(−10) = 0.9999546000… (> M2) | Patent US11085084 |
| M4 | DetectionProbability λ=0 | n=0 or d=0 | p = 0 | Poisson: 1−e⁰=0 |
| M5 | IsCtDnaDetected above LoD | n=15000, d=0.001, k=1, thr=0.95 | true (λ=15≥1, p≈1≥0.95) | Patent US11085084; Newman 2014 |
| M6 | IsCtDnaDetected below LoD (sub-molecule) | n=100, d=0.0001(0.01%), k=1 ⇒ λ=0.01 | false (λ<1; p=0.00995<0.95) | Newman 2014 (LoD); λ≥1 rule |
| M7 | CalculateTumorFraction = 2×mean VAF | two clonal het SNVs VAF 0.10 & 0.20 (mean 0.15) | TF = 0.30 | Antonello 2024 (TF=2·v) |
| M8 | CalculateTumorFraction clamps at 1.0 | mean VAF 0.6 ⇒ 2·0.6=1.2 | TF = 1.0 (clamped) | TF is a fraction ∈[0,1] |
| M9 | CalculateMeanVaf | reporters (alt 5/100, 30/100, 1/100) | mean = (0.05+0.30+0.01)/3 = 0.12 | Newman 2014 (per-reporter fraction) |
| M10 | HaploidGenomeEquivalents from ng | 1 ng | 1000/3.3 = 303.0303… GE | Devonshire 2014; Alcaide 2020 |
| M11 | HaploidGenomeEquivalents zero | 0 ng | 0 GE | linear, GE(0)=0 |
| M12 | DetectionProbability null variants → invalid d/n | d=−0.1 or d=1.1 | ArgumentOutOfRangeException | d is a fraction ∈[0,1] |
| M13 | CalculateTumorFraction null | null input | ArgumentNullException | domain guard |
| M14 | CalculateTumorFraction empty | empty input | ArgumentException | TF undefined with no variants |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | DetectionProbability negative n | n=−1 | ArgumentOutOfRangeException | GE count ≥ 0 |
| S2 | DetectionProbability k<1 | k=0 | ArgumentOutOfRangeException | reporters ≥ 1 |
| S3 | CalculateMeanVaf null/empty | null; empty | ArgumentNullException; ArgumentException | domain guards |
| S4 | CalculateTumorFraction VAF>0.5 | a variant VAF 0.6 | ArgumentOutOfRangeException (per-variant het cap) | het diploid SNV ≤0.5 |
| S5 | HaploidGenomeEquivalents negative | −1 ng | ArgumentOutOfRangeException | mass ≥ 0 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Detection probability monotone in d | p(d=0.002) > p(d=0.001) at fixed n,k | strictly greater | INV-3 |
| C2 | Detection probability bounded | p ∈ [0,1] for large λ | →1, never >1 | INV-1 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for ctDNA / tumor-fraction / liquid-biopsy tests: none exist. `OncologyAnalyzer` has no ctDNA methods prior to this unit. New canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CtDnaAnalysis_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M14, S1–S5, C1–C2 | ❌ Missing | brand-new unit; no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CtDnaAnalysis_Tests.cs` — all cases for this unit.
- **Remove:** none (no pre-existing ctDNA tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_CtDnaAnalysis_Tests.cs` | Canonical | 21 |

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
| 11 | M11 | ❌ Missing | Implemented | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented | ✅ Done |
| 14 | M14 | ❌ Missing | Implemented | ✅ Done |
| 15 | S1 | ❌ Missing | Implemented | ✅ Done |
| 16 | S2 | ❌ Missing | Implemented | ✅ Done |
| 17 | S3 | ❌ Missing | Implemented | ✅ Done |
| 18 | S4 | ❌ Missing | Implemented | ✅ Done |
| 19 | S5 | ❌ Missing | Implemented | ✅ Done |
| 20 | C1 | ❌ Missing | Implemented | ✅ Done |
| 21 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 21
**✅ Done:** 21 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | DetectionProbability_WorkedExample_15Molecules |
| M2 | ✅ | DetectionProbability_OneExpectedMolecule_Returns1MinusEInverse |
| M3 | ✅ | DetectionProbability_TenReporters_HigherThanSingle |
| M4 | ✅ | DetectionProbability_ZeroLambda_ReturnsZero |
| M5 | ✅ | IsCtDnaDetected_AboveLod_ReturnsTrue |
| M6 | ✅ | IsCtDnaDetected_BelowOneMolecule_ReturnsFalse |
| M7 | ✅ | CalculateTumorFraction_TwoClonalHetSnvs_ReturnsTwiceMeanVaf |
| M8 | ✅ | CalculateTumorFraction_MaxVafs_ClampedToOne (VAF 0.5+0.5 ⇒ 1.0; per-variant cap enforced) |
| M9 | ✅ | CalculateMeanVaf_ThreeReporters_ReturnsArithmeticMean |
| M10 | ✅ | HaploidGenomeEquivalents_OneNanogram_Returns303 |
| M11 | ✅ | HaploidGenomeEquivalents_Zero_ReturnsZero |
| M12 | ✅ | DetectionProbability_AlleleFractionOutOfRange_Throws |
| M13 | ✅ | CalculateTumorFraction_Null_Throws |
| M14 | ✅ | CalculateTumorFraction_Empty_Throws |
| S1 | ✅ | DetectionProbability_NegativeGenomeEquivalents_Throws |
| S2 | ✅ | DetectionProbability_ReporterCountBelowOne_Throws |
| S3 | ✅ | CalculateMeanVaf_NullAndEmpty_Throws |
| S4 | ✅ | CalculateTumorFraction_VafAboveHalf_Throws |
| S5 | ✅ | HaploidGenomeEquivalents_Negative_Throws |
| C1 | ✅ | DetectionProbability_IncreasingAlleleFraction_IsMonotone |
| C2 | ✅ | DetectionProbability_LargeLambda_BoundedByOne |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Detection-decision rule combines the source-exact probability p = 1 − e^(−ndk) with λ≥1 and a caller-supplied probability threshold (default 0.95). Only the boolean `IsCtDnaDetected` depends on the 0.95 default; `CtDnaDetectionProbability` is fully source-exact. | M5, M6 |

---

## 7. Open Questions / Decisions

1. The Avanzini et al. 2020 *Science Advances* primary (DOI 10.1126/sciadv.abc4308) was paywalled (HTTP 403) this session; the identical detection equations were retrieved verbatim from US Patent 11,085,084 and corroborated by the Pessoa 2023 worked molecule count (n·d = 15). The probability formula is therefore source-backed; no value depends on the unretrieved primary.
2. `AnalyzeFragmentSizeDistribution(bamFile)` from the Registry method list is **out of scope** for this unit: it requires BAM parsing infrastructure (not present) and is a fragmentomics analysis, not one of the quantitative ctDNA pieces specified. Recorded in the algorithm doc §5.3 "Not implemented".
