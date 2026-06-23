# Test Specification: ONCO-CHIP-001

**Test Unit ID:** ONCO-CHIP-001
**Area:** Oncology
**Algorithm:** Clonal Hematopoiesis (CHIP) Filtering for cfDNA Liquid Biopsy
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Steensma et al. (2015) *Blood* 126(1):9–16 — CHIP definition | 1 | https://doi.org/10.1182/blood-2015-03-631747 (PMC4624443) | 2026-06-15 |
| 2 | Genovese et al. (2014) *NEJM* 371(26):2477–2487 — CH driver genes | 1 | https://doi.org/10.1056/NEJMoa1409405 (PMC4290021) | 2026-06-15 |
| 3 | Razavi et al. (2019) *Nat Med* 25:1928–1937 — matched-WBC cfDNA filtering | 1 | https://doi.org/10.1038/s41591-019-0652-7 | 2026-06-15 |
| 4 | Arango-Argoty et al. (2025) *NPJ Precis Oncol* 9:147 — matched-WBC rule, top-3 genes | 3 | https://doi.org/10.1038/s41698-025-00921-w (PMC12092662) | 2026-06-15 |
| 5 | Bolton et al. (2020) *Nat Genet* 52(11):1219–1226 — strict matched-WBC origin rule | 1 | https://doi.org/10.1038/s41588-020-00710-0 (PMC7891089) | 2026-06-23 |

### 1.2 Key Evidence Points

1. CHIP = a somatic mutation in a gene recurrently mutated in hematologic malignancies, at VAF ≥ 2% (0.02) in blood, without diagnostic criteria for a hematologic malignancy — Steensma 2015 ("the mutant allele fraction must be ≥2% in the peripheral blood").
2. Canonical CH driver genes: DNMT3A, TET2, ASXL1, PPM1D, JAK2, SF3B1 (and TP53, SRSF2) — Genovese 2014 ("Four genes (DNMT3A, TET2, ASXL1, and PPM1D) had disproportionately high numbers of somatic mutations"); top-3 DNMT3A/TET2/ASXL1 — Arango-Argoty 2025.
3. In cfDNA liquid biopsy, CH variants are the dominant non-tumor confounder (81.6% controls / 53.2% cancer patients) and are identified by matched white-blood-cell (WBC) sequencing — a cfDNA variant also present in matched WBC is CH-derived, not tumor — Razavi 2019; Arango-Argoty 2025.
4. The gene+VAF heuristic flags *candidate* CHIP; matched-WBC subtraction is the definitive origin test (VAF–origin relationship "remains unclear") — Arango-Argoty 2025.
5. Strict matched-WBC origin call: a variant is WBC/CH origin when its matched-WBC VAF is ≥ 2% with ≥ 10 supporting reads AND its WBC VAF is ≥ 2× the tumour VAF (1.5× if the biopsy site is a lymph node); otherwise tumour-derived — Bolton 2020 (Methods: "a variant allele fraction of at least 2% and at least 10 supporting reads"; blood VAF "at least twice that in the tumor or 1.5 times … if the tumor biopsy site was a lymph node").

### 1.3 Documented Corner Cases

- Driver-gene mutation with VAF < 0.02 does NOT meet the CHIP definition — Steensma 2015.
- A cfDNA variant absent from matched WBC is retained as candidate tumor even in a CHIP gene — Razavi 2019.
- A cfDNA variant present in matched WBC is removed regardless of gene — Razavi 2019.

### 1.4 Known Failure Modes / Pitfalls

1. Omitting matched-WBC subtraction mis-calls CH variants as tumor (dominant false positive) — Razavi 2019.
2. Using a fixed VAF cutoff alone over-/under-calls origin because the VAF–origin relationship is unclear — Arango-Argoty 2025.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `IdentifyCHIPVariants(variants, chipGenes?, minVaf?)` | OncologyAnalyzer | Canonical | Gene+VAF heuristic (Steensma 2015) |
| `FilterCHIP(variants, whiteBloodCellVariants, ...)` | OncologyAnalyzer | Canonical | Matched-WBC subtraction (Razavi 2019) |
| `CallVariantOrigin(variants, whiteBloodCellObservations, wbcVafFold?, chipMinWbcVaf?, minWbcAltReads?)` | OncologyAnalyzer | Canonical | Strict matched-WBC origin call (Bolton 2020) |
| `IsCanonicalChipGene(gene, chipGenes?)` | OncologyAnalyzer | Internal | Case-insensitive gene membership |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | A variant is flagged CHIP ⟺ its gene ∈ CHIP set AND its VAF ≥ 0.02 | Yes | Steensma 2015 |
| INV-2 | VAF threshold is inclusive at 0.02 (≥, not >) | Yes | Steensma 2015 ("≥2%") |
| INV-3 | `FilterCHIP` output ⊆ input, preserving input order | Yes | deterministic-filter contract |
| INV-4 | A cfDNA variant present in matched WBC is removed regardless of gene | Yes | Razavi 2019 |
| INV-5 | A cfDNA variant absent from matched WBC is retained | Yes | Razavi 2019 |
| INV-6 | Gene membership comparison is case-insensitive | Yes | HGNC symbol convention (ASSUMPTION) |
| INV-7 | `CallVariantOrigin` calls CHIP ⟺ WBC VAF ≥ τ_w (0.02) AND WBC reads ≥ ρ (10) AND WBC VAF ≥ φ·tumour VAF (φ=2.0; 1.5 lymph node), at the matched locus | Yes | Bolton 2020 |
| INV-8 | `CallVariantOrigin` calls tumour when the locus is absent from matched WBC (even a CH driver gene) | Yes | Bolton 2020 (no WBC evidence ⇒ tumour) |
| INV-9 | `CallVariantOrigin` emits exactly one call per input variant, in input order | Yes | single-pass contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | CHIP gene at VAF ≥ 0.02 | DNMT3A, VAF 0.05 | flagged CHIP | Steensma 2015 |
| M2 | CHIP gene below threshold | DNMT3A, VAF 0.01 | not CHIP | Steensma 2015 |
| M3 | Non-CHIP gene high VAF | EGFR, VAF 0.30 | not CHIP | Steensma 2015 (driver-gene rule) |
| M4 | VAF exactly at 0.02 | TET2, VAF 0.02 | flagged CHIP (inclusive ≥) | Steensma 2015 ("≥2%") |
| M5 | All canonical genes recognized | each of DNMT3A,TET2,ASXL1,TP53,JAK2,SF3B1,SRSF2,PPM1D at VAF 0.10 | all flagged CHIP | Genovese 2014 / Steensma 2015 |
| M6 | Caller-supplied panel | gene "ABCX" with custom panel {ABCX} | flagged CHIP | Razavi/Arango-Argoty (caller panel) |
| M7 | Default panel excludes non-member | "ABCX" with default panel | not CHIP | Steensma 2015 set |
| M8 | FilterCHIP removes WBC-matched | EGFR variant present in matched WBC | removed | Razavi 2019 |
| M9 | FilterCHIP retains WBC-absent | EGFR variant absent from matched WBC | retained | Razavi 2019 |
| M10 | FilterCHIP removes CHIP-gene variant by gene+VAF even if not in WBC | DNMT3A VAF 0.05, not in WBC | removed (CHIP heuristic) | Steensma 2015 |
| M11 | FilterCHIP retains sub-threshold CHIP-gene variant absent from WBC | DNMT3A VAF 0.01, not in WBC | retained | Steensma 2015 |
| M12 | Mixed panel filtering | tumor EGFR(absent), CH DNMT3A(0.06), WBC-matched KRAS | only EGFR retained | Razavi 2019 + Steensma 2015 |
| OC1 | CallVariantOrigin: present in WBC at comparable VAF | tumour 0.10, WBC 0.30, 40 reads | Chip | Bolton 2020 |
| OC2 | CallVariantOrigin: absent from WBC (CH gene) | DNMT3A 0.40, no WBC observation | Tumor | Bolton 2020 (no WBC evidence) |
| OC3 | CallVariantOrigin: WBC reads below floor | tumour 0.10, WBC 0.30, 9 reads | Tumor | Bolton 2020 (≥ 10 reads) |
| OC4 | CallVariantOrigin: all thresholds at boundary | tumour 0.01, WBC 0.02, 10 reads | Chip (inclusive) | Bolton 2020 |
| OC5 | CallVariantOrigin: WBC VAF below 2% | tumour 0.005, WBC 0.015, 50 reads | Tumor | Bolton 2020 (VAF ≥ 2%) |
| OC6 | CallVariantOrigin: WBC VAF below fold | tumour 0.30, WBC 0.40, 40 reads | Tumor | Bolton 2020 (≥ 2×) |
| OC7 | CallVariantOrigin: lymph-node fold changes call | tumour 0.25, WBC 0.40, 40 reads | Tumor (2.0) / Chip (1.5) | Bolton 2020 (1.5× lymph node) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitive gene | "dnmt3a" lower-case, VAF 0.05 | flagged CHIP | INV-6 |
| S2 | WBC presence by alt-read evidence | locus with 0 alt reads in WBC ⇒ absent | retained | Wan 2020 alt-read rule |
| S3 | Custom minVaf | minVaf 0.10, DNMT3A VAF 0.05 | not CHIP | threshold configurable |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Order preservation | 3 retained variants | same relative order | INV-3 |
| C2 | CallVariantOrigin order/multiplicity | mixed CHIP + tumour variants | one call per input, input order | INV-9 |

### 4.4 Validation / Error Tests

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| V1 | null variants (IdentifyCHIPVariants) | null input | ArgumentNullException | repo convention |
| V2 | null variants (FilterCHIP) | null cfDNA input | ArgumentNullException | repo convention |
| V3 | null matched WBC | null WBC collection | ArgumentNullException | repo convention |
| V4 | minVaf out of (0,1] | minVaf 0 or 1.5 | ArgumentOutOfRangeException | domain |
| V5 | empty cfDNA input | empty list | empty result | repo convention |
| V6 | empty matched WBC | nothing to subtract; only gene/VAF heuristic applies | CH-gene variants removed, rest kept | Razavi 2019 |
| V7 | CallVariantOrigin invalid args | null variants/WBC, fold<1, chipMinWbcVaf∉(0,1], minWbcAltReads<1, empty variants | throws / empty result | domain + repo convention |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- New Test Unit (2026-06-15). No prior `FilterCHIP` / `IdentifyCHIPVariants` implementation or tests existed in `OncologyAnalyzer` or `tests/Seqeron/Seqeron.Genomics.Tests/`. Canonical test file: `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_FilterCHIP_Tests.cs`.
- 2026-06-23 limitation-fix: added the strict matched-WBC origin-calling method `CallVariantOrigin` (Bolton 2020) and its tests (OC1–OC9 region) to the same canonical file; the gene+VAF `FilterCHIP` heuristic is unchanged and remains the no-WBC fallback.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M12, S1–S3, C1, V1–V6 | ✅ Covered | Implemented at unit creation (2026-06-15) |
| OC1–OC7, C2, V7 | ❌ Missing | New strict-origin cases for `CallVariantOrigin` (2026-06-23) |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_FilterCHIP_Tests.cs` — all cases for both methods, `#region` per method.
- **Remove:** none (new unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_FilterCHIP_Tests.cs` | canonical | 35 |

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
| 9 | M9 | ❌ Missing | implemented | ✅ Done |
| 10 | M10 | ❌ Missing | implemented | ✅ Done |
| 11 | M11 | ❌ Missing | implemented | ✅ Done |
| 12 | M12 | ❌ Missing | implemented | ✅ Done |
| 13 | S1 | ❌ Missing | implemented | ✅ Done |
| 14 | S2 | ❌ Missing | implemented | ✅ Done |
| 15 | S3 | ❌ Missing | implemented | ✅ Done |
| 16 | C1 | ❌ Missing | implemented | ✅ Done |
| 17 | V1 | ❌ Missing | implemented | ✅ Done |
| 18 | V2 | ❌ Missing | implemented | ✅ Done |
| 19 | V3 | ❌ Missing | implemented | ✅ Done |
| 20 | V4 | ❌ Missing | implemented | ✅ Done |
| 21 | V5 | ❌ Missing | implemented | ✅ Done |
| 22 | V6 | ❌ Missing | implemented | ✅ Done |
| 23 | OC1 | ❌ Missing | implemented (2026-06-23) | ✅ Done |
| 24 | OC2 | ❌ Missing | implemented (2026-06-23) | ✅ Done |
| 25 | OC3 | ❌ Missing | implemented (2026-06-23) | ✅ Done |
| 26 | OC4 | ❌ Missing | implemented (2026-06-23) | ✅ Done |
| 27 | OC5 | ❌ Missing | implemented (2026-06-23) | ✅ Done |
| 28 | OC6 | ❌ Missing | implemented (2026-06-23) | ✅ Done |
| 29 | OC7 | ❌ Missing | implemented (2026-06-23) | ✅ Done |
| 30 | C2 | ❌ Missing | implemented (2026-06-23) | ✅ Done |
| 31 | V7 | ❌ Missing | implemented (2026-06-23) | ✅ Done |

**Total items:** 31
**✅ Done:** 31 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | IdentifyCHIPVariants_DriverGeneAboveThreshold_FlagsChip |
| M2 | ✅ | IdentifyCHIPVariants_DriverGeneBelowThreshold_NotChip |
| M3 | ✅ | IdentifyCHIPVariants_NonChipGeneHighVaf_NotChip |
| M4 | ✅ | IdentifyCHIPVariants_VafExactlyAtThreshold_FlagsChip |
| M5 | ✅ | IdentifyCHIPVariants_AllCanonicalGenes_AllFlagged |
| M6 | ✅ | IdentifyCHIPVariants_CallerSuppliedPanel_FlagsCustomGene |
| M7 | ✅ | IdentifyCHIPVariants_DefaultPanel_ExcludesNonMember |
| M8 | ✅ | FilterCHIP_VariantInMatchedWbc_Removed |
| M9 | ✅ | FilterCHIP_VariantAbsentFromWbc_Retained |
| M10 | ✅ | FilterCHIP_ChipGeneVariantNotInWbc_RemovedByHeuristic |
| M11 | ✅ | FilterCHIP_SubThresholdChipGeneAbsentFromWbc_Retained |
| M12 | ✅ | FilterCHIP_MixedPanel_OnlyTumorVariantRetained |
| S1 | ✅ | IdentifyCHIPVariants_LowercaseGene_FlagsChip |
| S2 | ✅ | FilterCHIP_WbcLocusZeroAltReads_TreatedAbsent |
| S3 | ✅ | IdentifyCHIPVariants_CustomMinVaf_ExcludesBelowCustom |
| C1 | ✅ | FilterCHIP_PreservesInputOrder |
| V1 | ✅ | IdentifyCHIPVariants_NullVariants_Throws |
| V2 | ✅ | FilterCHIP_NullVariants_Throws |
| V3 | ✅ | FilterCHIP_NullWbc_Throws |
| V4 | ✅ | IdentifyCHIPVariants_MinVafOutOfRange_Throws |
| V5 | ✅ | FilterCHIP_EmptyVariants_ReturnsEmpty |
| V6 | ✅ | FilterCHIP_EmptyWbc_AppliesGeneHeuristicOnly |
| OC1 | ✅ | CallVariantOrigin_PresentInWbcAtComparableVaf_CallsChip |
| OC2 | ✅ | CallVariantOrigin_AbsentFromWbc_CallsTumor |
| OC3 | ✅ | CallVariantOrigin_WbcReadsBelowFloor_CallsTumor |
| OC4 | ✅ | CallVariantOrigin_AllThresholdsExactlyAtBoundary_CallsChip |
| OC5 | ✅ | CallVariantOrigin_WbcVafBelowChipMinimum_CallsTumor |
| OC6 | ✅ | CallVariantOrigin_WbcVafBelowFold_CallsTumor |
| OC7 | ✅ | CallVariantOrigin_LymphNodeFold_ChangesCallVersusDefault |
| C2 | ✅ | CallVariantOrigin_MixedVariants_PreservesOrderAndCallsEach |
| V7 | ✅ | CallVariantOrigin_InvalidArguments_Throw |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Default canonical CHIP gene set {DNMT3A,TET2,ASXL1,TP53,JAK2,SF3B1,SRSF2,PPM1D}; case-insensitive | M5, M7, S1, INV-6 |
| 2 | "Present in matched WBC" = ≥ 1 alt read at the same locus (Wan 2020 convention), configurable | M8, S2 |

---

## 7. Open Questions / Decisions

1. The by-area Registry names class `LiquidBiopsyAnalyzer`; the actual repository places all ONCO units in `OncologyAnalyzer`. Decision: implement in `OncologyAnalyzer` (matches every sibling ONCO unit and the task directive). Noted as a checklist conflict.
2. Universal per-locus alt-read cutoff for WBC presence is assay-specific (Wan 2020); exposed as a configurable parameter rather than a hard constant.
