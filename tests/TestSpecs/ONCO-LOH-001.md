# Test Specification: ONCO-LOH-001

**Test Unit ID:** ONCO-LOH-001
**Area:** Oncology
**Algorithm:** Loss of Heterozygosity (LOH) detection and HRD-LOH genomic-scar score
**Status:** ☑ Validated (2026-06-16 — Stage A PASS, Stage B PASS-WITH-NOTES, CLEAN; see docs/Validation/reports/ONCO-LOH-001.md)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-16

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Abkevich et al. (2012), Br J Cancer 107(10):1776–1782 | 1 | https://doi.org/10.1038/bjc.2012.451 (PMID 23047548) | 2026-06-14 |
| 2 | scarHRD (sztup) `calc.hrd.R` reference implementation | 3 | https://raw.githubusercontent.com/sztup/scarHRD/master/R/calc.hrd.R | 2026-06-14 |
| 3 | scarHRD (sztup) `scar_score.R` / `scarHRD.md` (input format, 15 Mb cutoff) | 3 | https://github.com/sztup/scarHRD/blob/master/scarHRD.md | 2026-06-14 |
| 4 | oncoscanR (Christinat) `score_loh` doc | 3 | https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. HRD-LOH score = "the number of 15 Mb exceeding LOH regions which do not cover the whole chromosome." — scarHRD.md; Abkevich 2012 ("number of these [intermediate-size] regions").
2. A segment is LOH when minor allele copy number == 0 AND major allele copy number != 0 (`segSamp[,nB]==0 & segSamp[,nA]!=0`). — scarHRD `calc.hrd.R`.
3. Size filter is strict: `end - start > 15e6` (15,000,000 bp); length = end − start. — scarHRD `calc.hrd.R`.
4. Whole-chromosome exclusion: a chromosome where ALL segments have minor==0 is dropped from the count (`chrDel`). — scarHRD `calc.hrd.R`; Abkevich "< whole chromosome".
5. Input columns: chromosome, start, end, total CN, A (major) CN, B (minor) CN. — scarHRD `scar_score.R`.
6. Adjacent/overlapping LOH segments are merged before the size filter. — oncoscanR `score_loh`.

### 1.3 Documented Corner Cases

- Homozygous deletion (minor==0 & major==0) is not LOH (Evidence §Corner Cases / `nA != 0`).
- Segment of length exactly 15 Mb is not counted (strict `>`).
- Whole-chromosome LOH excluded.
- Heterozygous-retained segment (minor != 0) is not LOH.

### 1.4 Known Failure Modes / Pitfalls

1. Treating length as `end - start + 1` instead of `end - start` would shift the boundary — Evidence §Corner Cases (source: scarHRD `calc.hrd.R`).
2. Counting whole-chromosome LOH inflates the score (it does not correlate with HRD) — Abkevich 2012.
3. Counting homozygous deletions as LOH — scarHRD `nA != 0` clause.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DetectLOH(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Canonical | Returns qualifying HRD-LOH regions and the HRD-LOH count |
| `CalculateHrdLohScore(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Canonical | HRD-LOH count = number of qualifying regions (Abkevich/scarHRD) |
| `CalculateLOHFraction(IEnumerable<AlleleSpecificSegment>, chromosome)` | OncologyAnalyzer | Canonical | Length-weighted LOH fraction of one chromosome ∈ [0,1] |
| `IsLohSegment(AlleleSpecificSegment)` | OncologyAnalyzer | Internal | Tested indirectly via DetectLOH |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | HRD-LOH score ≥ 0 (count of regions) | Yes | scarHRD `nrow(segLOH)` ≥ 0 |
| INV-2 | 0 ≤ LOH_fraction ≤ 1 per chromosome | Yes | Registry invariant; LOH length ⊆ total length |
| INV-3 | A counted LOH region has minor CN == 0 and major CN != 0 | Yes | scarHRD `calc.hrd.R` |
| INV-4 | A counted LOH region has length > 15,000,000 bp (strict) | Yes | scarHRD `> sizelimit1`, `sizelimitLOH=15e6` |
| INV-5 | No counted LOH region lies on a whole-chromosome-LOH chromosome | Yes | scarHRD `chrDel`; Abkevich "< whole chromosome" |
| INV-6 | Count is independent of input segment order | Yes | Set-based count (per-chromosome aggregation) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Synthetic 7-segment dataset | Full Evidence dataset (chr1–5) | HRD-LOH score = 1 | scarHRD `calc.hrd` rule (Evidence Test Datasets) |
| M2 | LOH segment > 15 Mb counted | chr1 20 Mb, minor=0, major=1, plus a het segment on chr1 | counted (score contribution 1) | INV-3, INV-4 |
| M3 | Length exactly 15 Mb not counted | minor=0, major=1, length=15,000,000 | not counted | INV-4 (`> 15e6`) |
| M4 | Homozygous deletion not LOH | minor=0, major=0, length 30 Mb | not counted | INV-3 (`major != 0`) |
| M5 | Heterozygous-retained not LOH | minor=1, major=1, length 40 Mb | not counted | INV-3 (`minor == 0`) |
| M6 | Whole-chromosome LOH excluded | single chr with all segments minor=0, length 16 Mb | not counted | INV-5 (`chrDel`) |
| M7 | LOH on chr with a non-LOH segment IS counted | chr has one >15Mb LOH segment + one het segment | counted | INV-5 (not whole-chromosome) |
| M8 | CalculateLOHFraction partial | chr1: (0–20M minor0)+(20M–60M minor1) | 20M/60M = 0.3333333333 | INV-2, LOH length/total |
| M9 | CalculateLOHFraction none | chr2: single het 0–50M | 0.0 | INV-2 |
| M10 | CalculateLOHFraction full | chr3: single LOH 0–40M minor0 | 1.0 | INV-2 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Adjacent LOH merge crosses 15 Mb | two adjacent 8 Mb LOH segments same chr (with a non-LOH segment present) | merged → 16 Mb → counted as 1 | oncoscanR merge rule |
| S2 | Empty input | no segments | score 0; fraction 0 | empty domain |
| S3 | Null input | null enumerable | ArgumentNullException | input validation |
| S4 | Negative / zero-length segment | end ≤ start | ArgumentException | invalid segment |
| S5 | CalculateLOHFraction unknown chromosome | chromosome absent from segments | 0.0 | no covered length |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Order invariance | shuffle the M1 dataset | score unchanged = 1 | INV-6 |
| C2 | DetectLOH region payload | returned regions match counted segments | one region (chr1) with correct coords | DetectLOH contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- New unit: no prior LOH tests existed in `tests/Seqeron/Seqeron.Genomics.Tests/`. All planned cases start ❌ Missing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new |
| M2 | ❌ Missing | new |
| M3 | ❌ Missing | new |
| M4 | ❌ Missing | new |
| M5 | ❌ Missing | new |
| M6 | ❌ Missing | new |
| M7 | ❌ Missing | new |
| M8 | ❌ Missing | new |
| M9 | ❌ Missing | new |
| M10 | ❌ Missing | new |
| S1 | ❌ Missing | new |
| S2 | ❌ Missing | new |
| S3 | ❌ Missing | new |
| S4 | ❌ Missing | new |
| S5 | ❌ Missing | new |
| C1 | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectLOH_Tests.cs` — all DetectLOH / CalculateHrdLohScore / CalculateLOHFraction tests.
- **Remove:** (none — new unit)

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_DetectLOH_Tests.cs` | canonical (all cases) | 18 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | `DetectLOH_EvidenceDataset_ScoreIsOne` | ✅ Done |
| 2 | M2 | ❌ Missing | `DetectLOH_LohSegmentOver15Mb_IsCounted` | ✅ Done |
| 3 | M3 | ❌ Missing | `DetectLOH_SegmentExactly15Mb_IsNotCounted` | ✅ Done |
| 4 | M4 | ❌ Missing | `DetectLOH_HomozygousDeletion_IsNotLoh` | ✅ Done |
| 5 | M5 | ❌ Missing | `DetectLOH_HeterozygousRetained_IsNotLoh` | ✅ Done |
| 6 | M6 | ❌ Missing | `DetectLOH_WholeChromosomeLoh_IsExcluded` | ✅ Done |
| 7 | M7 | ❌ Missing | `DetectLOH_LohWithNonLohSegmentOnSameChromosome_IsCounted` | ✅ Done |
| 8 | M8 | ❌ Missing | `CalculateLOHFraction_PartialChromosome_ReturnsLengthWeightedFraction` | ✅ Done |
| 9 | M9 | ❌ Missing | `CalculateLOHFraction_NoLoh_ReturnsZero` | ✅ Done |
| 10 | M10 | ❌ Missing | `CalculateLOHFraction_FullLoh_ReturnsOne` | ✅ Done |
| 11 | INV-2 | ❌ Missing | `CalculateLOHFraction_MixedChromosome_IsWithinUnitInterval` | ✅ Done |
| 12 | S1 | ❌ Missing | `DetectLOH_AdjacentLohSegmentsMergeAcross15Mb_CountedAsOne` | ✅ Done |
| 13 | S2 | ❌ Missing | `DetectLOH_EmptyInput_ScoreZero` | ✅ Done |
| 14 | S3 | ❌ Missing | `DetectLOH_NullInput_Throws` | ✅ Done |
| 15 | S4 | ❌ Missing | `DetectLOH_NonPositiveLength_Throws` + `DetectLOH_NegativeCopyNumber_Throws` | ✅ Done |
| 16 | S5 | ❌ Missing | `CalculateLOHFraction_UnknownChromosome_ReturnsZero` | ✅ Done |
| 17 | C1 | ❌ Missing | `DetectLOH_ShuffledInput_ScoreUnchanged` | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `DetectLOH_EvidenceDataset_ScoreIsOne` |
| M2 | ✅ Covered | `DetectLOH_LohSegmentOver15Mb_IsCounted` |
| M3 | ✅ Covered | `DetectLOH_SegmentExactly15Mb_IsNotCounted` |
| M4 | ✅ Covered | `DetectLOH_HomozygousDeletion_IsNotLoh` |
| M5 | ✅ Covered | `DetectLOH_HeterozygousRetained_IsNotLoh` |
| M6 | ✅ Covered | `DetectLOH_WholeChromosomeLoh_IsExcluded` |
| M7 | ✅ Covered | `DetectLOH_LohWithNonLohSegmentOnSameChromosome_IsCounted` |
| M8 | ✅ Covered | `CalculateLOHFraction_PartialChromosome_ReturnsLengthWeightedFraction` |
| M9 | ✅ Covered | `CalculateLOHFraction_NoLoh_ReturnsZero` |
| M10 | ✅ Covered | `CalculateLOHFraction_FullLoh_ReturnsOne` |
| INV-2 | ✅ Covered | `CalculateLOHFraction_MixedChromosome_IsWithinUnitInterval` |
| S1 | ✅ Covered | `DetectLOH_AdjacentLohSegmentsMergeAcross15Mb_CountedAsOne` |
| S2 | ✅ Covered | `DetectLOH_EmptyInput_ScoreZero` |
| S3 | ✅ Covered | `DetectLOH_NullInput_Throws` |
| S4 | ✅ Covered | `DetectLOH_NonPositiveLength_Throws`, `DetectLOH_NegativeCopyNumber_Throws` |
| S5 | ✅ Covered | `CalculateLOHFraction_UnknownChromosome_ReturnsZero` |
| C1 | ✅ Covered | `DetectLOH_ShuffledInput_ScoreUnchanged` |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | LOH fraction = length-weighted fraction of a chromosome under LOH (Registry invariant; aggregation is a definitional choice, segment criterion is source-backed) | M8–M10, S5, INV-2 |
| 2 | Input is allele-specific copy-number segments (scarHRD `seg` shape), not raw VCF text | All |

---

## 7. Open Questions / Decisions

1. **Decision:** `DetectLOH` operates on allele-specific segments rather than `(tumorVcf, normalVcf)` text because the retrievable Abkevich/scarHRD/oncoscanR algorithm is defined over segmented allele-specific copy number; raw segmentation/BAF is upstream (ONCO-CNA-001). Recorded in Evidence Assumptions.
2. **Decision:** segment merging of adjacent same-state LOH is applied before the size filter, per oncoscanR; absent merging, two adjacent < 15 Mb LOH pieces could be wrongly missed.
