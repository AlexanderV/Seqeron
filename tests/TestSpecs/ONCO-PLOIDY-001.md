# Test Specification: ONCO-PLOIDY-001

**Test Unit ID:** ONCO-PLOIDY-001
**Area:** Oncology
**Algorithm:** Tumor Ploidy Estimation (length-weighted mean segment copy number) + Whole-Genome-Doubling detection
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Patchwork (Genome Biology) — verbatim ploidy definition | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4053982/ | 2026-06-14 |
| 2 | Van Loo et al. ASCAT (PNAS 2010) | 1 | https://doi.org/10.1073/pnas.1009843107 | 2026-06-14 |
| 3 | Bielski et al. (Nature Genetics 2018) — WGD | 1 | https://doi.org/10.1038/s41588-018-0165-1 | 2026-06-14 |
| 4 | facets-suite `copy-number-scores.R` `is_genome_doubled` (reference impl, PMID 30013179) | 3 | https://github.com/mskcc/facets-suite/blob/master/R/copy-number-scores.R | 2026-06-14 |

### 1.2 Key Evidence Points

1. Average tumour ploidy is "the average total copy number of all genomic segments weighted by segment length" → ψ = Σ(CN_i · L_i) / Σ(L_i) — Patchwork, PMC4053982.
2. Ploidy is reported on the n-scale (2n = diploid); ">2.7n" marks aneuploidy / near-triploid genomes — Van Loo et al. 2010, PNAS abstract.
3. WGD is called when the autosome-restricted fraction of genome with **major copy number ≥ 2** is strictly greater than 0.5: `frac_elevated_mcn > treshold` (treshold = 0.5) — facets-suite `is_genome_doubled` (PMID 30013179).
4. Major copy number `mcn = tcn - lcn` (total − minor); WGD uses the major allele CN ≥ 2, not total CN ≥ 2 — facets-suite `parse_segs`.

### 1.3 Documented Corner Cases

- Empty segment set → Σ(L) = 0, ploidy undefined (Patchwork weighted mean).
- Segment with Length ≤ 0 or negative copy number → invalid input.
- WGD threshold is strict (`>` 0.5): exactly 50% at major CN ≥ 2 is NOT doubled (facets-suite).
- WGD uses **major** CN: an all-1:1 genome (total CN 2) is NOT doubled (facets-suite `mcn >= 2`).

### 1.4 Known Failure Modes / Pitfalls

1. Using a plain (unweighted) mean of per-segment copy numbers instead of length-weighting it — Patchwork ("weighted by segment length").
2. Calling WGD on total CN ≥ 2 instead of major CN ≥ 2 (would mis-call balanced diploids) — facets-suite.
3. Using `≥ 0.5` instead of strict `> 0.5` for the fraction — facets-suite `> treshold`.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `EstimatePloidy(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Canonical | ψ = Σ(CN·L)/Σ(L), CN = Major+Minor |
| `DetectWholeGenomeDoubling(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Canonical | facets-suite rule: frac(major CN ≥ 2 by length) > 0.5 |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | ploidy > 0 for any non-empty valid segment set with at least one positive copy number | Yes | Patchwork weighted mean; registry invariant |
| INV-2 | a genome of pure 1:1 (total CN 2) segments has ploidy exactly 2.0 | Yes | n-scale 2n diploid (ASCAT/Patchwork) |
| INV-3 | ploidy is length-weighted: min(CN_i) ≤ ψ ≤ max(CN_i) | Yes | weighted mean lies within the value range (Patchwork) |
| INV-4 | WGD = true ⇔ (Σ length where major CN ≥ 2) / (Σ length) > 0.5 over the supplied segments | Yes | facets-suite `is_genome_doubled` |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Ploidy worked example | CN 2 (1:1)/4 (2:2)/3 (2:1), lengths 100/100/50 Mb | ψ = 750M/250M = 3.0 | Patchwork weighted mean |
| M2 | Pure diploid | all 1:1 segments | ψ = 2.0 exactly | n-scale 2n |
| M3 | Length weighting dominates | long 1:1 (300 Mb) + short 2:2 (10 Mb) | ψ = (2·300+4·10)/310 = 640/310 ≈ 2.0645 | "weighted by segment length" |
| M4 | Single segment | one 2:1 segment (total 3) | ψ = 3.0 | weighted mean of one value |
| M5 | Empty segments → reject | no segments | ArgumentException | Σ(L)=0 undefined |
| M6 | Invalid segment length → reject | End ≤ Start (Length ≤ 0) | ArgumentException | invalid input |
| M7 | Negative copy number → reject | Major or Minor < 0 | ArgumentException | invalid input |
| M8 | WGD 60% elevated → true | 60% of length at major CN ≥ 2, 40% at 1:1 | true | facets-suite > 0.5 |
| M9 | WGD exactly 50% → false | half length major CN ≥ 2, half 1:1 | false | strict `>` 0.5 |
| M10 | WGD 40% elevated → false | 40% at major CN ≥ 2 | false | facets-suite ≤ 0.5 |
| M11 | WGD all 1:1 (total 2) → false | every segment major CN = 1 | false | mcn >= 2 (not total) |
| M12 | WGD all major ≥ 2 → true | every segment 2:0/2:2 | true | frac = 1.0 > 0.5 |
| M13 | WGD empty/invalid → reject | empty set; Length ≤ 0; negative CN | ArgumentException | shared validation |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | WGD LOH counts as elevated | 2:0 segments (major 2, minor 0) over >50% | true | major CN, not heterozygosity |
| S2 | Ploidy with a CN-0 (homozygous deletion) segment | 0 (0:0)/4 (2:2) equal lengths | ψ = 2.0 | weighted mean includes zeros |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Near-triploid genome ploidy | mostly CN-3 segments | ψ ≈ 3 (>2.7n aneuploid direction) | Van Loo aneuploidy direction |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- New unit. No existing `EstimatePloidy` / `DetectWholeGenomeDoubling` methods in `OncologyAnalyzer.cs` and no `OncologyAnalyzer_EstimatePloidy_Tests.cs`. Reuses the existing `AlleleSpecificSegment` record (ONCO-LOH-001 / ONCO-HRD-001).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M13, S1–S2, C1 | ❌ Missing | brand-new unit, no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_EstimatePloidy_Tests.cs` — all ploidy + WGD tests.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_EstimatePloidy_Tests.cs` | canonical | 16 |

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
| 13 | M13 | ❌ Missing | implemented | ✅ Done |
| 14 | S1 | ❌ Missing | implemented | ✅ Done |
| 15 | S2 | ❌ Missing | implemented | ✅ Done |
| 16 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | EstimatePloidy_WorkedExample_ReturnsThree |
| M2 | ✅ Covered | EstimatePloidy_PureDiploid_ReturnsTwo |
| M3 | ✅ Covered | EstimatePloidy_LongDiploidShortAmplified_IsLengthWeighted |
| M4 | ✅ Covered | EstimatePloidy_SingleSegment_ReturnsItsTotalCopyNumber |
| M5 | ✅ Covered | EstimatePloidy_EmptySegments_Throws |
| M6 | ✅ Covered | EstimatePloidy_NonPositiveLength_Throws |
| M7 | ✅ Covered | EstimatePloidy_NegativeCopyNumber_Throws |
| M8 | ✅ Covered | DetectWholeGenomeDoubling_SixtyPercentElevated_ReturnsTrue |
| M9 | ✅ Covered | DetectWholeGenomeDoubling_ExactlyHalfElevated_ReturnsFalse |
| M10 | ✅ Covered | DetectWholeGenomeDoubling_FortyPercentElevated_ReturnsFalse |
| M11 | ✅ Covered | DetectWholeGenomeDoubling_AllBalancedDiploid_ReturnsFalse |
| M12 | ✅ Covered | DetectWholeGenomeDoubling_AllMajorElevated_ReturnsTrue |
| M13 | ✅ Covered | DetectWholeGenomeDoubling_EmptyOrInvalid_Throws |
| S1 | ✅ Covered | DetectWholeGenomeDoubling_LohSegments_CountAsElevated |
| S2 | ✅ Covered | EstimatePloidy_WithHomozygousDeletionSegment_IncludesZeros |
| C1 | ✅ Covered | EstimatePloidy_NearTriploidGenome_ExceedsAneuploidyDirection |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Per-segment total CN supplied as `AlleleSpecificSegment` (total = Major+Minor; length = End−Start) | EstimatePloidy, DetectWholeGenomeDoubling |
| 2 | WGD fraction denominator is the supplied segments' total length (no external chromosome-size table in scope) | DetectWholeGenomeDoubling |

---

## 7. Open Questions / Decisions

1. Registry lists `DetectWholeGenomeDoubling(ploidy)` (scalar). The authoritative facets-suite/Bielski WGD definition is the major-CN≥2 / >50%-of-genome rule, which requires segments, not a scalar ploidy. The canonical method therefore takes segments. Decision recorded; registry method-signature note updated in the algorithm doc.
</content>
