# Test Specification: ONCO-HRD-001

**Test Unit ID:** ONCO-HRD-001
**Area:** Oncology
**Algorithm:** Homologous Recombination Deficiency (HRD) composite genomic-scar score
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Telli ML et al. (2016), Clin Cancer Res 22(15):3764–3773 | 1 | https://pubmed.ncbi.nlm.nih.gov/26957554/ | 2026-06-14 |
| 2 | Stewart MD et al. (2022), Oncologist 27(3):167–174 (HRD review) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC8914493/ | 2026-06-14 |
| 3 | Birkbak NJ et al. (2012), Cancer Discov 2(4):366–375 (TAI) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC3806629/ | 2026-06-14 |
| 4 | Popova T et al. (2012), Cancer Res 72(21):5454–5462 (LST) | 1 | https://aacrjournals.org/cancerres/article/72/21/5454/576090/ | 2026-06-14 |
| 5 | oncoscanR `score_loh` reference impl (Abkevich LOH) | 3 | https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html | 2026-06-14 |
| 6 | scarHRD `calc.hrd.R` / `calc.ai_new.R` / `calc.lst.R` / `scar_score.R` (Sztupinszki 2018) | 3 | https://github.com/sztup/scarHRD | 2026-06-23 |
| 7 | UCSC Genome Browser cytoBand acen — hg38 + hg19 centromere coordinates | 5 | https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cytoBand | 2026-06-23 |
| 8 | NCBI GRC modeled centromeres (GRCh38 cross-verification) | 2 | https://www.ncbi.nlm.nih.gov/grc/human | 2026-06-23 |

### 1.2 Key Evidence Points

1. The combined HRD score is an **unweighted sum of LOH, TAI, and LST scores** — Telli 2016 (verbatim). Corroborated by Stewart 2022 ("gLOH + TAI + LST").
2. **HRD-high cutoff is ≥ 42** (inclusive) — Telli 2016 ("HRD score ≥42"); Stewart 2022 ("a cutoff of 42").
3. LOH component: LOH regions > 15 Mb and < whole chromosome (Abkevich; oncoscanR `score_loh`).
4. TAI component: allelic-imbalance (major≠minor) regions extending to a sub-telomere but not crossing the centromere (Birkbak 2012; scarHRD `calc.ai_new`: AI==2 at a chromosome end on one arm relative to the centromere → telomeric AI==1; ≥1 Mb segments only; single-segment imbalance is whole-chr AI==3, not telomeric).
5. LST component: chromosomal breaks between two adjacent ≥ 10 Mb regions after smoothing < 3 Mb regions, gap between the pair < 3 Mb, counted per arm (Popova 2012; scarHRD `calc.lst`).
6. Centromere coordinates (p-arm end = `chrominfo[i,2]`, q-arm start = `chrominfo[i,3]`) are the UCSC cytoBand acen bands, embedded for GRCh38 + GRCh37, cross-verified vs NCBI GRC modeled centromeres.

### 1.3 Documented Corner Cases

- Boundary: score exactly 42 → HRD-high; 41 → HRD-negative (Telli 2016, inclusive cutoff).
- Near-diploid / low-signal tumours produce a small sum; all-zero components → 0 → HRD-negative.
- Component counts are non-negative event counts; negative inputs are invalid.

### 1.4 Known Failure Modes / Pitfalls

1. Using a strict `> 42` instead of `≥ 42` would misclassify exactly-42 samples — Telli 2016 ("≥42").
2. Weighting any component would violate the "unweighted sum" definition — Telli 2016.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateHRDScore(int loh, int tai, int lst)` | OncologyAnalyzer | Canonical | Unweighted sum LOH+TAI+LST |
| `ClassifyHRDStatus(int score)` | OncologyAnalyzer | Canonical | ≥ 42 → HrdHigh |
| `DetectHRD(HrdComponents)` | OncologyAnalyzer | Canonical | End-to-end sum + classify |
| `DetectHRD(IEnumerable<AlleleSpecificSegment>, int tai, int lst)` | OncologyAnalyzer | Canonical | Derives LOH from segments (scarHRD `calc.hrd`); TAI/LST caller-supplied |
| `CalculateHrdTaiScore(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)` | OncologyAnalyzer | Canonical | Derives HRD-TAI (scarHRD `calc.ai_new`, even-ploidy path) |
| `CalculateHrdLstScore(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)` | OncologyAnalyzer | Canonical | Derives HRD-LST (scarHRD `calc.lst`, `chr.arm='no'`) |
| `DetectHRD(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)` | OncologyAnalyzer | Canonical | Derives ALL THREE components (LOH+TAI+LST) end-to-end (scarHRD `sum_HRD0`) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | HRD score = LOH + TAI + LST (unweighted) | Yes | Telli 2016 |
| INV-2 | Sum is order-independent (commutative) | Yes | Telli 2016 ("unweighted sum") |
| INV-3 | Status is HrdHigh iff score ≥ 42 | Yes | Telli 2016 |
| INV-4 | `DetectHRD(c).Score == CalculateHRDScore(c.Loh,c.Tai,c.Lst)` and Status matches `ClassifyHRDStatus(Score)` | Yes | composition contract |
| INV-5 | Component counts and score are ≥ 0 | Yes | counts are event counts (Birkbak/Popova/Abkevich) |
| INV-6 | `DetectHRD(segments,tai,lst).Components.Loh == DetectLOH(segments).Score` (LOH derived, not supplied) | Yes | scarHRD `calc.hrd`; ONCO-LOH-001 |
| INV-7 | TAI and LST are ≥ 0 event counts derived from segments; sex chromosomes contribute 0 to both | Yes | scarHRD `calc.ai_new` / `calc.lst` (autosome-only) |
| INV-8 | `DetectHRD(segments).Components == (DetectLOH(segments).Score, CalculateHrdTaiScore(segments), CalculateHrdLstScore(segments))` | Yes | scarHRD `sum_HRD0` |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Score sum | CalculateHRDScore(20,15,12) | 47 | Telli 2016 (unweighted sum) |
| M2 | Score sum 2 | CalculateHRDScore(5,4,3) | 12 | Telli 2016 |
| M3 | Boundary high | ClassifyHRDStatus(42) | HrdHigh | Telli 2016 ("≥42") |
| M4 | Below boundary | ClassifyHRDStatus(41) | HrdNegative | Telli 2016 |
| M5 | Above boundary | ClassifyHRDStatus(100) | HrdHigh | Telli 2016 |
| M6 | Boundary via sum | CalculateHRDScore(14,14,14)=42 → HRD-high | score 42, HrdHigh | Telli 2016 |
| M7 | Just-below via sum | CalculateHRDScore(14,13,14)=41 → HRD-negative | score 41, HrdNegative | Telli 2016 |
| M8 | End-to-end high | DetectHRD(20,15,12) | Score 47, HrdHigh, components preserved | Telli 2016 |
| M9 | End-to-end negative | DetectHRD(5,4,3) | Score 12, HrdNegative | Telli 2016 |
| M10 | Segment-driven LOH derivation | DetectHRD(LOH-dataset, tai=25, lst=16) derives LOH=1 | Loh=1, Score 42, HrdHigh | scarHRD `calc.hrd` (ONCO-LOH-001 dataset = 1); Telli 2016 sum |
| M11 | Two-path consistency | DetectHRD(segments,4,3) == DetectHRD(HrdComponents(DetectLOH(segments).Score,4,3)) | equal; Score 8, HrdNegative | INV-6 |
| M12 | No LOH segments | DetectHRD(balanced segments, tai=10, lst=5) | Loh=0, Score 15, HrdNegative | scarHRD `calc.hrd`; Telli 2016 |
| M13 | TAI both telomeric arms | chr1: p-terminal imbalance (end<cen start) + q-terminal imbalance (start>cen end) | TAI=2 | Birkbak 2012; scarHRD `calc.ai_new` |
| M14 | TAI interstitial | chr1 imbalanced interior segment, balanced ends | TAI=0 | scarHRD `calc.ai_new` (AI=2 interstitial) |
| M15 | TAI crossing centromere | first imbalanced segment whose end ≥ cen start | TAI=0 | Birkbak 2012 ("not cross the centromere") |
| M16 | TAI whole-chromosome | single imbalanced segment spanning chr1 | TAI=0 | scarHRD `calc.ai_new` (AI=3) |
| M17 | TAI sub-1 Mb dropped | < 1 Mb terminal imbalanced fragment + balanced rest | TAI=0 | scarHRD `min.size=1e6` |
| M18 | LST adjacent large pair | two adjacent ≥10 Mb p-arm segments, diff state, gap<3 Mb | LST=1 | Popova 2012; scarHRD `calc.lst` |
| M19 | LST one side <10 Mb | 40 Mb / 5 Mb / 65 Mb chain (middle ≥3 Mb, <10 Mb) | LST=0 | Popova 2012 ("each ≥10 Mb") |
| M20 | LST 3 Mb smoothing | 40 Mb / 2 Mb / 48 Mb (middle <3 Mb smoothed) | LST=1 | scarHRD smoothing while-loop |
| M21 | LST single segment | one segment on chr1 | LST=0 | scarHRD (`nrow<2 → next`) |
| M22 | LST q-arm transition | two adjacent ≥10 Mb q-arm segments past cen end | LST=1 | scarHRD `calc.lst` q.arm block |
| M23 | All-derived sum | DetectHRD(segments) derives LOH=0,TAI=2,LST=1 | Score 3, HrdNegative | scarHRD `sum_HRD0`; Telli 2016 |
| M24 | All-derived consistency | DetectHRD(segments) == DetectHRD(HrdComponents(DetectLOH, CalcTai, CalcLst)) | equal | INV-8 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Near-diploid zero | DetectHRD(0,0,0) | Score 0, HrdNegative | low-signal edge case |
| S2 | Negative component | CalculateHRDScore(-1,0,0) | ArgumentOutOfRangeException | counts ≥ 0 |
| S3 | Negative score | ClassifyHRDStatus(-1) | ArgumentOutOfRangeException | score ≥ 0 |
| S4 | Zero score classify | ClassifyHRDStatus(0) | HrdNegative | well-defined low signal |
| S5 | Null segments | DetectHRD(null, 0, 0) | ArgumentNullException | segment-driven overload guard |
| S6 | Negative TAI/LST (segment path) | DetectHRD(segments, -1, 0) / (segments, 0, -1) | ArgumentOutOfRangeException | counts ≥ 0 (Birkbak/Popova) |
| S7 | TAI only q-arm | chr1: balanced first, q-terminal imbalance only | TAI=1 | one telomeric side |
| S8 | TAI sex chromosome | chrX imbalanced terminal segments | TAI=0 | autosome-only table |
| S9 | TAI GRCh37 table | chr1 q-terminal at 130 Mb under hg19 (cen end 128.9 Mb) | TAI=1 | hg19 cytoBand acen |
| S10 | LST sex chromosome | chrX adjacent large segments | LST=0 | scarHRD excludes X/Y |
| S11 | TAI/LST null + empty | null → throw; empty → 0 | guards | input validation |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Order invariance | CalculateHRDScore over permutations of (20,15,12) | all equal 47 | INV-2 commutativity |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Original unit (2026-06-14): no prior HRD tests/implementation; M1–M9, S1–S4, C1 implemented in `OncologyAnalyzer_CalculateHRDScore_Tests.cs`.
- Segment-driven LOH derivation (2026-06-23): the `DetectHRD(segments, tai, lst)` overload adds M10–M12, S5, S6 to the same canonical fixture. LOH derivation reuses `DetectLOH` (ONCO-LOH-001, already scarHRD-verified).
- Segment-driven TAI + LST derivation (2026-06-23, this fix): adds `CalculateHrdTaiScore`, `CalculateHrdLstScore`, and the all-derived `DetectHRD(segments, ReferenceGenome)` overload, with M13–M24 and S7–S11 in the same canonical fixture. Centromere coordinates embedded from UCSC cytoBand acen (GRCh38 + GRCh37).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing → ✅ | new (2026-06-14) |
| M2 | ❌ Missing → ✅ | new (2026-06-14) |
| M3 | ❌ Missing → ✅ | new (2026-06-14) |
| M4 | ❌ Missing → ✅ | new (2026-06-14) |
| M5 | ❌ Missing → ✅ | new (2026-06-14) |
| M6 | ❌ Missing → ✅ | new (2026-06-14) |
| M7 | ❌ Missing → ✅ | new (2026-06-14) |
| M8 | ❌ Missing → ✅ | new (2026-06-14) |
| M9 | ❌ Missing → ✅ | new (2026-06-14) |
| M10 | ❌ Missing → ✅ | new (2026-06-23, segment-driven LOH) |
| M11 | ❌ Missing → ✅ | new (2026-06-23, two-path consistency) |
| M12 | ❌ Missing → ✅ | new (2026-06-23, no-LOH segments) |
| S1 | ❌ Missing → ✅ | new (2026-06-14) |
| S2 | ❌ Missing → ✅ | new (2026-06-14) |
| S3 | ❌ Missing → ✅ | new (2026-06-14) |
| S4 | ❌ Missing → ✅ | new (2026-06-14) |
| S5 | ❌ Missing → ✅ | new (2026-06-23, null segments) |
| S6 | ❌ Missing → ✅ | new (2026-06-23, negative TAI/LST) |
| C1 | ❌ Missing → ✅ | new (2026-06-14) |
| M13 | ❌ Missing → ✅ | new (2026-06-23, TAI both arms) |
| M14 | ❌ Missing → ✅ | new (2026-06-23, TAI interstitial) |
| M15 | ❌ Missing → ✅ | new (2026-06-23, TAI crossing centromere) |
| M16 | ❌ Missing → ✅ | new (2026-06-23, TAI whole-chromosome) |
| M17 | ❌ Missing → ✅ | new (2026-06-23, TAI sub-1 Mb dropped) |
| M18 | ❌ Missing → ✅ | new (2026-06-23, LST adjacent large pair) |
| M19 | ❌ Missing → ✅ | new (2026-06-23, LST one side <10 Mb) |
| M20 | ❌ Missing → ✅ | new (2026-06-23, LST 3 Mb smoothing) |
| M21 | ❌ Missing → ✅ | new (2026-06-23, LST single segment) |
| M22 | ❌ Missing → ✅ | new (2026-06-23, LST q-arm) |
| M23 | ❌ Missing → ✅ | new (2026-06-23, all-derived sum) |
| M24 | ❌ Missing → ✅ | new (2026-06-23, all-derived consistency) |
| S7 | ❌ Missing → ✅ | new (2026-06-23, TAI only q-arm) |
| S8 | ❌ Missing → ✅ | new (2026-06-23, TAI sex chr) |
| S9 | ❌ Missing → ✅ | new (2026-06-23, TAI GRCh37 table) |
| S10 | ❌ Missing → ✅ | new (2026-06-23, LST sex chr) |
| S11 | ❌ Missing → ✅ | new (2026-06-23, TAI/LST null+empty) |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateHRDScore_Tests.cs` — all HRD cases.
- **Remove:** nothing (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_CalculateHRDScore_Tests.cs` | Canonical HRD fixture | 42 |

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
| 13 | S4 | ❌ Missing | Implemented | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented | ✅ Done |
| 15 | M10 | ❌ Missing | Implemented | ✅ Done |
| 16 | M11 | ❌ Missing | Implemented | ✅ Done |
| 17 | M12 | ❌ Missing | Implemented | ✅ Done |
| 18 | S5 | ❌ Missing | Implemented | ✅ Done |
| 19 | S6 | ❌ Missing | Implemented | ✅ Done |
| 20 | M13 | ❌ Missing | Implemented | ✅ Done |
| 21 | M14 | ❌ Missing | Implemented | ✅ Done |
| 22 | M15 | ❌ Missing | Implemented | ✅ Done |
| 23 | M16 | ❌ Missing | Implemented | ✅ Done |
| 24 | M17 | ❌ Missing | Implemented | ✅ Done |
| 25 | M18 | ❌ Missing | Implemented | ✅ Done |
| 26 | M19 | ❌ Missing | Implemented | ✅ Done |
| 27 | M20 | ❌ Missing | Implemented | ✅ Done |
| 28 | M21 | ❌ Missing | Implemented | ✅ Done |
| 29 | M22 | ❌ Missing | Implemented | ✅ Done |
| 30 | M23 | ❌ Missing | Implemented | ✅ Done |
| 31 | M24 | ❌ Missing | Implemented | ✅ Done |
| 32 | S7 | ❌ Missing | Implemented | ✅ Done |
| 33 | S8 | ❌ Missing | Implemented | ✅ Done |
| 34 | S9 | ❌ Missing | Implemented | ✅ Done |
| 35 | S10 | ❌ Missing | Implemented | ✅ Done |
| 36 | S11 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 36
**✅ Done:** 36 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | CalculateHRDScore_StandardComponents_ReturnsUnweightedSum |
| M2 | ✅ | CalculateHRDScore_SmallComponents_ReturnsUnweightedSum |
| M3 | ✅ | ClassifyHRDStatus_ExactCutoff_ReturnsHrdHigh |
| M4 | ✅ | ClassifyHRDStatus_OneBelowCutoff_ReturnsHrdNegative |
| M5 | ✅ | ClassifyHRDStatus_WellAboveCutoff_ReturnsHrdHigh |
| M6 | ✅ | DetectHRD_ComponentsSummingTo42_IsHrdHigh |
| M7 | ✅ | DetectHRD_ComponentsSummingTo41_IsHrdNegative |
| M8 | ✅ | DetectHRD_HighComponents_ReturnsScoreAndHrdHigh |
| M9 | ✅ | DetectHRD_LowComponents_ReturnsScoreAndHrdNegative |
| S1 | ✅ | DetectHRD_NearDiploidZeroComponents_IsHrdNegative |
| S2 | ✅ | CalculateHRDScore_NegativeComponent_Throws |
| S3 | ✅ | ClassifyHRDStatus_NegativeScore_Throws |
| S4 | ✅ | ClassifyHRDStatus_ZeroScore_ReturnsHrdNegative |
| C1 | ✅ | CalculateHRDScore_ComponentOrderPermuted_ReturnsSameSum |
| M10 | ✅ | DetectHRD_FromSegments_DerivesLohAndSumsToScore |
| M11 | ✅ | DetectHRD_FromSegments_MatchesDetectLohPlusComponentsPath |
| M12 | ✅ | DetectHRD_FromSegmentsWithNoLoh_DerivesZeroLoh |
| S5 | ✅ | DetectHRD_FromSegments_NullSegments_Throws |
| S6 | ✅ | DetectHRD_FromSegments_NegativeTaiOrLst_Throws |
| M13 | ✅ | CalculateHrdTaiScore_BothTelomericArmsImbalanced_CountsTwo |
| M14 | ✅ | CalculateHrdTaiScore_InterstitialImbalance_NotCounted |
| M15 | ✅ | CalculateHrdTaiScore_FirstSegmentCrossingCentromere_NotCounted |
| M16 | ✅ | CalculateHrdTaiScore_SingleImbalancedSegment_WholeChromosomeNotCounted |
| M17 | ✅ | CalculateHrdTaiScore_SubMegabaseTerminalSegment_Dropped |
| M18 | ✅ | CalculateHrdLstScore_TwoAdjacentLargeArmSegments_CountsOne |
| M19 | ✅ | CalculateHrdLstScore_OneNeighbourBelow10Mb_NotCounted |
| M20 | ✅ | CalculateHrdLstScore_ShortSegmentSmoothed_ExposesTransition |
| M21 | ✅ | CalculateHrdLstScore_SingleSegment_NotCounted |
| M22 | ✅ | CalculateHrdLstScore_QArmTransition_CountsOne |
| M23 | ✅ | DetectHRD_AllDerivedFromSegments_SumsThreeComponents |
| M24 | ✅ | DetectHRD_AllDerived_MatchesStandaloneComponentDerivations |
| S7 | ✅ | CalculateHrdTaiScore_OnlyQArmTelomericImbalanced_CountsOne |
| S8 | ✅ | CalculateHrdTaiScore_SexChromosome_Excluded |
| S9 | ✅ | CalculateHrdTaiScore_GRCh37CentromereTable_Used |
| S10 | ✅ | CalculateHrdLstScore_SexChromosome_Excluded |
| S11 | ✅ | CalculateHrdTaiScore_Null/Empty + CalculateHrdLstScore_Null/Empty + DetectHRD_AllDerived_NullSegments |

All in-scope cases ✅ (36 of 36).

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Even-ploidy / standard allelic-imbalance path for TAI: AI present ⟺ major ≠ minor (scarHRD `calc.ai_new` default even/diploid path `seg[,7]==seg[,8]`). `AlleleSpecificSegment` lacks the ASCAT per-sample ploidy / aberrant-cell-fraction columns scarHRD uses to re-derive AI on odd-ploidy chromosomes, so that branch is not reproduced. This is the dominant path and matches Birkbak's "regions of allelic imbalance". | `CalculateHrdTaiScore` |

The prior caller-supplied-TAI/LST assumption is **RESOLVED**: the per-chromosome centromere coordinates are the UCSC cytoBand `acen` regions, retrieved as citable text and cross-verified vs the NCBI GRC modeled-centromere table, then embedded for GRCh38 + GRCh37. TAI and LST are now derived from segments.

---

## 7. Open Questions / Decisions

1. **Decision (2026-06-23):** All three HRD components — LOH, TAI, LST — are now derived end-to-end from allele-specific segments. `DetectHRD(segments)` (genome-parameterised, default GRCh38) computes `DetectLOH` (Abkevich 2012 / scarHRD `calc.hrd`), `CalculateHrdTaiScore` (Birkbak 2012 / scarHRD `calc.ai_new`), and `CalculateHrdLstScore` (Popova 2012 / scarHRD `calc.lst`), then sums them (`sum_HRD0`) and classifies at ≥42 (Telli 2016). The earlier blocker — scarHRD's binary `chrominfo` centromere table — is resolved by embedding the UCSC cytoBand `acen` coordinates (citable, cross-verified vs NCBI GRC). Caller-supplied TAI/LST remains available via the `DetectHRD(segments, tai, lst)` overload for externally computed components.
2. **Verification:** scarHRD's published bundled example (`test1.small.seqz`) yields HRD-LOH=1, TAI=2, LST=0, HRD-sum=3, but it is a Sequenza file requiring ASCAT ploidy/cellularity preprocessing absent from `AlleleSpecificSegment`, so it is recorded but not reproduced end-to-end. TAI/LST expected values are instead derived from the verbatim `calc.ai_new` / `calc.lst` logic and the embedded centromere coordinates (deterministic worked examples M13–M24).
