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

### 1.2 Key Evidence Points

1. The combined HRD score is an **unweighted sum of LOH, TAI, and LST scores** — Telli 2016 (verbatim). Corroborated by Stewart 2022 ("gLOH + TAI + LST").
2. **HRD-high cutoff is ≥ 42** (inclusive) — Telli 2016 ("HRD score ≥42"); Stewart 2022 ("a cutoff of 42").
3. LOH component: LOH regions > 15 Mb and < whole chromosome (Abkevich; oncoscanR `score_loh`).
4. TAI component: allelic-imbalance regions extending to a sub-telomere but not crossing the centromere (Birkbak 2012).
5. LST component: chromosomal breaks between adjacent ≥ 10 Mb regions after filtering < 3 Mb (Popova 2012).

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

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Near-diploid zero | DetectHRD(0,0,0) | Score 0, HrdNegative | low-signal edge case |
| S2 | Negative component | CalculateHRDScore(-1,0,0) | ArgumentOutOfRangeException | counts ≥ 0 |
| S3 | Negative score | ClassifyHRDStatus(-1) | ArgumentOutOfRangeException | score ≥ 0 |
| S4 | Zero score classify | ClassifyHRDStatus(0) | HrdNegative | well-defined low signal |
| S5 | Null segments | DetectHRD(null, 0, 0) | ArgumentNullException | segment-driven overload guard |
| S6 | Negative TAI/LST (segment path) | DetectHRD(segments, -1, 0) / (segments, 0, -1) | ArgumentOutOfRangeException | counts ≥ 0 (Birkbak/Popova) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Order invariance | CalculateHRDScore over permutations of (20,15,12) | all equal 47 | INV-2 commutativity |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Original unit (2026-06-14): no prior HRD tests/implementation; M1–M9, S1–S4, C1 implemented in `OncologyAnalyzer_CalculateHRDScore_Tests.cs`.
- Segment-driven LOH derivation (2026-06-23): the new `DetectHRD(segments, tai, lst)` overload adds M10–M12, S5, S6 to the same canonical fixture. LOH derivation reuses `DetectLOH` (ONCO-LOH-001, already scarHRD-verified).

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

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateHRDScore_Tests.cs` — all HRD cases.
- **Remove:** nothing (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_CalculateHRDScore_Tests.cs` | Canonical HRD fixture | 21 |

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

**Total items:** 19
**✅ Done:** 19 | **⛔ Blocked:** 0 | **Remaining:** 0

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

All in-scope cases ✅ (19 of 19).

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | TAI and LST remain caller-supplied. Their faithful derivation (scarHRD `calc.ai_new` / `calc.lst`) requires the exact per-build centromere/telomere `chrominfo` table, shipped only as binary `R/sysdata.rda` and not retrievable as a verifiable numeric table; TAI's telomeric classification and LST's p/q-arm split are sensitive to those coordinates. Per the conditional guard they are left caller-supplied rather than approximated. LOH IS now derived (source-backed). | `DetectHRD(segments,tai,lst)` |

---

## 7. Open Questions / Decisions

1. **Decision (2026-06-23):** HRD-LOH is now derived end-to-end from allele-specific segments via `DetectHRD(segments, tai, lst)` → `DetectLOH` (Abkevich 2012 / scarHRD `calc.hrd`, no centromere table needed). TAI and LST stay caller-supplied: their derivation depends on scarHRD's exact binary centromere `chrominfo` table, which could not be retrieved/cross-verified in this session (Evidence §scarHRD point 4). Deriving them from an unverified centromere table would not reproduce scarHRD, so per the conditional guard they are not approximated. The sum and the 42 cutoff remain source-verified (Telli 2016, Stewart 2022).
