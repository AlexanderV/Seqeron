# Test Specification: ONCO-HRD-001

**Test Unit ID:** ONCO-HRD-001
**Area:** Oncology
**Algorithm:** Homologous Recombination Deficiency (HRD) composite genomic-scar score
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

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

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | HRD score = LOH + TAI + LST (unweighted) | Yes | Telli 2016 |
| INV-2 | Sum is order-independent (commutative) | Yes | Telli 2016 ("unweighted sum") |
| INV-3 | Status is HrdHigh iff score ≥ 42 | Yes | Telli 2016 |
| INV-4 | `DetectHRD(c).Score == CalculateHRDScore(c.Loh,c.Tai,c.Lst)` and Status matches `ClassifyHRDStatus(Score)` | Yes | composition contract |
| INV-5 | Component counts and score are ≥ 0 | Yes | counts are event counts (Birkbak/Popova/Abkevich) |

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

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Near-diploid zero | DetectHRD(0,0,0) | Score 0, HrdNegative | low-signal edge case |
| S2 | Negative component | CalculateHRDScore(-1,0,0) | ArgumentOutOfRangeException | counts ≥ 0 |
| S3 | Negative score | ClassifyHRDStatus(-1) | ArgumentOutOfRangeException | score ≥ 0 |
| S4 | Zero score classify | ClassifyHRDStatus(0) | HrdNegative | well-defined low signal |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Order invariance | CalculateHRDScore over permutations of (20,15,12) | all equal 47 | INV-2 commutativity |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for HRD; no prior `CalculateHRDScore`/`ClassifyHRDStatus`/`DetectHRD` implementation. Searched `tests/Seqeron/Seqeron.Genomics.Tests/` (no `*HRD*` file) and `OncologyAnalyzer.cs` (no HRD members) — both absent before this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing → ✅ | new |
| M2 | ❌ Missing → ✅ | new |
| M3 | ❌ Missing → ✅ | new |
| M4 | ❌ Missing → ✅ | new |
| M5 | ❌ Missing → ✅ | new |
| M6 | ❌ Missing → ✅ | new |
| M7 | ❌ Missing → ✅ | new |
| M8 | ❌ Missing → ✅ | new |
| M9 | ❌ Missing → ✅ | new |
| S1 | ❌ Missing → ✅ | new |
| S2 | ❌ Missing → ✅ | new |
| S3 | ❌ Missing → ✅ | new |
| S4 | ❌ Missing → ✅ | new |
| C1 | ❌ Missing → ✅ | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateHRDScore_Tests.cs` — all HRD cases.
- **Remove:** nothing (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_CalculateHRDScore_Tests.cs` | Canonical HRD fixture | 14 |

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

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

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

All in-scope cases ✅ (14 of 14).

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Component counts (LOH/TAI/LST) are supplied as already-computed integers; computing them from raw segments is out of scope (ONCO-LOH-001/ONCO-CNA-001). API-shape only — sum and cutoff are source-backed. | §2 method signatures |

---

## 7. Open Questions / Decisions

1. **Decision:** Per the unit NOTE, this unit implements the retrievable composite-sum + 42-threshold classification taking the three component counts as input; per-segment component computation is deferred to dependent units. No open correctness questions remain — the sum and the 42 cutoff are both source-verified (Telli 2016, Stewart 2022).
