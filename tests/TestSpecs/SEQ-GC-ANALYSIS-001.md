# Test Specification: SEQ-GC-ANALYSIS-001

**Test Unit ID:** SEQ-GC-ANALYSIS-001
**Area:** Composition
**Algorithm:** Comprehensive GC Analysis (`GcSkewCalculator.AnalyzeGcContent`)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wikipedia: GC-content (cites Brock/Madigan & Martinko 2003) | 4 | https://en.wikipedia.org/wiki/GC-content | 2026-06-14 |
| 2 | Wikipedia: GC skew (cites Lobry 1996, Grigoriev 1998) | 4 | https://en.wikipedia.org/wiki/GC_skew | 2026-06-14 |
| 3 | Biopython Bio.SeqUtils v1.84 (`GC_skew`, `gc_fraction`) | 3 | https://biopython.org/docs/1.84/api/Bio.SeqUtils.html | 2026-06-14 |
| 4 | Cuemath: Population Variance (formula + worked example) | 4 | https://www.cuemath.com/data/population-variance/ | 2026-06-14 |

### 1.2 Key Evidence Points

1. GC content % = (G+C)/(A+T+G+C)×100 — Source 1 (Brock via Wikipedia).
2. GC skew = (G−C)/(G+C), range [−1,+1]; G+C=0 → 0 — Sources 2, 3 (Lobry; Biopython GC_skew).
3. AT skew = (A−T)/(A+T) — sibling SEQ-ATSKEW-001 (Charneski 2011 / Lobry 1996).
4. Population variance σ² = Σ(xᵢ−μ)²/N; worked example {12,13,12,14,19} → 6.8 — Source 4.
5. Windowed metrics are computed over sliding windows ("multiple windows along the sequence") — Source 3.

### 1.3 Documented Corner Cases

- Window with no G/C → GC skew 0 (zero-division handled) — Source 3.
- Pure-G window → skew +1; pure-C window → skew −1 — Source 2.
- Sequence shorter than window → no full window → empty windowed output → window-derived variance 0 (implementation contract).

### 1.4 Known Failure Modes / Pitfalls

1. Using sample variance (÷N−1) instead of population variance (÷N) changes the value — Source 4.
2. Counting ambiguous bases in skew — must be ignored (only G,C) — Source 3.
3. Reporting GC content as fraction vs percentage — repository convention is percentage (×100) — Sources 1, 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `AnalyzeGcContent(DnaSequence, windowSize, stepSize)` | GcSkewCalculator | **Canonical** | Aggregates GC content, GC skew, AT skew, windowed profiles, population variances. |
| `AnalyzeGcContent(string, windowSize, stepSize)` | GcSkewCalculator | **Delegate** | New string overload; same core, smoke-tested for delegation + null/empty. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | OverallGcSkew ∈ [−1, +1] | Yes | Source 2 (range −1..+1) |
| INV-2 | OverallAtSkew ∈ [−1, +1] | Yes | sibling SEQ-ATSKEW-001 |
| INV-3 | OverallGcContent ∈ [0, 100] | Yes | Source 1 (percentage of total bases) |
| INV-4 | GcContentVariance ≥ 0 and GcSkewVariance ≥ 0 | Yes | Source 4 (variance is a sum of squares ÷N) |
| INV-5 | SequenceLength equals input length; windowed-list count = number of full windows = ⌊(n−w)/step⌋+1 (0 if n<w) | Yes | Source 3 (sliding windows) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | OverallGcContent | "GGGCCAT": G+C=5 of 7 | 71.42857142857143 | Source 1 |
| M2 | OverallGcSkew | "GGGCCAT": (3−2)/5 | 0.2 | Sources 2,3 |
| M3 | OverallAtSkew | "GGGCCAT": A=1,T=1 | 0.0 | SEQ-ATSKEW-001 |
| M4 | GcSkewVariance | "GGCC", w=2,step=2 → windows GG(+1),CC(−1) | 1.0 | Source 4 + windowed dataset |
| M5 | GcContentVariance | "GGCC", w=2,step=2 → 100,100 | 0.0 | Source 4 |
| M6 | Population variance correctness | windowed GC% over `AAAGGGCCCTTT` w=3,step=3 (0,100,100,0) | mean 50, var = (2500·4)/4 = 2500 | Source 4 |
| M7 | Windowed counts & boundaries | "ACGTACGTAC" w=4,step=2 → 4 windows; first WindowStart=0,WindowEnd=3,Position=2 | 4 windows; boundaries exact | Source 3, INV-5 |
| M8 | Pure-G overall skew bound | "GGGG" | OverallGcSkew = +1, GcContent=100 | Source 2 |
| M9 | No G/C overall skew | "ATATAT" | OverallGcSkew=0, GcContent=0, AtSkew=0 | Source 3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Sequence shorter than window | "ACGT", w=10 | empty windowed lists; variances 0; overall scalars still computed | windowing contract |
| S2 | Null DnaSequence | AnalyzeGcContent((DnaSequence)null) | ArgumentNullException | parity w/ siblings |
| S3 | Empty/null string overload | AnalyzeGcContent("") and (string)null | zero result, empty windows, length 0 | parity w/ siblings |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Overload equivalence | string vs DnaSequence on same seq | identical OverallGcContent/Skew/AtSkew/variances/window counts | delegation |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for `AnalyzeGcContent`. Sibling fixtures exist for the same class: `GcSkewCalculator_CalculateAtSkew_Tests.cs`, `GcSkewCalculator_PredictReplicationOrigin_Tests.cs`. No prior coverage of `AnalyzeGcContent` / `GcAnalysisResult`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M9 | ❌ Missing | new unit, no tests existed |
| S1–S3 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculator_AnalyzeGcContent_Tests.cs` — all M/S/C cases for `AnalyzeGcContent`.
- **Remove:** nothing (no prior tests for this method).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| GcSkewCalculator_AnalyzeGcContent_Tests.cs | canonical (this unit) | 13 |

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
| 10 | S1 | ❌ Missing | implemented | ✅ Done |
| 11 | S2 | ❌ Missing | implemented | ✅ Done |
| 12 | S3 | ❌ Missing | implemented | ✅ Done |
| 13 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact 71.42857142857143 within 1e-10 |
| M2 | ✅ | exact 0.2 |
| M3 | ✅ | exact 0.0 |
| M4 | ✅ | exact 1.0 (population variance of {+1,−1}) |
| M5 | ✅ | exact 0.0 |
| M6 | ✅ | exact 2500 (population variance of {0,100,100,0}) |
| M7 | ✅ | 4 windows; first window boundaries exact |
| M8 | ✅ | skew +1, GC% 100 |
| M9 | ✅ | skew 0, GC% 0, AT skew 0 |
| S1 | ✅ | empty windowed lists, variances 0, scalars computed |
| S2 | ✅ | ArgumentNullException |
| S3 | ✅ | zero result, empty windows, length 0 |
| C1 | ✅ | overload equivalence |

**Total in-scope cases:** 13 — ✅ count: 13.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | GC content reported as percentage (×100) per Brock convention, not Biopython [0,1] fraction | M1, M5, M6, M8, M9 |
| 2 | "Variability" = population variance Σ(x−μ)²/N (not sample n−1) | M4, M5, M6 |

---

## 7. Open Questions / Decisions

1. Decided: add a `string` overload of `AnalyzeGcContent` for API parity with the other `GcSkewCalculator` methods (all of which expose both `DnaSequence` and `string`). No change to the numerical contract.
2. None remaining on numerical behavior — all formulas are source-backed.
