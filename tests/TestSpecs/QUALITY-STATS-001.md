# Test Specification: QUALITY-STATS-001

**Test Unit ID:** QUALITY-STATS-001
**Area:** Quality
**Algorithm:** Quality Statistics (Phred summary statistics; Q20/Q30 percentages)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Ewing & Green (1998), Genome Research 8(3) (via Phred Wikipedia) | 1 (via 4) | https://doi.org/10.1101/gr.8.3.186 ; https://en.wikipedia.org/wiki/Phred_quality_score | 2026-06-13 |
| 2 | Illumina — Sequencing Quality Scores | 2 | https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html | 2026-06-13 |
| 3 | Newcastle Univ. ASK — Variance and Standard Deviation | 1 | https://www.mas.ncl.ac.uk/ask/numeracy-maths-statistics/statistics/descriptive-statistics/variance-and-standard-deviation.html | 2026-06-13 |
| 4 | Math is Fun — Median | 4 | https://www.mathsisfun.com/median.html | 2026-06-13 |
| 5 | Cock et al. (2010), NAR 38(6) | 1 | https://doi.org/10.1093/nar/gkp1137 | 2026-06-13 |

### 1.2 Key Evidence Points

1. Phred score relates to error probability by `Q = -10 log10 P`; Q30 = 1-in-1000 = 99.9% accuracy — Illumina; Phred Wikipedia (Ewing & Green 1998).
2. %≥Q30 is the percentage of bases with quality score ≥ 30 (inclusive); Q30 is the NGS benchmark metric — Illumina.
3. Population standard deviation σ = √((1/N) Σ(xᵢ−μ)²) (÷N, not ÷(N−1)) — Newcastle.
4. Median = middle of sorted list (odd); average of the two central values (even) — Math is Fun.
5. Phred+33 decode score = ord(c) − 33 (0–93); Phred+64 score = ord(c) − 64 (0–62) — Cock et al. (2010).

### 1.3 Documented Corner Cases

- Even-length input averages the two central order statistics for the median; odd-length takes the single middle (Math is Fun).
- Q20/Q30 thresholds are inclusive (≥) — a base exactly at the threshold is counted (Illumina table).
- Single observation: σ = 0; mean = median = min = max (descriptive-statistics identities).

### 1.4 Known Failure Modes / Pitfalls

1. Using sample std dev (÷(N−1)) instead of population (÷N) — Newcastle (population applies: the quality string is the complete observed set).
2. Exclusive `>` instead of inclusive `≥` for Q30/Q20 — Illumina table maps Q30 itself to 99.9% (counted).
3. Wrong median branch for even-length input (taking a single element rather than averaging the pair) — Math is Fun.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateStatistics(string, encoding)` | QualityScoreAnalyzer | Canonical | mean/median/min/max/stddev/Q20/Q30 over decoded scores |
| `CalculateQ30Percentage(string, encoding)` | QualityScoreAnalyzer | Canonical | new; = PercentAboveQ30 (%bases ≥ Q30) |
| `CalculateStatistics(IEnumerable<string>, encoding)` | QualityScoreAnalyzer | Delegate | multi-read aggregation; smoke only |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | MinQuality ≤ MeanQuality ≤ MaxQuality for any non-empty input | Yes | Order-statistic / mean bound (standard) |
| INV-2 | StandardDeviation ≥ 0; = 0 iff all scores equal | Yes | Newcastle (σ = √variance ≥ 0) |
| INV-3 | 0 ≤ PercentAboveQ30 ≤ PercentAboveQ20 ≤ 100 (Q30 ⊆ Q20) | Yes | Illumina (≥30 ⇒ ≥20) |
| INV-4 | `CalculateQ30Percentage(s)` == `CalculateStatistics(s).PercentAboveQ30` | Yes | Q30 = %bases ≥ Q30 (Illumina) |
| INV-5 | Statistics depend only on decoded scores, not on the encoding used to obtain them | Yes | Cock et al. (2010) score invariance |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Mean/min/max on "5?I" | Phred33 decode 20,30,40 | Mean=30.0, Min=20, Max=40 | Phred decode (Cock); arithmetic mean |
| M2 | Population std dev on "5?I" | σ=√(200/3) | 8.16496580927726 | Newcastle (÷N) |
| M3 | Median odd on "5?I" | n=3 middle | 30.0 | Math is Fun (odd) |
| M4 | Median even on "5II?" | sorted 20,30,40,40 → (30+40)/2 | 35.0 | Math is Fun (even) |
| M5 | %≥Q30 inclusive on "5?I" | 2 of 3 bases ≥30 | 66.66666666666667 | Illumina (≥30) |
| M6 | %≥Q20 inclusive on "5?I" | 3 of 3 ≥20 | 100.0 | Illumina (≥20) |
| M7 | `CalculateQ30Percentage` == PercentAboveQ30 | same input "5?I" | both 66.66666666666667 | INV-4 / Illumina |
| M8 | TotalBases / BasesAboveQ30 counts | "5?I" | TotalBases=3, BasesAboveQ30=2, BasesAboveQ20=3 | Illumina thresholds |
| M9 | Encoding invariance | "5?I"(P33) vs equivalent P64 chars decode same scores | identical Mean & PercentAboveQ30 | Cock et al. (2010) INV-5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Single base "I" | Q40 | Mean=Median=Min=Max=40, StdDev=0 | zero-spread identity |
| S2 | All ≥ Q30 | "?I" (30,40) | PercentAboveQ30=100.0 | upper percentage boundary |
| S3 | None ≥ Q30 | "5" (20) | PercentAboveQ30=0.0, PercentAboveQ20=100.0 | lower Q30 boundary, INV-3 |
| S4 | Q30 boundary inclusivity | "?" (exactly 30) | counted in Q30 (100.0) | Illumina ≥ semantics |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Multi-read delegate | two equal reads "5?I","5?I" | Mean=30.0, PercentAboveQ30≈66.667 | smoke for IEnumerable overload |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/QualityScoreAnalyzerTests.cs` — broad legacy fixture for the class; checked for `CalculateStatistics`/Q30 coverage.
- `tests/Seqeron/Seqeron.Genomics.Tests/QualityScoreAnalyzer_ParseQualityString_Tests.cs` — QUALITY-PHRED-001 (decode only; out of scope here).
- `FastqParser*` Q30 tests are for a different class (`FastqParser`), not this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M9 | ❌ Missing | no canonical evidence-based stats tests for QualityScoreAnalyzer.CalculateStatistics |
| S1–S4 | ❌ Missing | — |
| C1 | ❌ Missing | — |
| `CalculateQ30Percentage` | ❌ Missing | method does not yet exist |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/QualityScoreAnalyzer_CalculateStatistics_Tests.cs` — all M/S/C cases for this unit.
- **Remove:** nothing. Legacy `QualityScoreAnalyzerTests.cs` covers other methods (trim, mask, filter) and is left untouched; it has no exact-value CalculateStatistics tests that would duplicate this unit.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| QualityScoreAnalyzer_CalculateStatistics_Tests.cs | Canonical for QUALITY-STATS-001 | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented (.Within(1e-10)) | ✅ Done |
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
| 14 | C1 | ❌ Missing | Implemented (delegate smoke) | ✅ Done |
| 15 | Edge: empty | ❌ Missing | Implemented (zeroed result + Q30%=0) | ✅ Done |
| 16 | Edge: null | ❌ Missing | Implemented (zeroed result, no throw) | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | CalculateStatistics_Phred33Triplet_MeanMinMax |
| M2 | ✅ | CalculateStatistics_Phred33Triplet_PopulationStdDev |
| M3 | ✅ | CalculateStatistics_OddCount_MedianMiddle |
| M4 | ✅ | CalculateStatistics_EvenCount_MedianAveragesTwoCentral |
| M5 | ✅ | CalculateStatistics_Phred33Triplet_PercentAboveQ30 |
| M6 | ✅ | CalculateStatistics_Phred33Triplet_PercentAboveQ20 |
| M7 | ✅ | CalculateQ30Percentage_MatchesStatisticsPercentAboveQ30 |
| M8 | ✅ | CalculateStatistics_Phred33Triplet_Counts |
| M9 | ✅ | CalculateStatistics_EncodingInvariance_SameScoresSameStats |
| S1 | ✅ | CalculateStatistics_SingleBase_ZeroStdDev |
| S2 | ✅ | CalculateStatistics_AllAboveQ30_HundredPercent |
| S3 | ✅ | CalculateStatistics_NoneAboveQ30_ZeroPercent |
| S4 | ✅ | CalculateStatistics_ExactlyQ30_CountedInclusive |
| C1 | ✅ | CalculateStatistics_MultiReadDelegate_AggregatesScores |
| Edge: empty | ✅ | CalculateStatistics_Empty_ReturnsZeroedResult |
| Edge: null | ✅ | CalculateStatistics_Null_ReturnsZeroedResult |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Empty/null quality string → zeroed `QualityStatistics` (TotalBases=0) / Q30%=0 rather than throwing (non-correctness-affecting; no value invented) | Empty/null edge tests |

---

## 7. Open Questions / Decisions

1. None. Population std dev (÷N), inclusive Q20/Q30 (≥), and even/odd median rules are all source-backed; empty/null shape is the documented repository contract.
