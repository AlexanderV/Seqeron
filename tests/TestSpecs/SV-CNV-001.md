# Test Specification: SV-CNV-001

**Test Unit ID:** SV-CNV-001
**Area:** StructuralVar
**Algorithm:** Read-Depth Copy Number Variation Detection
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Yoon S, Xuan Z, Makarov V, Ye K, Sebat J (2009). Sensitive and accurate detection of CNVs using read depth of coverage. Genome Research 19(9):1586–1592 | 1 | https://doi.org/10.1101/gr.092981.109 (PMC2752127) | 2026-06-13 |
| 2 | CNVkit `cnvlib/call.py` (etal/cnvkit, master); Talevich et al. 2016 PLoS Comput Biol 12(4):e1004873 | 3 | https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py · https://doi.org/10.1371/journal.pcbi.1004873 | 2026-06-13 |
| 3 | CNVkit "Calling copy number gains and losses" documentation | 3 | https://cnvkit.readthedocs.io/en/stable/calling.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. Read depth in non-overlapping fixed-size windows is proportional to copy number; "a linear relationship between coverage and copy number" — Yoon et al. 2009.
2. Read depth is "counting the number of mapped reads in 100-bp windows" — Yoon et al. 2009 (windowing).
3. log2 ratio = log2(observed copies / ploidy) = log2(observed depth / reference depth) — CNVkit `_log2_ratio_to_absolute` docstring.
4. Absolute copy number (pure/round method) = `ref_copies * 2^log2`, rounded to nearest integer; for diploid "the absolute copy number is calculated as 2 * 2^(log2 value)" — CNVkit `call.py` + docs.
5. Anchors: log2(1/2) = −1.0 ⇒ CN 1; log2(2/2) = 0 ⇒ CN 2; log2(3/2) = 0.585 ⇒ CN 3; log2(2) = 1.0 ⇒ CN 4 — CNVkit docs / call.py.
6. Normalisation baseline = overall median of windows (`m`) — Yoon et al. GC-correction equation `r_i' = r_i·m/m_GC`.
7. Copy number is physically ≥ 0 (CNVkit `max(0.0, ncopies)`); NaN/undefined-signal windows default to neutral / no-call.

### 1.3 Documented Corner Cases

- Zero-depth window ⇒ log2(0) undefined (−∞) ⇒ no-call, not a finite call (Yoon RD = read count; CNVkit treats unusable signal as neutral/no-call).
- NaN log2 ⇒ replaced with neutral copy number (CNVkit `absolute_threshold`).
- Copy number clamped to ≥ 0 (CNVkit).

### 1.4 Known Failure Modes / Pitfalls

1. Taking log2 of a zero-RD window yields −∞ / undefined — must be guarded (Yoon RD definition; CNVkit no-call) .
2. Choosing the wrong reference baseline shifts the log2=0 anchor and miscalls every window — baseline = overall median per Yoon et al.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DetectCNV(depthData, window, ...)` | StructuralVariantAnalyzer | Canonical | Windowed RD → log2 ratio → integer copy number → per-window CopyNumberSegment |
| `SegmentCopyNumber(logRatios, ...)` | StructuralVariantAnalyzer | Delegate | Converts a sequence of per-window log2 ratios into copy-number calls and merges adjacent equal-CN runs; thin wrapper over the same log2→CN conversion used by `DetectCNV` |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | A window whose mean RD equals the reference RD has log2 ratio 0 and copy number = ploidy (2). | Yes | CNVkit CN = 2·2^0 = 2; log2(2/2)=0 |
| INV-02 | Copy number = round(ploidy · 2^log2) where log2 = log2(windowMeanRD / referenceRD). | Yes | CNVkit `_log2_ratio_to_absolute_pure`; diploid docs |
| INV-03 | Copy number is a non-negative integer (≥ 0). | Yes | CNVkit `max(0.0, ncopies)` |
| INV-04 | Copy number is monotonically non-decreasing in window mean RD (more depth ⇒ ≥ copy number). | Yes | CN = 2·2^log2 strictly increasing in log2 (Yoon linear coverage↔CN) |
| INV-05 | Depth data is partitioned into non-overlapping windows of the given size; a trailing partial window is dropped (windows are fixed-size read-count windows). | Yes | Yoon et al. (fixed-size windows) |
| INV-06 | A zero-depth window produces a no-call (excluded), not a finite copy-number call. | Yes | Yoon RD=0 ⇒ undefined ratio; CNVkit no-call |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Neutral window | Window mean RD = reference RD (100) | log2 = 0.0; CopyNumber = 2; not flagged as CNV | CNVkit CN=2·2^0=2 (Src 2/3); INV-01/INV-02 |
| M2 | Single-copy loss | Window mean RD = 50, reference 100 | log2 = −1.0; CopyNumber = 1 (Deletion) | CNVkit "single-copy loss is log2(1/2) = -1.0"; round(2·2^−1)=1 |
| M3 | Single-copy gain | Window mean RD = 150, reference 100 | log2 ≈ 0.585; CopyNumber = 3 (Duplication) | CNVkit "single-copy gain ... log2(3/2)=0.585"; round(2·2^0.585)=3 |
| M4 | Amplification | Window mean RD = 200, reference 100 | log2 = 1.0; CopyNumber = 4 (Duplication) | CNVkit diploid CN=2·2^1=4 |
| M5 | Windowing | 8 positions, window=4, two windows | Two segments; each ProbeCount = 4; means over each window | Yoon et al. windowed read counts; INV-05 |
| M6 | Reference = overall median | Mixed-depth windows, no baseline given | Reference RD = overall median of window means; the median-depth window is CN 2 | Yoon overall-median baseline `m`; INV-01 (ASSUMPTION A1) |
| M7 | Zero-depth window no-call | A window of all-zero depth | Excluded from output (no finite call) | Yoon RD=0 undefined; CNVkit no-call; INV-06 |
| M8 | SegmentCopyNumber log2→CN + merge | log2 ratios [0,0,−1,−1] | Two segments: CN 2 (len 2) then CN 1 (len 2) | CNVkit CN=round(2·2^log2); same conversion as DetectCNV |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Non-negativity & monotonicity (property) | Random non-negative depth arrays | All CopyNumber ≥ 0; CN non-decreasing as window mean RD increases (sorted) | INV-03/INV-04 |
| S2 | Explicit baseline override | Supply referenceDepth=50; window mean 50 | log2 = 0; CN = 2 | Reference-against-profile (CNVkit); overrides ASSUMPTION A1 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty input | depthData empty | Empty result | Trivial defined behaviour |
| C2 | Null input | depthData null | ArgumentNullException | Input-validation contract |
| C3 | Window larger than data | window > length | Empty result (no full window) | INV-05 (partial window dropped) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `StructuralVariantAnalyzer.cs` (Annotation project) contains `SegmentCopyNumber(probes,...)` (probe-based, different signature), `IdentifyCNVs(segments,...)`, `CreateSegment(...)` — none is the Registry-canonical `DetectCNV(depthData, window)`, which does not exist.
- No test file exercises `DetectCNV`. Sibling unit tests: `StructuralVariantAnalyzer_DetectSVs_Tests.cs`, `StructuralVariantAnalyzer_FindBreakpoints_Tests.cs`. Legacy `StructuralVariantAnalyzerTests.cs` exists but does not cover read-depth CNV calling.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 Neutral window | ✅ Covered | Implemented Phase 7 |
| M2 Single-copy loss | ✅ Covered | Implemented Phase 7 |
| M3 Single-copy gain | ✅ Covered | Implemented Phase 7 |
| M4 Amplification | ✅ Covered | Implemented Phase 7 |
| M5 Windowing | ✅ Covered | Implemented Phase 7 |
| M6 Reference = overall median | ✅ Covered | Implemented Phase 7 |
| M7 Zero-depth no-call | ✅ Covered | Implemented Phase 7 |
| M8 SegmentCopyNumber merge | ✅ Covered | Implemented Phase 7 |
| S1 Non-negativity/monotonicity (property) | ✅ Covered | Implemented Phase 7 |
| S2 Explicit baseline override | ✅ Covered | Implemented Phase 7 |
| C1 Empty input | ✅ Covered | Implemented Phase 7 |
| C2 Null input | ✅ Covered | Implemented Phase 7 |
| C3 Window larger than data | ✅ Covered | Implemented Phase 7 |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/StructuralVariantAnalyzer_DetectCNV_Tests.cs` — all SV-CNV-001 cases (`DetectCNV`, `SegmentCopyNumber` smoke).
- **Remove:** nothing — no pre-existing `DetectCNV` tests exist; legacy `StructuralVariantAnalyzerTests.cs` is untouched (out of scope).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `StructuralVariantAnalyzer_DetectCNV_Tests.cs` | Canonical SV-CNV-001 fixture | 17 (13 planned cases + 4 branch-closure tests: invalid window, all-zero/no-baseline, NaN no-call merge) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented neutral-window test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented loss test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented gain test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented amplification test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented windowing test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented overall-median baseline test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented zero-depth no-call test | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented SegmentCopyNumber merge test | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented non-negativity + monotonicity property test | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented explicit-baseline test | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented empty-input test | ✅ Done |
| 12 | C2 | ❌ Missing | Implemented null-input test | ✅ Done |
| 13 | C3 | ❌ Missing | Implemented window>length test | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | `DetectCNV_NeutralWindow_CopyNumberTwo` |
| M2 | ✅ | `DetectCNV_HalfDepthWindow_CopyNumberOneDeletion` |
| M3 | ✅ | `DetectCNV_OnePointFiveDepthWindow_CopyNumberThreeDuplication` |
| M4 | ✅ | `DetectCNV_DoubleDepthWindow_CopyNumberFour` |
| M5 | ✅ | `DetectCNV_EightPositionsWindowFour_TwoSegments` |
| M6 | ✅ | `DetectCNV_NoBaseline_UsesOverallMedianReference` |
| M7 | ✅ | `DetectCNV_ZeroDepthWindow_ExcludedAsNoCall` |
| M8 | ✅ | `SegmentCopyNumber_LogRatios_MergesAdjacentEqualCopyNumber` |
| S1 | ✅ | `DetectCNV_NonNegativeDepths_CopyNumbersNonNegativeAndMonotone` |
| S2 | ✅ | `DetectCNV_ExplicitBaseline_OverridesMedian` |
| C1 | ✅ | `DetectCNV_EmptyInput_ReturnsEmpty` |
| C2 | ✅ | `DetectCNV_NullInput_Throws` |
| C3 | ✅ | `DetectCNV_WindowLargerThanData_ReturnsEmpty` |

**✅ count:** 13 of 13 in-scope cases.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | Reference (diploid baseline) RD defaults to the overall median of non-zero window means (Yoon overall-median `m`). | M6, INV-01 |
| A2 | Diploid ploidy (2) is the copy-number baseline (CNVkit diploid conversion). | M1–M4, INV-01/INV-02 |

---

## 7. Open Questions / Decisions

1. None. Both assumptions are source-anchored (Yoon overall-median baseline; CNVkit diploid `CN = 2·2^log2`) and an explicit baseline parameter is provided to override A1.
