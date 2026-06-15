# Validation Report: SV-CNV-001 — Read-Depth Copy Number Variation Detection

- **Validated:** 2026-06-15   **Area:** StructuralVar
- **Canonical method(s):** `StructuralVariantAnalyzer.DetectCNV(IReadOnlyList<int> depthData, int windowSize, double? referenceDepth, string chromosome)`; delegate `SegmentCopyNumber(IEnumerable<double> logRatios, string chromosome)`; private `LogRatioToCopyNumber`, `OverallMedianNonZero`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (one rounding defect found and fixed in-session)

## Stage A — Description

### Sources opened this session (retrieved, not trusted by label)

1. **Yoon S, Xuan Z, Makarov V, Ye K, Sebat J (2009)** — *Sensitive and accurate detection of copy number variants using read depth of coverage*, Genome Research 19(9):1586–1592 (open access PMC2752127). Fetched https://pmc.ncbi.nlm.nih.gov/articles/PMC2752127/ this session. Confirms verbatim:
   - "a linear relationship between coverage and copy number" (read depth ∝ copy number).
   - "RD was measured by counting the number of mapped reads in 100-bp windows, assigning each read only once by its start position." (windowed read counting).
   - GC-correction equation `r_i' = r_i × (m / m_GC)` with **m = overall median read count across all windows genome-wide** — the source for the default reference baseline (ASSUMPTION A1).
2. **CNVkit `cnvlib/call.py`** (etal/cnvkit, master). Fetched https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py this session. Confirms verbatim:
   - `_log2_ratio_to_absolute_pure`: `ncopies = ref_copies * 2**log2_ratio`, docstring `.. math :: n = r*2^v`.
   - Negative clamp `ncopies = max(0.0, ncopies)` (impure path) → CN ≥ 0.
   - `absolute_threshold` replaces NaN log2 with the neutral reference copy number.
   - **Final integer CN in `do_call` = `.round().astype("int")`** — i.e. NumPy `ndarray.round()`.
3. **CNVkit calling docs** https://cnvkit.readthedocs.io/en/stable/calling.html (fetched). Confirms verbatim: "a single-copy gain … has a copy ratio of 3/2. In log2 scale, this is log2(3/2) = 0.585, and a single-copy loss is log2(1/2) = -1.0."

### Formula check
- log2 ratio = `log2(windowMeanRD / referenceRD)` = log2(observed/reference) — matches CNVkit `_log2_ratio_to_absolute` docstring `log2_ratio = log2(ncopies/ploidy)`.
- Absolute CN = `round(ploidy · 2^log2)`, ploidy = 2 (diploid) — matches `n = r·2^v` with r = ref_copies = 2, and `.round()`.
- Anchors hand-verified against the sources: log2(2/2)=0→CN2; log2(1/2)=−1→CN1; log2(3/2)=0.5849625→CN3; log2(2)=1→CN4. All exact.

### Edge-case semantics
- Zero-depth window: log2(0)=−∞ undefined → no-call (excluded). Sourced (Yoon RD=read count; CNVkit no-call for unusable signal).
- NaN log2 → no-call (CNVkit replaces with neutral; here `SegmentCopyNumber` drops it / breaks the run). Acceptable, documented.
- CN ≥ 0 (CNVkit `max(0.0, ncopies)`).
- Trailing partial window dropped (fixed-size windows, Yoon).
- Default reference = overall median of non-zero window means (Yoon `m`, ASSUMPTION A1; overridable via parameter).

### Independent cross-check (hand + numpy)
| meanRD | log2(RD/100) | 2·2^log2 | CN (round-half-to-even, numpy) |
|--------|--------------|----------|--------------------------------|
| 100 | 0.000000 | 2.0 | 2 |
| 50 | −1.000000 | 1.0 | 1 |
| 150 | +0.584963 | 3.0 | 3 |
| 200 | +1.000000 | 4.0 | 4 |
| 25 | −2.000000 | **0.5** | **0** (numpy.round(0.5)=0) |
| 400 | +2.000000 | 8.0 | 8 |

All non-midpoint anchors agree with the TestSpec. **Stage A: PASS** — biology and maths trace to retrieved primary + reference-implementation sources.

## Stage B — Implementation

- **Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs:858` (`DetectCNV`), `:870` (`DetectCnvIterator`), `:937` (`SegmentCopyNumber` overload), `:1011` (`LogRatioToCopyNumber`), `:1025` (`OverallMedianNonZero`).
- **Formula realised correctly:** windowing, mean RD, `Math.Log2(meanRd/reference)`, `DiploidPloidy * Math.Pow(2, logRatio)`, non-negative clamp, zero-depth/no-call guard, overall-median baseline — all match the validated description.

### Defect found & fixed (FR-SV-CNV-001-01)
`LogRatioToCopyNumber` used `Math.Round(copies, MidpointRounding.AwayFromZero)`. CNVkit's `do_call` uses NumPy `ndarray.round()`, which is **round-half-to-even** (banker's): `numpy.round(0.5)=0`, `numpy.round(2.5)=2`. The two disagree at every exact half-integer copy number:
- meanRD 25 vs ref 100 → copies = 2·2^(−2) = **0.5**. Code returned CN **1**; CNVkit/NumPy returns **0**.
- copies 2.5 → code 3, reference 2; copies 1.5 → code 2, reference 2 (agrees); etc.

**Fix:** changed to `MidpointRounding.ToEven` (matches NumPy and also the sibling probe-based `CreateSegment`, which already used the `Math.Round` ToEven default). Comment updated to cite the NumPy round-half-to-even basis. No other caller depends on `LogRatioToCopyNumber`.

### Cross-verification table recomputed vs code (after fix)
S1 array means [25,50,100,200,400] vs ref 100 → CN sequence **[0,1,2,4,8]** — now matches the sourced numpy values exactly (test `DetectCNV_NonNegativeDepths_…`). Midpoint cases [0.5→0, 1.5→2, 2.5→2, 3.5→4] match (test `SegmentCopyNumber_HalfIntegerCopyNumber_RoundsHalfToEven`).

### Variant/delegate consistency
`DetectCNV` and `SegmentCopyNumber(IEnumerable<double>)` share `LogRatioToCopyNumber`; both now round-half-to-even. The probe-based `SegmentCopyNumber(probes,…)`/`CreateSegment` already used ToEven — consistent.

### Test quality audit (HARD gate)
- **Sourced, not code-echo:** M1–M4 assert exact sourced log2/CN anchors; would fail a wrong implementation. M3 asserts `Math.Log2(1.5)` = the exact value the source rounds to 0.585. PASS.
- **Weak-assertion defect (fixed):** original S1 asserted only `>=0` and non-decreasing — a monotone-but-wrong implementation (e.g. CN = window index) would have passed. **Strengthened** to assert the exact sourced sequence `[0,1,2,4,8]` (still retaining the INV-03/INV-04 property checks). This strengthened assertion is what exposed the rounding defect.
- **Added** midpoint round-half-to-even test (`[TestCase]` 0.5/1.5/2.5/3.5) to lock the NumPy-equivalent behaviour and prevent regression. log2 derived in-test as `Math.Log2(target/2)` so the reconstruction is exact.
- **Coverage:** every public method/overload exercised (`DetectCNV`, both `SegmentCopyNumber` overloads via the delegate path), all Stage-A branches (neutral/loss/gain/amp, windowing, median baseline, explicit baseline, zero-depth no-call, NaN no-call, empty, null, partial-window, all-zero, invalid window size, half-integer rounding). No skipped/ignored/weakened assertions remain.
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6493`; `dotnet build` 0 errors; no new warnings in the changed files.

**Stage B: PASS-WITH-NOTES** — code realises the validated formula; one rounding divergence from the CNVkit reference found and completely fixed in-session, tests strengthened to sourced exact values.

## Verdict & follow-ups
- **Stage A:** PASS. **Stage B:** PASS-WITH-NOTES. **End-state:** ✅ CLEAN (defect fully fixed; build + full suite green).
- **Defect logged:** FR-SV-CNV-001-01 (round-half-away-from-zero vs CNVkit/NumPy round-half-to-even at half-integer copy numbers) — fixed.
- No remaining limitations.
