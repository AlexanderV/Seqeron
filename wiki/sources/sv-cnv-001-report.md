---
type: source
title: "Validation report: SV-CNV-001 (read-depth CNV — windowed depth → log2 ratio → integer copy number)"
tags: [validation, structural-variant, governance]
doc_path: docs/Validation/reports/SV-CNV-001.md
sources:
  - docs/Validation/reports/SV-CNV-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: SV-CNV-001

The two-stage **validation write-up** for test unit **SV-CNV-001** (Read-Depth Copy
Number Variation Detection — windowed read depth → log2 ratio → integer copy number →
del/dup), validated 2026-06-15 (Area: StructuralVar). This is the *report* artifact that
feeds one row of the [[validation-ledger]]; it records the validator's **verdict** on both
the algorithm description and the shipped code. The two-stage methodology is the
[[validation-protocol]]; the copy-number algorithm itself is summarized in
[[read-depth-cnv-segmentation]]. Distinct from the pre-implementation
[[sv-cnv-001-evidence]] artifact.

## Verdict

**Stage A: PASS · Stage B: PASS-WITH-NOTES · State: ✅ CLEAN.** One rounding defect
(FR-SV-CNV-001-01) found and fully fixed in-session; tests strengthened to the sourced
exact values. Full unfiltered suite **6493 passed / 0 failed**, `dotnet build` 0 errors,
no new warnings in the changed files. No remaining limitations.

Canonical methods: `StructuralVariantAnalyzer.DetectCNV(IReadOnlyList<int> depthData,
int windowSize, double? referenceDepth, string chromosome)`; delegate
`SegmentCopyNumber(IEnumerable<double> logRatios, string chromosome)`; private
`LogRatioToCopyNumber`, `OverallMedianNonZero`.

## Stage A — description (algorithm faithfulness)

- Sources retrieved and read this session (not trusted by label): **Yoon et al. 2009**
  (Genome Research 19(9), PMC2752127) — the RD∝CN linear relationship, 100-bp windowed
  read counting (each read assigned once by start position), and the GC-correction
  equation `r_i' = r_i × (m / m_GC)` with `m` = overall median across all windows (the
  source for the default baseline, ASSUMPTION A1); **CNVkit `cnvlib/call.py`** — the
  `n = r·2^v` log2→absolute conversion, `max(0.0, ncopies)` non-negative clamp, NaN→neutral
  replacement, and final integer CN via NumPy `ndarray.round()`; **CNVkit calling docs** —
  the log2(3/2)=0.585 single-copy-gain / log2(1/2)=−1.0 single-copy-loss anchors.
- Formula confirmed: `log2 ratio = log2(windowMeanRD / referenceRD)`; absolute
  `CN = round(ploidy · 2^log2)` with ploidy = 2 (diploid). Anchors hand-verified: log2(2/2)=0→CN2,
  log2(1/2)=−1→CN1, log2(3/2)=0.5849625→CN3, log2(2)=1→CN4 — all exact.
- Edge cases defined and sourced: zero-depth window → log2(0)=−∞ → no-call (excluded);
  NaN log2 → no-call; CN clamped ≥ 0; trailing partial window dropped; default reference =
  overall median of non-zero window means (overridable).
- **Stage A: PASS** — biology and maths trace to retrieved primary + reference-implementation
  sources.

## Stage B — implementation (code review + cross-check)

- Code path: `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs`
  — `DetectCNV` (:858), `DetectCnvIterator` (:870), `SegmentCopyNumber` overload (:937),
  `LogRatioToCopyNumber` (:1011), `OverallMedianNonZero` (:1025). Windowing, mean RD,
  `Math.Log2`, `DiploidPloidy * Math.Pow(2, logRatio)`, non-negative clamp, no-call guard,
  and overall-median baseline all match the validated description.
- **Defect FR-SV-CNV-001-01 (found & fixed):** `LogRatioToCopyNumber` used
  `Math.Round(copies, MidpointRounding.AwayFromZero)`, but CNVkit's `do_call` uses NumPy
  `ndarray.round()` = **round-half-to-even** (banker's). They disagree at every exact
  half-integer CN — e.g. meanRD 25 vs ref 100 → copies 0.5, code returned CN 1 while
  CNVkit/NumPy returns 0. **Fix:** changed to `MidpointRounding.ToEven` (also matching the
  sibling probe-based `CreateSegment`); comment updated to cite the NumPy basis. No other
  caller depends on `LogRatioToCopyNumber`.
- Cross-verification after fix: array means [25,50,100,200,400] vs ref 100 → CN [0,1,2,4,8]
  (matches sourced numpy values); midpoints [0.5→0, 1.5→2, 2.5→2, 3.5→4] match. Both
  `SegmentCopyNumber` overloads and `CreateSegment` now round-half-to-even consistently.
- Test-quality audit (HARD gate) PASS: a weak S1 assertion (only `>=0` + non-decreasing,
  which a monotone-but-wrong impl would pass) was **strengthened** to assert the exact
  sourced sequence `[0,1,2,4,8]` — this is what exposed the rounding defect. Added a midpoint
  round-half-to-even `[TestCase]` (0.5/1.5/2.5/3.5) with in-test-derived log2 to lock the
  NumPy-equivalent behaviour. Every public method/overload and every Stage-A branch exercised;
  no skipped/weakened assertions remain.
- **Stage B: PASS-WITH-NOTES** — code realises the validated formula; one rounding divergence
  from the CNVkit reference found and completely fixed in-session.

## Findings

- **State ✅ CLEAN.** Defect **FR-SV-CNV-001-01** (round-half-away-from-zero vs
  CNVkit/NumPy round-half-to-even at half-integer copy numbers) logged and **fixed**.
- No remaining limitations; build + full suite green (6493 / 0).

See the report at `docs/Validation/reports/SV-CNV-001.md`.
