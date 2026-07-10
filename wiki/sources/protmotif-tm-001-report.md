---
type: source
title: "Validation report: PROTMOTIF-TM-001 (transmembrane helix prediction — Kyte-Doolittle hydropathy sliding window)"
tags: [validation, protein, governance]
doc_path: docs/Validation/reports/PROTMOTIF-TM-001.md
sources:
  - docs/Validation/reports/PROTMOTIF-TM-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: PROTMOTIF-TM-001

The two-stage **validation write-up** for test unit **PROTMOTIF-TM-001** (Transmembrane
Helix Prediction via the Kyte-Doolittle hydropathy sliding window), validated 2026-06-16
in the ProteinMotif area. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's **verdict** on both the algorithm
description and the shipped code. The method itself is summarized in
[[transmembrane-helix-prediction]]; the two-stage methodology is the
[[validation-protocol]]. Distinct from the pre-implementation
[[protmotif-tm-001-evidence]] artifact.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS (after in-session fix) · State: ✅ CLEAN.**
Full unfiltered suite **6579 passed / 0 failed**, `dotnet build` 0 errors, changed files
warning-free (the 4 NUnit2007 warnings are pre-existing in an unrelated test file). One
defect (an off-by-one segment-End coordinate) was found and completely fixed in-session,
with code, tests, and doc/Evidence/TestSpec all corrected to the sourced value.

Canonical method: `ProteinMotifFinder.PredictTransmembraneHelices(string, int windowSize=19,
double threshold=1.6)` (private helper `CalculateHydropathyProfile`). See report at
`docs/Validation/reports/PROTMOTIF-TM-001.md`.

## Stage A — description (algorithm faithfulness)

- **Parameters and scale fully source-confirmed.** Davidson College (Kyte-Doolittle
  background) and QIAGEN CLC (protein hydrophobicity) confirm window size **19** and
  threshold **1.6** ("peaks with scores greater than 1.6 using a window size of 19"). The
  20-value Kyte-Doolittle scale (I +4.5 … R −4.5) matches the implementation's
  `HydropathyScale` dictionary exactly, per a WebSearch of the full table.
- **Formula:** the profile point is the arithmetic mean of the window's per-residue values,
  `P(i) = (1/w)·Σ h(s[j])` — matches Davidson ("average of all the hydrophobicity scores in
  that window") and Biopython `protein_scale` with edge weight 1.0.
- **Edge cases** are sourced and consistent: sequences shorter than `w` → no segment;
  non-standard residues (X/B/Z/*) have no scale value and are excluded from the mean;
  `windowSize ≤ 0` / null / empty → empty.
- **Independent cross-check:** the validator recomputed the full 22-point profile in Python
  for `D×10 L×20 D×10` (40 res, w=19, T=1.6): above-threshold profile indices 5..16; the
  union of passing windows covers residues **5..34**, peak 3.8.
- **NOTE (the one non-source-prescribed item):** the segment-boundary mapping is not uniquely
  prescribed (score-at-first-residue per Davidson vs score-at-midpoint per QIAGEN). Stage A
  fixes the correct sourced End = `lastProfileIndex + windowSize − 1` (the last residue
  covered by any passing window), correcting the prior `+ windowSize` in the description,
  which named a residue outside every passing window.

## Stage B — implementation (code review + cross-check)

- **Code path:** `ProteinMotifFinder.cs:687-786`. Scale dictionary (751-757) matches the
  sourced KD table; profile (764-786) is the arithmetic mean over scored residues with
  non-standard residues excluded; scan (704-744) opens a run at the first `P(i) ≥ T` and
  closes at the first `P(i) < T`.
- **Defect found and FIXED (off-by-one End).** The original close logic reported
  `end = lastPassingIndex + windowSize` as the 0-based inclusive End; for M1 this yields 35,
  but the last residue covered by any passing window is 34 (last window starts at 16, covers
  16..34). Root cause: `end` was used as a half-open span boundary but reported as if
  inclusive without subtracting 1. Fix reports the last covered residue
  `lastPassingIndex + windowSize − 1` (clamped to length−1), preserving the inclusive-span
  filter.
- **Cross-verification table (recomputed vs fixed code, all PASS):** M1 `D×10 L×20 D×10` →
  (5, **34**, 3.8); M2 `D×40` → empty; M3 `L×19` → (0,18,3.8); M4 I/V/R×19 →
  (0,18,4.5)/(0,18,4.2)/none; M5–M7 null/""/`L×18` → empty; S1 `L×9 + X + L×9` → score 3.8
  (X excluded); S2 `D×10 A×20 D×10, T=2.0` → empty; S4 `windowSize=0` → empty.
- **Variant/delegate consistency:** single public method; MCP wrapper
  `AnalysisTools.PredictTransmembraneHelices` forwards; no `*Fast` variant. The same-file
  disorder predictor uses a different (midpoint) convention and is out of scope.
- **Test-quality audit (HARD gate) PASS:** M1/S3 Ends were code-echoes of the buggy
  convention (35) → rewritten to the sourced 34. Exact-equality assertions
  (`.Within(1e-10)`); the property test was tightened from `End ≤ length` to `End < length`
  + `Start ≤ End` (INV-02). All 12 TestSpec cases present.

## Findings

- **FIXED** — off-by-one End coordinate (reported a residue outside every passing window);
  code + tests + doc/Evidence/TestSpec corrected to `lastProfileIndex + windowSize − 1`.
- **Stage A NOTE** — the segment-end coordinate is the only non-source-prescribed item; now
  defined consistently across code/doc/Evidence/TestSpec as the last residue covered by a
  passing window. Parameters (window 19, threshold 1.6), the KD scale (20 values), and the
  arithmetic-mean profile are all fully and independently source-confirmed.
- **End-state ✅ CLEAN** — full suite 6579/0, build 0 errors.
