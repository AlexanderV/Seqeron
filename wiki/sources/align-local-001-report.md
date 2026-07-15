---
type: source
title: "Validation report: ALIGN-LOCAL-001 (Smith–Waterman local alignment)"
tags: [validation, alignment, governance]
doc_path: docs/Validation/reports/ALIGN-LOCAL-001.md
sources:
  - docs/Validation/reports/ALIGN-LOCAL-001.md
source_commit: b7e2c1eeb773db02af22541c7e80e6eb7019780c
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ALIGN-LOCAL-001

The two-stage **validation write-up** for test unit **ALIGN-LOCAL-001** (Local
Alignment, Smith–Waterman), validated 2026-06-24 as a re-validation in fresh context.
This is the *report* artifact that feeds one row of the [[validation-ledger]]; it records
the validator's **verdict** on both the algorithm description and the shipped code. The
algorithm itself is summarized in [[local-alignment-smith-waterman]]; the two-stage
methodology is the [[validation-protocol]].

## Verdict

**Stage A: PASS · Stage B: PASS · State: ✅ CLEAN.** No defects found, no code changes.
`SequenceAligner_LocalAlign_Tests` 7/7 green; build green. The local core/traceback path is
logically identical to the previously-validated version (last functionally touched by the
`e19a8a02` "O(n) traceback" perf work).

## Stage A — description (algorithm faithfulness)

- Canonical method: `SequenceAligner.LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)`
  (plus a `string` delegate) → `LocalAlignCore` → `TracebackLocal`.
- Checked against live Wikipedia "Smith–Waterman algorithm" + Smith & Waterman (1981):
  the **zero-floor** linear-gap recurrence `H(i,j)=max(0, diag+s, up−W₁, left−W₁)`; first
  row/column initialized to 0 (no end-gap penalty); optimum = **max cell anywhere** (not the
  corner); traceback from that max cell, **stopping at the first 0 cell**.
- Worked example `TGTTACGG` vs `GGTTGACTA` (+3/−3/−2) → **score 13**, alignment
  `GTT-AC` / `GTTGAC`. The validator **hand-recomputed the full DP matrix** cell-for-cell
  (max 13 at H[6,7]) and the swapped case (gap moves to seq2 via the "up" branch). No divergences.

## Stage B — implementation (code review + cross-check)

- Code path in `SequenceAligner.cs`: `LocalAlignCore` (L315–348) + `TracebackLocal` (L350–398);
  default `ScoringMatrix.SimpleDna`; linear gap supplied as the negative `GapExtend`.
- Evidence the formula is realised: zero floor at L335 `Math.Max(0, …)`; max-cell optimum
  tracked at L337–342 (does **not** take F(m,n)); traceback-to-zero via
  `while (i>0 && j>0 && score[i,j] > 0)` (L359) with diag→up→left tie priority.
- Cross-verification table recomputed vs code, **7/7 pass** — M1 Wikipedia (13), M2 swapped
  (13, gap on seq2), M3 string overload, M4 empty→`AlignmentResult.Empty`, M5 null→
  `ArgumentNullException`, S1 identical (24, full match), S2 all-mismatch (0, empty).
- `GapOpen` is intentionally unused (no affine term); string overload uppercases then routes
  to the same core; M3 confirms delegate parity.

## Findings

- **None.** No defects, no code changes required. State ✅ CLEAN.
