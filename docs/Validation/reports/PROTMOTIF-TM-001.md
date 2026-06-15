# Validation Report: PROTMOTIF-TM-001 — Transmembrane Helix Prediction (Kyte-Doolittle hydropathy)

- **Validated:** 2026-06-16   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.PredictTransmembraneHelices(string, int windowSize=19, double threshold=1.6)` (private helper `CalculateHydropathyProfile`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS (after in-session fix of an off-by-one End coordinate)

## Stage A — Description

### Sources opened & what they confirm (retrieved this session, 2026-06-16)
- **Davidson College — Kyte-Doolittle background** (https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm), WebFetch: *"a window size of 19 is needed"*; *"Transmembrane regions are identified by peaks with scores greater than 1.6 using a window size of 19"*; window value = arithmetic mean of the scores in the window, sliding one residue at a time; the window average is *"assign[ed] … to the first amino acid in the window."* Confirms window 19, threshold 1.6, arithmetic-mean profile.
- **QIAGEN CLC — Protein hydrophobicity** (https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/2200/index.php?manual=BE_Protein_hydrophobicity.html), WebFetch: *"Large window sizes of 19-21 are well suited for finding transmembrane domains if the values calculated are above 1.6"*; the window average is *"plotted … at the central position of the window."* Confirms parameters; documents the **alternative midpoint plotting convention**.
- **Kyte-Doolittle scale (one-letter values)** — WebSearch returned the full table: I +4.5, V +4.2, L +3.8, F +2.8, C +2.5, M +1.9, A +1.8, G −0.4, T −0.7, S −0.8, W −0.9, Y −1.3, P −1.6, H −3.2, E −3.5, Q −3.5, D −3.5, N −3.5, K −3.9, R −4.5. **Matches the implementation's `HydropathyScale` dictionary exactly** (all 20 residues).
- **TM α-helix length** — a single helix needs ≈18–21 residues (≈1.5 Å rise/residue) to cross the ≈30 Å bilayer; justifies the 19-residue window and a minimum-span filter equal to the window width.

### Formula check
`P(i) = (1/w) · Σ_{j=i}^{i+w-1} h(s[j])` — arithmetic mean over the window. Matches Davidson ("average of all the hydrophobicity scores in that window") and Biopython `protein_scale` with edge weight 1.0 (every weight 1 ⇒ plain mean). ✅

### Edge-case semantics
Window undefined below `w` residues ⇒ no segment (Davidson). Non-standard residues (X/B/Z/*) have no scale value ⇒ excluded from the mean (Biopython scale coverage). Both sourced and standard. ✅

### Independent cross-check (hand computation this session)
For `D×10 L×20 D×10` (40 res, w=19, T=1.6) I recomputed the full 22-point profile in Python from the scale (D=−3.5, L=3.8):
- above-threshold profile indices = **5..16**; first=5, last=16; peak = 3.8 (any all-L window).
- The **first** passing window starts at residue 5; the **last** passing window starts at residue 16 and covers residues **16..34**. The union of all passing windows' residue coverage = residues **5..34**.

### Findings / divergences (Stage A)
- **Boundary-mapping convention is not uniquely prescribed by the sources.** Two documented conventions exist: score-at-first-residue (Davidson) and score-at-midpoint (QIAGEN/most tools). Either way, the biologically meaningful *segment* is the set of residues lying within at least one above-threshold window — i.e. `[firstProfileIndex, lastProfileIndex + windowSize − 1]`. Stage A therefore fixes the **correct sourced End = `lastProfileIndex + windowSize − 1`** (last covered residue), not `+ windowSize`. The description (algorithm doc/Evidence/TestSpec) previously stated `+ windowSize`, which names a residue *outside* every passing window — corrected this session. Stage A = PASS-WITH-NOTES (parameters/scale/profile fully sourced; the single non-source-prescribed coordinate is now defined consistently with the union-of-windows semantics).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs:687-786`.
- Scale dictionary (lines 751-757): matches the sourced KD table exactly. ✅
- Profile (lines 764-786): arithmetic mean over scored residues in each window; non-standard residues excluded from sum and count (`count`-normalised); short/null/empty/`windowSize≤0` ⇒ empty. ✅
- Scan (lines 704-744): opens a run at first `P(i) ≥ T`, tracks peak, closes at first `P(i) < T`; same for trailing run.

### Defect found and fixed (off-by-one End)
The original close logic computed `end = (firstFailingIndex − 1) + windowSize = lastPassingIndex + windowSize` and reported it as the **0-based inclusive** End. For M1 this yields **35**, but the last residue covered by any passing window is **34** (last passing window starts at 16, covers 16..34). Residue 35 is a hydrophilic `D` not inside any above-threshold window — a genuine off-by-one in the reported coordinate.

Root cause: the variable `end` was being used as a *half-open* boundary for the span filter (`end − start` = inclusive residue count) but then reported *as if inclusive* without subtracting 1. Fix keeps the inclusive-span filter and reports the last **covered** residue:
```
int lastCoveredResidue = (lastPassingIndex) + windowSize - 1;
if (lastCoveredResidue - start + 1 >= MinTransmembraneHelixLength)
    yield return (start, Math.Min(lastCoveredResidue, length - 1), maxScore);
```
Verified the span filter still passes a single-window region (span = windowSize = 19 = MinLength): M3/M4/S1 → (0,18). M1 → (5,34). Legacy `ProteinMotifFinderTests` cases unaffected (their Ends were clamped to length−1, masking the off-by-one; new value equals old after clamp; Multi test asserts Start only).

### Cross-verification table (recomputed vs fixed code)
| Case | Input | Sourced expected (Start, End, Score) | Code output |
|------|-------|--------------------------------------|-------------|
| M1 | D×10 L×20 D×10 | (5, **34**, 3.8) — hand-computed profile, union of passing windows | (5, 34, 3.8) ✅ |
| M2 | D×40 | none (D=−3.5 < 1.6) | empty ✅ |
| M3 | L×19 | (0, 18, 3.8) — single window, clamped to len−1 | (0, 18, 3.8) ✅ |
| M4 | I×19 / V×19 / R×19 | (0,18,4.5) / (0,18,4.2) / none (−4.5<1.6) | matches ✅ |
| M5/M6/M7 | null / "" / L×18 | empty | empty ✅ |
| S1 | L×9 + X + L×9 | (…, 3.8) — X excluded, mean over 18 L | score 3.8 ✅ |
| S2 | D×10 A×20 D×10, T=2.0 | none (peak 1.8 < 2.0) | empty ✅ |
| S4 | windowSize=0 | empty | empty ✅ |

### Variant/delegate consistency
Single public method; MCP wrapper `AnalysisTools.PredictTransmembraneHelices` simply forwards. No `*Fast` variant. The disorder predictor in the same file uses a *different* (midpoint) convention `windowSize/2`, which is a separate unit and out of scope.

### Test quality audit (HARD gate)
- **Code-echo:** none. M1/M3/M4 assert exact sourced (Start, End, Score) derived from the hand-computed profile + sourced scale, not from code output. M1/S3 End were **wrong** (35, a code-echo of the buggy convention) → rewritten to the sourced value **34**.
- **No green-washing:** exact equality on all numeric assertions (`.Within(1e-10)` for floats); no ranges/`Greater`/`Contains` where an exact value is known. Tightened the property test `PredictTransmembraneHelices_Positions_WithinBounds` from `End ≤ length` (allowed an OOB index, masking the bug) to `End < length` + `Start ≤ End` (INV-02).
- **Coverage:** all 12 TestSpec cases (M1–M7, S1–S4, C1) present; every Stage-A branch exercised (single segment, all-hydrophilic, single window, scale reproduction, null/empty/short, non-standard residue, custom threshold, lowercase, windowSize≤0, invariants). Note: the `MinTransmembraneHelixLength == windowSize` span filter cannot reject anything at the default window (a single passing window already spans exactly 19); it is effectively a guard, documented, not separately testable without a >window MinLength.
- **Honest green:** full unfiltered suite **Failed: 0, Passed: 6579**; `dotnet build` 0 errors; changed files warning-free (the 4 NUnit2007 warnings are pre-existing in unrelated `ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects (Stage B)
- **FIXED** off-by-one End coordinate (reported a residue outside every passing window). Code + tests + doc/Evidence/TestSpec corrected to `lastProfileIndex + windowSize − 1`.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES — parameters (window 19, threshold 1.6), KD scale (20 values), and arithmetic-mean profile are fully and independently source-confirmed; the only non-source-prescribed item (segment-end coordinate) is now defined as the last residue covered by a passing window, consistent across code/doc/Evidence/TestSpec.
- **Stage B:** PASS — code realises the validated formula; the one defect (off-by-one End) was completely fixed in-session, tests locked to sourced values, full suite green.
- **End-state:** ✅ CLEAN.
