# Validation Report: SEQ-GCSKEW-001 — GC Skew Analysis

- **Validated:** 2026-06-24   **Area:** Sequence Composition
- **Canonical method(s):** `GcSkewCalculator.CalculateGcSkew(...)`, `CalculateWindowedGcSkew(...)`, `CalculateCumulativeGcSkew(...)`, `CalculateAtSkew(...)`, `PredictReplicationOrigin(...)`, `AnalyzeGcContent(...)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculatorTests.cs`, `GcSkewCalculator_AnalyzeGcContent_Tests.cs`, `GcSkewCalculator_CalculateAtSkew_Tests.cs`, `GcSkewCalculator_PredictReplicationOrigin_Tests.cs`, `Properties/GcSkewProperties.cs`, `ConventionCompatibility_OptIn_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Re-validation context (commit 6e900e92)

This unit was reset because commit `6e900e92` ("opt-in Biopython/VCF compatibility modes")
touched `GcSkewCalculator.cs`. The diff was inspected in full:

- The change added a single new optional parameter `bool fraction = false` to the two
  `AnalyzeGcContent(...)` overloads and threaded it through `AnalyzeGcContentCore`,
  `CalculateWindowedGcContentCore`, and the private `CalculateGcContent` helper.
- `fraction == false` (the default) reproduces the prior behaviour exactly: GC **content**
  reported as a percentage GC% = (G+C)/(A+T+G+C)·100. `fraction == true` divides by 100,
  i.e. reports GC content in [0,1] per Biopython `gc_fraction` (which "produces a float value
  between 0 and 1").
- **The `fraction` flag affects GC *content* only. It does NOT touch any skew path.**
  `CalculateGcSkewCore`, `CalculateAtSkewCore`, windowed/cumulative skew, and
  `PredictReplicationOrigin` are byte-for-byte unchanged by the commit.
- No new "ambiguity" overload was added to *this* class by the commit; the Biopython
  ambiguity modes (`GcAmbiguityMode {Remove,Ignore,Weighted}`) were added to
  `SequenceExtensions.CalculateGcFraction` (a different unit, SEQ-GC-001).

Net effect on SEQ-GCSKEW-001: a purely additive, default-preserving opt-in on the
composite `AnalyzeGcContent` method. Defaults unchanged; the prior CLEAN verdict still holds.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "GC skew"** (fetched 2026-06-24) confirms verbatim:
  - "GC skew = (G - C)/(G + C)" and "AT skew = (A − T)/(A + T)".
  - Range "from −1, which corresponds to G = 0 or A = 0, to +1, which corresponds to T = 0
    or C = 0."
  - Cumulative: "the maximum value of the cumulative skew corresponds to the terminal, and
    the minimum value corresponds to the origin of replication."
  - "In contrast to Lobry's earlier paper, recent implementations of GC skew flips the
    original definition, redefining it to be: GC skew = (G − C)/(G + C)." — confirms the
    modern convention the code adopts, and that the min=origin/max=terminus pairing is the
    correct one for it.
- **Biopython `Bio.SeqUtils.GC_skew`** (fetched 2026-06-24): formula "(G-C)/(G+C)";
  "Returns 0 for windows without any G/C by handling zero division errors." `gc_fraction`
  "produces a float value between 0 and 1" — confirms both the zero-division → 0 convention
  and the [0,1] semantics for the opt-in `fraction:true` mode.
- **Lobry (1996) MBE 13:660–665** / **Grigoriev (1998) NAR 26:2286–2290** — original strand
  asymmetry observation; cumulative GC-skew method for oriC/ter detection. Lobry's original
  sign was (C−G)/(C+G); modern standard (used here) is (G−C)/(G+C).

### Formula check
- GC skew = (G − C)/(G + C) — matches Wikipedia exactly; sign correct (not reversed).
- AT skew = (A − T)/(A + T) — matches Wikipedia exactly.
- Cumulative GC skew = running sum of per-window skews — matches Grigoriev/Wikipedia.
- Origin = minimum, Terminus = maximum of cumulative skew — correct pairing for (G−C)/(G+C).
- Div-by-zero (G+C=0) → 0 — matches Biopython.
- GC content opt-in: percentage (default) vs [0,1] fraction (`fraction:true`) — matches
  Biopython `gc_fraction` range.

### Edge-case semantics check
- Empty / all-A,T (G+C=0): skew = 0 (Biopython ZeroDivisionError → 0.0). Defined & sourced.
- All-G → +1; All-C → −1; equal G,C → 0. Sourced to formula.
- Windowed on empty / window > length: empty result (no windows fit). Sourced.
- Null DnaSequence → ArgumentNullException; window/step ≤ 0 → ArgumentOutOfRangeException.

### Independent cross-check (hand-computed numbers)
- GGGGC: (4−1)/(4+1) = **0.6** ✓
- GCCC: (1−3)/(1+3) = **−0.5** ✓ (Biopython calc_gc_skew)
- GGGGCCCC, w=4, step=4: GGGG→**+1.0**, CCCC→**−1.0** ✓ (Biopython GC_skew → [1.0,−1.0])
- ATGCATGC, w=4: each window G=1,C=1 → **0.0, 0.0** ✓
- AAAAAAAA, w=4: G+C=0 → **0.0, 0.0** ✓
- AT skew AAAAT: (4−1)/(4+1) = **0.6** ✓
- Origin worked example (G×50 + C×100 + G×50, w=10): cumulative min **−5** at center **145**
  → PredictedOrigin = 145. Terminus mirror → max **+5** at **145**. ✓
- **Opt-in fraction cross-check** (GCATGCAT, w=4): default OverallGcContent = **50.0** (%),
  `fraction:true` → **0.5** ([0,1]); skew values identical in both. ✓ (matches Biopython
  `gc_fraction` and the `ConventionCompatibility_OptIn_Tests`).

### Findings / divergences
- **NOTE (sign convention, documented & consistent):** modern (G−C)/(G+C), not Lobry's
  original (C−G)/(C+G). Spec states this explicitly; self-consistent with min=origin/max=ter.
  → PASS-WITH-NOTES on this account only.
- TestSpec §2 lists the new opt-in only implicitly; the `fraction` parameter is documented in
  code XML-docs and covered by `ConventionCompatibility_OptIn_Tests`. Non-blocking.

## Stage B — Implementation

### Code path reviewed
`GcSkewCalculator.cs`:
- `CalculateGcSkewCore` (38–45): `(gCount − cCount)/total`, guarded `total > 0 ? … : 0`. Correct.
- `CalculateAtSkewCore` (197–206): `(aCount − tCount)/total`, same guard. Correct.
- `CalculateWindowedGcSkewCore` (85–101): `for (i=0; i+windowSize<=len; i+=stepSize)`,
  position = `i + windowSize/2` (window center). Correct.
- `CalculateCumulativeGcSkewCore` (139–157): step = windowSize (non-overlapping tiling),
  `cumulative += skew`. Correct.
- `PredictReplicationOriginCore` (255–286): per-base cumulative skew (G:+1, C:−1, A/T:0),
  tracks first global min/max prefix index; min→origin, max→terminus; empty → zero
  prediction `IsSignificant=false`. Matches Grigoriev/Rosalind BA1F convention.
- `AnalyzeGcContentCore` (338–366) + `CalculateGcContent` (391–417): the `fraction` flag
  scales GC content by 1.0 (fraction) or 100.0 (percent, default). Skew values unaffected.
  Verified the flag is NOT plumbed into any skew computation.

### Formula realised correctly?
Yes. (G−C)/(G+C) with correct sign; AT (A−T)/(A+T); cumulative running sum; window center
position; div-by-zero guarded for both GC and AT. Opt-in `fraction` correctly limited to GC
content with default preserved.

### Cross-verification table recomputed vs code (all pass)
| Input | Expected | Test | Result |
|-------|----------|------|--------|
| GGGGC | 0.6 | CalculateGcSkew_MoreG | ✓ |
| CCCCC | −1.0 | CalculateGcSkew_MoreC | ✓ |
| GGGGG | +1.0 | CalculateGcSkew_AllG | ✓ |
| GCGC | 0.0 | CalculateGcSkew_EqualGC | ✓ |
| AAATTT / "" | 0 | NoGC / EmptySequence | ✓ |
| GCCC | −0.5 | BiopythonCrossVerification_GCCC | ✓ |
| GGGGCCCC w=4 | [1.0,−1.0] | CorrectSkewValues | ✓ |
| ATGCATGC w=4 | [0.0,0.0] | CrossVerification_ATGCATGC | ✓ |
| AAAAAAAA w=4 | [0.0,0.0] | CrossVerification_AllA | ✓ |
| AAAAT (AT skew) | 0.6 | CalculateAtSkew_MoreA | ✓ |
| Origin G50C100G50 w=10 | pos=145, skew=−5 | FindsMinimum | ✓ |
| Terminus C50G100C50 w=10 | pos=145, skew=+5 | FindsMaximum | ✓ |
| GCATGCAT content % / frac | 50.0 / 0.5 | ConventionCompatibility (AnalyzeGcContent fraction) | ✓ |

### Variant/delegate consistency
String and DnaSequence overloads route to the same `*Core` methods; string overloads
upper-case input (case-insensitive, per Biopython). `AnalyzeGcContent` reuses the same skew
cores. New `fraction` parameter defaults false → identical to pre-commit behaviour. Consistent.

### Numerical robustness
Integer counts, double division; all div-by-zero paths guarded. Range invariant verified by
FsCheck (`GcSkewProperties`, 100 runs each).

### Test quality audit
NUnit methods + 3 FsCheck properties assert exact sourced values (0.6, −0.5, ±1.0, pos 145,
±5, 50.0/0.5), not tautologies; edge cases (empty/null/G+C=0/window>length/zero-step) covered.
The opt-in `fraction` surface is locked by `ConventionCompatibility_OptIn_Tests`
(default % vs [0,1], windowed and overall). GcSkew + ConventionCompatibility filter: 173/173.

### Findings / defects
None. Implementation faithfully realises the validated description; the commit 6e900e92 opt-in
is additive and default-preserving.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — formula/conventions match authoritative sources; modern
  (G−C)/(G+C) convention adopted by design (differs from Lobry's original sign), internally
  consistent with min=origin/max=terminus; new `fraction:true` mode matches Biopython [0,1].
- **Stage B: PASS** — code computes (G−C)/(G+C) and (A−T)/(A+T) with correct signs, correct
  cumulative sum, correct window centers, guarded div-by-zero; `fraction` flag confined to GC
  content; all worked examples recomputed and match.
- **State: CLEAN** — no defect found, no code changes required. Defaults unchanged after
  commit 6e900e92. GcSkew + ConventionCompatibility filter: 173/173 passed.
