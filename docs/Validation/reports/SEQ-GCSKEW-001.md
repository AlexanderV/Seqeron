# Validation Report: SEQ-GCSKEW-001 — GC Skew Analysis

- **Validated:** 2026-06-12   **Area:** Sequence Composition
- **Canonical method(s):** `GcSkewCalculator.CalculateGcSkew(...)`, `CalculateWindowedGcSkew(...)`, `CalculateCumulativeGcSkew(...)`, `CalculateAtSkew(...)`, `PredictReplicationOrigin(...)`, `AnalyzeGcContent(...)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculatorTests.cs`, `tests/Seqeron/Seqeron.Genomics.Tests/Properties/GcSkewProperties.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "GC skew"** (fetched 2026-06-12) confirms verbatim:
  - GC skew = (G − C)/(G + C)
  - AT skew = (A − T)/(A + T)
  - Cumulative skew "takes advantage of the sum of the adjacent windows from an arbitrary start … the maximum value of the cumulative skew corresponds to the terminal, and the minimum value corresponds to the origin of replication."
  - Range −1 to +1; "Positive GC skew represents richness of G over C and the negative GC skew represents richness of C over G."
  - Leading strand: positive GC skew / negative AT skew; lagging strand: the reverse.
- **Lobry (1996), Mol. Biol. Evol. 13:660–665** — first report of strand compositional asymmetry in *E. coli*, *B. subtilis*, *H. influenzae*. Note: Lobry's original convention was the **sign-flipped** form (C−G)/(C+G). The modern standard (Wikipedia, Grigoriev) is (G−C)/(G+C), which is what the spec and code adopt. This is an internally consistent choice, not an error — see Findings.
- **Grigoriev (1998), Nucleic Acids Res. 26:2286–2290** — cumulative GC skew method; with the (G−C)/(G+C) definition, global minimum = origin (oriC), global maximum = terminus (ter).
- **Biopython Bio.SeqUtils.GC_skew** — formula (G−C)/(G+C), returns 0.0 on ZeroDivisionError (G+C=0), case-insensitive.

### Formula check
- GC skew = (G − C)/(G + C) — matches Wikipedia exactly. **Sign is correct (not reversed).**
- AT skew = (A − T)/(A + T) — matches Wikipedia exactly.
- Cumulative GC skew = running sum of per-window skews — matches Grigoriev/Wikipedia.
- Origin = minimum of cumulative skew; Terminus = maximum — matches Wikipedia/Grigoriev under the (G−C)/(G+C) convention.
- Div-by-zero (G+C=0) → 0 — matches Biopython convention.

### Edge-case semantics check
- Empty / all-A,T (G+C=0): skew = 0 (Biopython ZeroDivisionError → 0.0). Defined and sourced.
- All-G → +1; All-C → −1; equal G,C → 0. Sourced to formula.
- Windowed on empty / window > length: empty result (no windows fit). Sourced to Biopython `GC_skew("") → []`.

### Independent cross-check (hand-computed numbers)
- GGGGC: (4−1)/(4+1) = **0.6** ✓
- GCCC: (1−3)/(1+3) = **−0.5** ✓ (Biopython calc_gc_skew)
- GGGGCCCC, w=4, step=4: GGGG→(4−0)/4 = **+1.0**, CCCC→(0−4)/4 = **−1.0** ✓ (Biopython GC_skew → [1.0, −1.0])
- ATGCATGC, w=4: each window G=1,C=1 → **0.0, 0.0** ✓
- AAAAAAAA, w=4: G+C=0 → **0.0, 0.0** ✓
- AT skew AAAAT: (4−1)/(4+1) = **0.6** ✓
- **Origin detection worked example** (`PredictReplicationOrigin_FindsMinimum`): G×50 + C×100 + G×50, window=10 (cumulative step = window). Windows 0–4 are G (skew +1, cumulative climbs +1…+5); windows 5–14 are C (skew −1, cumulative drops +4…−5); last C-window at i=140, center = 140 + 10/2 = **145**, cumulative = +5 − 10 = **−5** = global minimum → PredictedOrigin = 145, OriginSkew = −5. Matches spec/test exactly.
- **Terminus worked example** (mirror, C×50+G×100+C×50): maximum +5 at center **145** → PredictedTerminus = 145. Matches.

### Findings / divergences
- **NOTE (sign convention, documented & consistent):** The spec/code use the modern (G−C)/(G+C) convention, not Lobry's original (C−G)/(C+G). The spec explicitly states this (Evidence row for Lobry 1996 and Grigoriev 1998). The convention is self-consistent: min=origin / max=terminus is the correct pairing for (G−C)/(G+C). PASS-WITH-NOTES on this account only.
- **NOTE (doc drift, non-blocking):** `ALGORITHMS_CHECKLIST_V2.md` line 1543 lists a method `FindOriginOfReplication(sequence)`; the actual method is `PredictReplicationOrigin`. Documentation-only naming drift in the checklist; the TestSpec uses the correct name and the code/tests are consistent. No functional impact.

## Stage B — Implementation

### Code path reviewed
`GcSkewCalculator.cs`:
- `CalculateGcSkewCore` (line 38–45): `(gCount - cCount) / total`, guarded `total > 0 ? … : 0`. Correct formula, correct sign, div-by-zero guarded.
- `CalculateAtSkewCore` (184–191): `(aCount - tCount) / total`, same guard. Correct.
- `CalculateWindowedGcSkewCore` (85–101): windows `for (i=0; i+windowSize<=len; i+=stepSize)`; position = `i + windowSize/2` (window center, per Grigoriev). Correct.
- `CalculateCumulativeGcSkewCore` (139–157): step = windowSize (non-overlapping tiling), `cumulative += skew` running sum. Correct.
- `PredictReplicationOrigin` (203–229): `MinBy(CumulativeGcSkew)` → origin, `MaxBy` → terminus. Empty → default zeros. Correct mapping.

### Formula realised correctly?
Yes. Confirmed (G−C)/(G+C) with the correct (non-reversed) sign; cumulative sum correct; window bounds correct; div-by-zero guarded for both GC and AT skew.

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

### Variant/delegate consistency
String and DnaSequence overloads both route to the same `*Core` methods; string overload upper-cases input (case-insensitive, per Biopython). `AnalyzeGcContent` reuses the same core skew functions. Consistent.

### Numerical robustness
Integer counts, double division. No overflow on stated ranges; all div-by-zero paths guarded. Range invariant verified by FsCheck (`GcSkew_InRange`, `AtSkew_InRange`, `Complement_NegatesGcSkew`, 100 runs each).

### Test quality audit
44 NUnit methods + 3 FsCheck properties assert exact sourced values (0.6, −0.5, ±1.0, position 145, cumulative ±5), not tautologies. Edge cases (empty, null, G+C=0, window>length, zero/negative window/step) covered. Tests are deterministic.

### Findings / defects
None. Implementation faithfully realises the validated description.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — formula and conventions match authoritative sources; modern (G−C)/(G+C) convention is explicitly adopted (differs from Lobry's original sign, by design) and is internally consistent with min=origin/max=terminus.
- **Stage B: PASS** — code computes (G−C)/(G+C) with correct sign, correct cumulative sum, correct window handling, guarded div-by-zero; all worked examples recomputed and match.
- **State: CLEAN** — no defect found. No code changes required. Build green; GcSkew filter 44/44; full suite 4461/4461.
- Optional non-blocking follow-up (outside this unit's code): fix the stale method name `FindOriginOfReplication` → `PredictReplicationOrigin` in `ALGORITHMS_CHECKLIST_V2.md` line 1543.
