# Validation Report: PROTMOTIF-LC-001 — Low-Complexity Region Detection (SEG)

- **Validated:** 2026-06-16   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.FindLowComplexityRegions(string, int windowSize, double triggerComplexity, double extensionComplexity)`; internal `CalculateSegComplexity(ReadOnlySpan<char>)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES (two test-coverage defects found and fixed in-session)

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

- **NCBI `ncbi-seg` man page** (https://manpages.ubuntu.com/manpages/focal/man1/ncbi-seg.1.html, WebFetched 2026-06-16) — verbatim defaults: trigger window length W "[Default 12]"; trigger complexity K1 "in units of bits … maximum value is 4.322 (log[base 2]20) for amino acid sequences [Default 2.2]"; extension complexity K2 "[Default 2.5]", "Only values greater than K1 are effective in extending triggered windows." Confirms W=12, K1=2.2, K2=2.5, complexity in **bits/residue**, range 0…log₂20.
- **GCG/Rothlab SEG docs** (https://rothlab.ucdavis.edu/genhelp/seg.html, WebFetched) — "Seg identifies segments having a complexity equal to or less than the cutoff **in bits/residue**" (K1, stage 1); extends segments "≤ the cutoff in bits/residue" (K2, stage 2). "If 20 different characters were distributed randomly … each character would add 4.322 bits (log(base 2) 20)." Confirms the cutoff metric is a per-residue bits measure with max log₂20.
- **Weizmann GCG SEG doc** (https://bip.weizmann.ac.il/education/materials/gcg/seg.html, WebFetched) — LOWcut=K1, HIGHcut=K2, WINdow default 12, cutoffs in bits/residue.
- **Pei & Grishin 2005, Bioinformatics 21(2):160** (WebSearch result text) — "two-pass algorithm" with "default parameters W = 12, K2(1) = 2.2 bit and K2(2) = 2.5 bits"; "SEG uses a complexity measure based on Shannon entropy"; trigger windows of length W with complexity ≤ K(1), extension ≤ K(2).
- **SeqComplex `SeqComplex.pm`** (https://raw.githubusercontent.com/caballero/SeqComplex/master/SeqComplex.pm, WebFetched) — `sub ce` (Shannon entropy): `ce -= r * log_k(2, r)` for each symbol with count>0, r=count/tot → K = −Σ pᵢ·log₂ pᵢ; guarded `if tot>1`. `sub cwf` (Wootton-Federhen compositional form): `(up − dw)/tot` with `up = Σ log_N(W)`, `dw = Σ_b log_N(n_b)`, i.e. the (1/L)·log_N(L!/Πnᵢ!)-style measure normalized by alphabet base N.

### Formula check

The implementation computes per-window K = −Σᵢ pᵢ·log₂(pᵢ), pᵢ = nᵢ/L, in bits/residue.

**Important nuance (the eq.(1) vs eq.(3) distinction).** Wootton & Federhen (1993) define two complexity quantities: the Shannon-entropy form (their eq.(1), the "complexity state vector" entropy, the bits/residue measure compared to K1/K2) and the exact multinomial compositional complexity P₀ (their eq.(3), the log-factorial form, NCBI `s_LnPerm`/`lnfact[]`). The **trigger/extension cutoffs K1/K2** are applied to the **bits/residue Shannon-entropy measure** — this is what every retrieved operational source states (man page "bits, max log₂20"; GCG "bits/residue"; Pei & Grishin "complexity measure based on Shannon entropy"; SeqComplex `ce`). The multinomial P₀ (eq.(3)) is used by NCBI SEG only in the **optional pass-2 local optimization** that selects the P₀-minimizing sub-segment within a raw segment. The repo's Evidence/TestSpec loosely label the entropy measure "the bits/residue form of eq.(3)"; strictly it is closer to eq.(1), but the *operational* claim (K1/K2 are compared against the window's Shannon entropy in bits/residue, range 0…log₂20) is **correct and matches all sources**. → recorded as a NOTE, not a defect.

### Edge-case semantics check

- Sequence shorter than window → no complete trigger window → empty result (man page). ✓
- Homopolymer window → single symbol p=1, −1·log₂1 = 0 (lowest complexity). ✓
- K2 must exceed K1 to extend ("Only values greater than K1 are effective"). ✓
- Bound 0 ≤ K ≤ log₂20 for amino-acid windows. ✓

### Independent cross-check (numbers, computed this session via Python `math.log2`)

| Window composition (L=12) | K = −Σ pᵢ·log₂ pᵢ |
|---------------------------|-------------------|
| 12×A | 0.0 |
| 11×A,1×B | 0.41381685030363 |
| 6×A,6×B | 1.0 |
| 10×A,2×B | 0.65002242164835 |
| 12 distinct | log₂12 = 3.58496250072116 |
| max (20 aa) | log₂20 = 4.32192809488736 |

All match the implementation/test expected values exactly.

### Stage A findings / divergences

- **NOTE:** the entropy measure is the Shannon-entropy (eq.(1)) bits/residue form, correctly identified by units and range; Evidence's "form of eq.(3)" label is imprecise but the operational definition is right.
- **NOTE (documented simplification):** the optional P₀-minimized pass-2 local optimization (eq.(3) multinomial) is intentionally omitted; region boundaries are the window-run span, not the P₀-optimal sub-segment. This is disclosed in the algorithm doc §5.3 and limitations §6.2 — a faithful simplification, not an overclaim.

**Stage A verdict: PASS-WITH-NOTES.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs:1080–1159`.

- `CalculateSegComplexity` (1140–1159): `stackalloc int[char.MaxValue+1]` counts; sums −p·log₂(p) over each distinct residue once (slot cleared after consumption). Realises K = −Σ pᵢ·log₂ pᵢ exactly. For 12×A returns 0. ✓
- `FindLowComplexityRegions` (1080–1132): validates `windowSize>0` (else `ArgumentOutOfRangeException`); null/empty/`<W` → `yield break`; upper-cases input; precomputes per-window K; scans maximal runs of windows with K ≤ K2; emits a run as `(runStart, runEnd+W−1, minK)` only if at least one window has K ≤ K1 (triggered). Realises the two-pass trigger/extension rule and 0-based inclusive span. ✓

### Cross-verification table recomputed vs code (independent Python)

| Case | Independent expectation | Code behaviour |
|------|------------------------|----------------|
| M1 12×A complexity | 0.0 | 0.0 ✓ |
| M2 11A1B | 0.413817 | match ✓ |
| M3 6A6B | 1.0 | match ✓ |
| M4 10A2B | 0.650022 | match ✓ |
| M5 12 distinct | log₂12=3.584963 | match ✓ |
| M6 poly-Q (flank+Q×20+flank) | 1 region, span [6,37], minK=0 | match ✓ |
| M7 20 distinct | min window K=3.585>K2 → empty | empty ✓ |
| M9 two tracts (G×14+spacer+S×14) | regions [0,20] & [20,39], minK=0 | match ✓ |
| **M6b** (new) `AAABBBCCDDEE` | K=2.2925 ∈ (K1,K2] → not triggered → empty | empty ✓ |

The M9 boundary chain was re-derived window-by-window: window 9 K=2.2925 ≤ K2, window 10 K=2.6175 > K2 → poly-G span [0,20]; window 20 K=2.2925 ≤ K2 → poly-S span [20,39]. Matches code and the test's locked boundaries.

### Variant/delegate consistency

Single public entry point + one internal helper; defaults (W=12, K1=2.2, K2=2.5) verified equal to the explicit SEG defaults (M8). No `*Fast`/instance variants.

### Numerical robustness

`Math.Log2` on p∈(0,1]; p=0 terms skipped; homopolymer → 0 exactly. No div-by-zero (L=windowSize>0). `stackalloc int[char.MaxValue+1]` (~256 KB) — fine for the protein alphabet; any char is a valid index after upper-casing.

### Test quality audit (HARD gate)

Existing tests M1–M5 assert exact externally-sourced entropy literals (would fail a natural-log or normalized impl) — strong. M7/M8/M9 assert exact emptiness/defaults/boundaries — strong. Edge cases (null/empty, `<W`, non-positive window throw, bounds, case-insensitivity, determinism) covered.

**Two defects found in the test set** (no code defect):

1. **M6 boundary gap** — Evidence MUST-test #5 requires the poly-Q region with *correct inclusive boundaries*, but M6 asserted only `Count==1` and `Complexity==0.0`, not Start/End. Boundaries [6,37] are externally computable. **Fixed:** added Start=6, End=37 assertions (independently derived).
2. **K1-trigger suppression gap (the load-bearing one)** — no test exercised the two-pass *trigger* rule: a run that stays ≤K2 but never reaches ≤K1 must NOT be emitted. A deliberately-wrong implementation that emitted any run ≤K2 (ignoring K1) would have passed the entire suite — exactly the "passes against a wrong impl" defect the gate forbids. **Fixed:** added `FindLowComplexityRegions_WindowAboveK1WithinK2_NotTriggered_ReturnsEmpty` using `"AAABBBCCDDEE"` (counts 3,3,2,2,2 → K=2.292481, independently computed, ∈ (K1=2.2, K2=2.5]); with SEG defaults this window is in the extension band but never triggers → empty result. This locks the K1 trigger semantics.

No assertion was weakened, no tolerance widened, no test skipped. Expected values trace to the −Σp·log₂p formula computed this session, not to code output.

### Stage B findings / defects

- Two **test-coverage** defects (above), both fixed in-session with sourced expectations. No implementation defect — the code realises the validated formula and two-pass rule correctly.

**Stage B verdict: PASS-WITH-NOTES.**

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** (eq.(1) vs eq.(3) labeling nuance; documented P₀ pass-2 simplification — both honest, source-grounded).
- **Stage B: PASS-WITH-NOTES** (two test-coverage defects fixed; no code defect).
- **Test-quality gate:** PASS after fix — exact sourced values, all Stage-A branches covered (including the previously-untested K1 trigger-suppression path), full unfiltered suite green.
- **End-state: ✅ CLEAN.** `dotnet build` 0 errors; full unfiltered `dotnet test` = **6583 passed, 0 failed**.

### Build/test evidence

```
dotnet build … : 0 Error(s) (4 pre-existing warnings, unrelated to this unit)
dotnet test  … : Passed! - Failed: 0, Passed: 6583, Skipped: 0, Total: 6583
```
