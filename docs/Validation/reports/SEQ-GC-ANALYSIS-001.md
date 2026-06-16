# Validation Report: SEQ-GC-ANALYSIS-001 — Comprehensive GC Analysis

- **Validated:** 2026-06-16   **Area:** Composition
- **Canonical method(s):** `GcSkewCalculator.AnalyzeGcContent(DnaSequence, windowSize, stepSize)` and the
  `string` overload (`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs:310-401`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

This unit aggregates five independently-defined quantities over a DNA sequence: overall GC content (%),
overall GC skew, overall AT skew, per-window GC-skew / GC-content profiles (sliding window), and the
population variance of each per-window profile. There is no single canonical "GC analysis" paper; the
unit is a composition of standard formulas, each validated separately below.

## Stage A — Description

### Sources opened this session (all fetched live, 2026-06-16)

| # | Source | URL | What it confirmed |
|---|--------|-----|-------------------|
| 1 | Wikipedia: GC-content (cites Madigan & Martinko, *Brock Biology of Microorganisms*) | https://en.wikipedia.org/wiki/GC-content | GC% = `(G+C)/(A+T+G+C) × 100%`; numerator G+C, denominator = all four bases; expressed as a percentage (×100). |
| 2 | Wikipedia: GC skew (cites Lobry 1996, Grigoriev 1998) | https://en.wikipedia.org/wiki/GC_skew | `GC skew = (G−C)/(G+C)`; range −1…+1; **+1 ⇔ C=0**, **−1 ⇔ G=0**; positive = G-rich, negative = C-rich. |
| 3 | Biopython `Bio.SeqUtils` v1.84 (`GC_skew`, `gc_fraction`) | https://biopython.org/docs/1.84/api/Bio.SeqUtils.html | `GC_skew = (G−C)/(G+C)` "for multiple windows along the sequence"; **"Returns 0 for windows without any G/C by handling zero division errors"**; `gc_fraction` returns a float in [0,1] (fraction, not %). |
| 4 | Cuemath: Population Variance (formula + worked example) | https://www.cuemath.com/data/population-variance/ | `σ² = Σ(xᵢ−μ)²/n` (divide by N, **not** N−1); worked example {12,13,12,14,19} → μ=14, Σ(dev²)=34, σ²=34/5=**6.8**. |
| 5 | PLOS Genetics — Charneski et al. 2011 (via WebSearch) + Wikipedia GC skew | https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002283 | `AT skew = (A−T)/(A+T)` (companion to GC skew; Lobry 1996 method). |

### Formula check (each vs its source)

| Quantity | Code formula | Source formula | Match |
|----------|--------------|----------------|-------|
| Overall GC content | `(G+C)/len × 100` (`CalculateGcContent`, line 386-391; only A/T/G/C present because input is validated/uppercased) | (G+C)/(A+T+G+C)×100 — Source 1 | ✅ (denominator = sequence length = total bases, valid because DnaSequence permits only ACGT) |
| Overall GC skew | `(G−C)/(G+C)`, 0 if G+C=0 (`CalculateGcSkewCore`, line 38-45) | (G−C)/(G+C); 0 on zero division — Sources 2,3 | ✅ |
| Overall AT skew | `(A−T)/(A+T)`, 0 if A+T=0 (`CalculateAtSkewCore`, line 197-206) | (A−T)/(A+T) — Source 5 | ✅ |
| Windowed profiles | sliding window, start `i`, end `i+w−1`, position `i+w/2`, step `stepSize` (lines 85-101, 364-380) | "multiple windows along the sequence" — Source 3 | ✅ |
| Variance | `Σ(x−μ)²/N` (`CalculateVariance`, line 396-401) | population variance Σ(x−μ)²/N — Source 4 | ✅ |

### Definitions & conventions

- **Percentage vs fraction:** code reports GC content ×100 (percentage, Source 1 / Brock convention), not
  Biopython's [0,1] fraction. Documented assumption; both are standard; the unit pins the exact percentage so
  the convention is locked. PASS-with-note level, not a defect.
- **Ambiguity handling:** only G/C count for skew, only A/T for AT skew, only A/T/G/C for content. DnaSequence
  validates input to ACGT, so non-ACGT never reaches the core. Consistent with Source 3 (Biopython ignores
  non-G/C).
- **Coordinate base:** window `WindowStart` 0-based, `WindowEnd` 0-based inclusive (`i+w−1`), `Position` =
  midpoint `i + w/2`. Internally consistent and standard.
- **Population vs sample variance:** population (÷N). Documented assumption (the windows *are* the entire
  population of windows for the sequence). Sample variance (÷N−1) would change M4 from 1.0 to 2.0 — the test
  discriminates this.

### Edge-case semantics (all sourced)

- G+C=0 → GC skew 0 (Source 3, zero-division). ✅
- A+T=0 → AT skew 0 (same convention). ✅
- No G/C → GC% 0 (numerator 0, Source 1). ✅
- Pure-G window → skew +1 (Source 2, "+1 ⇔ C=0"); pure-C window → skew −1 (Source 2, "−1 ⇔ G=0"). ✅
- Sequence shorter than window → no full window → empty profiles → window-variances 0; overall scalars still
  computed over the whole sequence (implementation windowing contract; the loop `i+w ≤ n` simply never
  enters). Reasonable and documented.

### Independent cross-check (hand-computed from the sourced formulas, not from code output)

| Dataset | Quantity | Sourced computation | Value |
|---------|----------|---------------------|-------|
| `GGGCCAT` (G3 C2 A1 T1, n7) | GC% | (3+2)/7×100 | 71.42857142857143 |
| `GGGCCAT` | GC skew | (3−2)/(3+2) | 0.2 |
| `GGGCCAT` | AT skew | (1−1)/(1+1) | 0.0 |
| `GGCC` w2 s2 → GG(+1),CC(−1) | GcSkewVariance | μ=0, ((1)²+(−1)²)/2 | 1.0 |
| `GGCC` → 100,100 | GcContentVariance | var{100,100} | 0.0 |
| `AAAGGGCCCTTT` w3 s3 → 0,100,100,0 | GcContentVariance | μ=50, (2500·4)/4 | 2500 |
| `ACGTACGTAC` (n10) w4 s2 | window count | ⌊(10−4)/2⌋+1 | 4 |
| `GGGG` | GC skew / GC% | (4−0)/4 / 4/4×100 | +1 / 100 |
| `CCCC` | GC skew / GC% | (0−4)/4 / 4/4×100 | −1 / 100 |
| `ATATAT` (A3 T3) | GC skew / GC% / AT skew | 0 / 0 / (3−3)/6 | 0 / 0 / 0 |
| `AAAT` (A3 T1) | AT skew | (3−1)/(3+1) | 0.5 |
| {12,13,12,14,19} | population variance (anchor) | 34/5 | 6.8 |

Every non-trivial expected value in the test fixture traces to one of these sourced computations.

### Invariants

INV-1…INV-3 (skew ∈ [−1,1], GC% ∈ [0,100]) are true by the formula ranges (Sources 1,2). INV-4 (variance ≥ 0)
true (sum of squares ÷ N). INV-5 (window count = ⌊(n−w)/step⌋+1, 0 if n<w) matches the loop `i+w ≤ n` — verified
by hand on M7 (4) and S1 (0). All genuine.

**Stage A findings:** none. All formulas, conventions and edge cases match the retrieved authoritative
sources. **Verdict: PASS.**

## Stage B — Implementation

### Code path reviewed

`AnalyzeGcContent(DnaSequence)` (line 310) → `AnalyzeGcContentCore` (line 336), delegating to
`CalculateGcContent` (386), `CalculateGcSkewCore` (38), `CalculateAtSkewCore` (197),
`CalculateWindowedGcSkewCore` (85), `CalculateWindowedGcContentCore` (364), `CalculateVariance` (396).
The `string` overload (line 325) uppercases and delegates to the same core; null/empty → zero result.

- **Formula realised correctly:** confirmed line-by-line against the Stage-A table. `CalculateGcContent`
  does not itself uppercase, but every caller passes already-uppercased text (`DnaSequence.Sequence` is
  normalised in the ctor, `DnaSequence.cs:30`; the string overload calls `ToUpperInvariant()`), so case is
  handled. Verified by the C1b lowercase test.
- **Edge cases in code:** zero-division guards present (`total > 0 ? … : 0`) in both skew cores; variance
  guarded for empty input; windowed cores never enter the loop when `n < w` (so empty list → variance 0 via
  the `Count > 0 ? … : 0` ternary). All Stage-A edge cases handled.
- **Variant consistency:** string and DnaSequence overloads share `AnalyzeGcContentCore` — equivalence
  proven by C1.
- **Numerical robustness:** double arithmetic; no overflow on stated ranges; div-by-zero guarded.

### Cross-verification table recomputed vs code

The full unfiltered suite was run (`dotnet test`, 6609 passed / 0 failed). Every value in the Stage-A
cross-check table is asserted by a test and passes against the actual code.

### Test quality audit (HARD gate)

Fixture: `tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculator_AnalyzeGcContent_Tests.cs`.

- **Sourced expectations, not code echoes:** every MUST/SHOULD assertion uses an exact value derived from the
  external formula (71.42857142857143, 0.2, 1.0, 2500, ±1, 0.5, …), each within 1e-10. The comments even call
  out the discriminating alternatives a wrong implementation would produce (sample variance 2.0 vs population
  1.0; fraction vs percentage; `(G−C)/n` vs `(G−C)/(G+C)`). No `Greater`/`AtLeast`/range stand-ins where an
  exact value is known. No green-washing.
- **Branch coverage strengthened this session (test-quality, no algorithm defect):**
  the documented **lower bound** "pure-C window → skew −1" (Evidence corner case; Source 2 "−1 ⇔ G=0") was
  only exercised *indirectly* inside the M4 variance and never asserted directly; and **AT skew was only ever
  tested at 0** (M3, M9), so a sign-flip `(T−A)/(A+T)` or a wrong-denominator `(A−T)/n` would have passed.
  Added two exact-sourced tests:
  - `AnalyzeGcContent_PureCytosine_SkewIsMinusOneContentIsHundred` — `CCCC` → GC skew −1, GC% 100.
  - `AnalyzeGcContent_AdenineRich_AtSkewIsPositiveHalf` — `AAAT` (A3 T1) → AT skew (3−1)/(3+1)=+0.5.
- **Cover all the logic:** both overloads exercised (DnaSequence canonical + string delegate + lowercase);
  population variance pinned (M4/M5/M6) with the sample-variance trap; window geometry (count/start/end/
  position) pinned (M7); short-sequence, null-DnaSequence (throws), null/empty-string (zero result) all
  covered (S1/S2/S3); both skew bounds now pinned directly (M8 +1, M8b −1); AT-skew sign+formula pinned
  (M3b). Every public overload and every Stage-A branch is exercised.
- **Honest green:** full unfiltered `dotnet test` = **Passed 6609, Failed 0** (was 6607 before the 2 added
  tests). `dotnet build` = 0 errors; the 4 build warnings are pre-existing `NUnit2007` notes in an unrelated
  file (`ApproximateMatcher_EditDistance_Tests.cs`) — the changed fixture builds warning-free.

**Result of test-quality gate: PASS** (two coverage gaps closed with exact sourced values; suite honest-green).

### Findings / defects

No algorithm defect. One test-coverage finding (documented lower-bound + AT-skew-direction gap), fixed in
this session — logged as a FIXED-NOW test-quality item in the findings register.

## Verdict & follow-ups

- **Stage A: PASS.** **Stage B: PASS.** **End-state: ✅ CLEAN.**
- Algorithm fully functional; all formulas independently confirmed against retrieved authoritative sources;
  tests now pin both skew bounds and the AT-skew formula/sign with exact sourced values.
- Note (not a defect): GC content uses the percentage (×100) convention rather than Biopython's [0,1]
  fraction — a documented, locked units choice.
