# Validation Report: PAT-APPROX-003 — Best Match and Frequency Analysis

- **Validated:** 2026-06-15   **Area:** Matching (Approximate Pattern Matching / Frequent Words with Mismatches)
- **Canonical method(s):** `ApproximateMatcher.FindFrequentKmersWithMismatches(sequence, k, d)`,
  `ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, maxMismatches)`,
  `ApproximateMatcher.FindBestMatch(sequence, pattern)` (plus the supporting `FindWithMismatches`
  used by M4 and the `GenerateNeighbors` recursion).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (one test-quality defect found and fixed in-session)

## Stage A — Description

### Sources opened this session
- **ROSALIND BA1I — Frequent Words with Mismatches** (https://rosalind.info/problems/ba1i/),
  retrieved via WebFetch 2026-06-15. Confirms verbatim:
  - Count_d(Text, Pattern) = "the total number of occurrences of Pattern in Text with at most d mismatches."
  - Worked example: Count_1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = **4**, windows AACAA, ATAAA, AAACA, AAAGA.
  - Sample: Text=ACGTTGCATGTCGCATGATGCATGAGAGCT, k=4, d=1 → output **GATG ATGC ATGT** (all ties returned).
  - The most-frequent k-mer need not occur exactly in Text (counting is over the Hamming ball).
- **ROSALIND BA1H — Approximate Occurrences of a Pattern** (https://rosalind.info/problems/ba1h/),
  retrieved via WebFetch 2026-06-15. Confirms verbatim:
  - Approximate occurrence ⇔ HammingDistance(Pattern, Pattern') ≤ d.
  - Sample: Pattern=ATTCTGGA, Text=CGCCCGAATCC…ACGCTCC (99 nt), d=3 → positions **6 7 26 27 78** (0-based);
    therefore Count_3 = 5.
- **ROSALIND BA1N — d-Neighborhood** (cited in Evidence): neighborhood = all k-mers with Hamming distance ≤ d,
  includes the pattern itself. This matches `GenerateNeighbors` (identity included; `< d` branch allows the
  first character to vary, else it is fixed) — the standard Compeau & Pevzner `Neighbors(Pattern, d)` recursion.

### Formula check
- Count_d = count of equal-length windows with HammingDistance ≤ d. Matches BA1I definition.
- Approximate occurrence position predicate matches BA1H definition exactly.
- Frequent-words-with-mismatches tallies each window's full d-neighborhood and returns all maxima — matches BA1I.

### Edge-case semantics
- d=0 degenerates: Neighbors(P,0)={P} ⇒ Count_0 = exact count and FrequentWords = exact frequent k-mer.
  Sourced (BA1I/BA1N) and implemented (`GenerateNeighbors` returns only `pattern` when d==0). PASS.
- Pattern not an exact substring is valid (matching over the Hamming ball). Sourced. PASS.
- Ties: all maxima returned. Sourced (sample has 3). PASS.
- `FindBestMatch` "best single match" and its leftmost tie-break are **not** defined by any
  Rosalind/textbook problem — recorded as an API convention (ASSUMPTION-1). The *distance value* and
  *set of minimal windows* are fully sourced via the Hamming definition; only the choice among equal-distance
  windows is convention. Acceptable as a documented, deterministic convention.

### Independent cross-check (hand computation, not the repo code)
Recomputed with an independent Python script this session:

| Case | Independent value | Source value |
|------|-------------------|--------------|
| BA1H positions | [6, 7, 26, 27, 78], count 5 | 6 7 26 27 78 (BA1H sample) ✓ |
| Count_1 worked example | 4 (AACAA@0, ATAAA@9, AAACA@11, AAAGA@19) | 4 (BA1I) ✓ |
| BA1I frequent words | max 5, {ATGC, ATGT, GATG} | GATG ATGC ATGT (BA1I) ✓ |
| d=0 on AAAAAA, k=4 | {AAAA: 3} | INV-1 degenerate ✓ |
| HD(ACGT, TTTT) | 3 | Hamming def (BA1H) ✓ |
| ACGTACGA vs TTTT (best) | all windows dist 3, leftmost pos 0 | INV-4/INV-5 ✓ |

All sourced expected values reproduced by independent computation.

### Findings / divergences
None. Description is mathematically correct and matches the cited primary sources. Stage A = PASS.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs`
- `FindWithMismatches` (lines 40–87): sliding window, per-window Hamming early-out at `mismatches <= maxMismatches`,
  upper-cases input, yields position + matched window + mismatch positions. Correct O(n·m).
- `CountApproximateOccurrences` (296–299): `.Count()` over `FindWithMismatches`. Correct Count_d.
- `FindFrequentKmersWithMismatches` (312–348): tallies `GenerateNeighbors(window, d)` per window, returns all
  k-mers tying the max. `GenerateNeighbors` (353–385) is the standard recursion, identity included. Correct BA1I.
- `FindBestMatch` (248–288): scans equal-length windows, keeps strictly-smaller distance (so leftmost minimum),
  short-circuits on distance 0. Correct minimum-Hamming with leftmost tie-break.

### Formula realised correctly?
Yes. Verified against hand computation: all six cross-check rows above match the code's behaviour
(the existing tests M1–M8/S1–S3 plus my new tie-break test all pass).

### Cross-verification table recomputed vs code
The full PAT-APPROX-003 fixture (16 tests after the fix) passes, including the BA1I sample (M1/M2),
BA1H count and positions (M3/M4), Count_1 worked example (M5), Count_0 (M6), best-match exact/non-exact
(M7/M8), d=0 degenerate (S1), tie-break (S2 + new S2b), case-insensitivity (S3), and contract cases (C1–C4).

### Variant/delegate consistency
`CountApproximateOccurrences` delegates to `FindWithMismatches`; `FindWithMismatches` DnaSequence overloads
delegate to the string overload — consistent. `FindBestMatch` uses the same `HammingDistance` as the rest.

### Numerical robustness
Integer counts; no overflow on stated ranges (k ≤ 12, d ≤ 3 practical bound). Empty/oversized-pattern guarded.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** every MUST/SHOULD expected value (BA1I set & count 5, BA1H positions & count 5,
  Count_1=4, Count_0=2, d=0 AAAA×3, HD-based distances) is now traced to BA1H/BA1I retrieved this session and
  to an independent hand computation — not to the implementation's output. PASS.
- **No green-washing:** assertions use exact `Is.EqualTo` on values/sets/positions; no Greater/AtLeast/ranges,
  no widened tolerances, no skips. PASS.
- **Defect found — INV-5 tie-break under-tested (test-quality defect):** the original S2 test
  (`ACGTACGA` vs `TTTT`) cannot distinguish "leftmost minimum" from "any minimum": every window has distance 3,
  so the first window scanned (pos 0) is trivially a minimum under the strict `<` comparison — the test would
  still pass against an implementation that returned *any* minimal window. INV-5 (the leftmost convention,
  the one assumption in this unit) was therefore not actually exercised.
  **Fix (in-session):** added `FindBestMatch_LaterTiedMinimum_ReturnsLeftmostMinimumNotFirstWindow`
  using `GAATAAAT` vs `AAAA` (window distances 2,1,1,1,1): position 0 is *not* a minimum and must be skipped,
  positions 1–4 tie at the minimum, so returning position 1 (`AATA`) is correct only if the leftmost minimum
  is chosen. Expected values hand-computed independently. The original S2 was kept (it is correct, just weak).
- **Coverage:** all three canonical methods plus `FindWithMismatches` exercised; all Stage-A branches
  (neighbors/d>0, d=0 degenerate, ties, exact, non-exact, leftmost tie-break) and contract cases
  (empty seq, empty pattern, oversized pattern, invalid k/d) covered.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6536` (baseline 6535 + 1 new test);
  `dotnet build` 0 errors (the 4 warnings are pre-existing in `ApproximateMatcher_EditDistance_Tests.cs`,
  a different unit, untouched here). PASS.

### Findings / defects
- One test-quality defect (INV-5 leftmost tie-break not genuinely tested) — **fixed in-session** by adding a
  discriminating test. No implementation defect; the code already returns the leftmost minimum correctly.

## Verdict & follow-ups
- **Stage A:** PASS — description matches BA1H/BA1I/BA1N primary sources and independent hand computation.
- **Stage B:** PASS-WITH-NOTES — implementation correct; one weak (non-discriminating) tie-break test
  upgraded with a sourced, discriminating case.
- **End-state:** CLEAN — defect completely fixed; full suite green, build clean.
- **Refs:** PAT-APPROX-003
