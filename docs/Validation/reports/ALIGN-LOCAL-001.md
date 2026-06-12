# Validation Report: ALIGN-LOCAL-001 — Local Alignment (Smith–Waterman)

- **Validated:** 2026-06-12   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` (canonical) and `LocalAlign(string, string, ScoringMatrix?)` (delegate), both routing to `LocalAlignCore` → `TracebackLocal` in `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "Smith–Waterman algorithm"** (https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm, accessed 2026-06-12). Confirms the recurrence, the zero floor, first-row/first-column = 0 initialization, max-cell optimum, traceback-to-zero, the linear/affine gap models, and the worked example. The page cites Smith, T.F. & Waterman, M.S. (1981) "Identification of common molecular subsequences", *J. Mol. Biol.* 147:195–197 as the primary reference.

### Formula check
Wikipedia general recurrence (verbatim):

H(i,j) = max{ H(i-1,j-1) + s(a_i,b_j),  max_{k≥1}{H(i-k,j) − W_k},  max_{l≥1}{H(i,j-l) − W_l},  0 }

- **Zero floor:** the `0` term — negative scores are set to 0. This is THE difference from Needleman–Wunsch and the reason the alignment is local.
- **Initialization:** "first row and first column are set to 0" (H(k,0)=H(0,l)=0). No end-gap penalty.
- **Gap model:** linear `W_k = k·W_1` or affine `W_k = u·k + v`. Under the **linear** model the two `max_{k≥1}` terms collapse to the adjacent-cell form, giving the standard:

  **H(i,j) = max( 0, H(i-1,j-1)+s(a_i,b_j), H(i-1,j) − W_1, H(i,j-1) − W_1 )**

- **Optimal local score:** the **maximum cell anywhere in the matrix** (not the corner F(m,n)).
- **Traceback:** starts at that maximum cell; follows the source of each score; **ends at the first cell whose score is 0**.

All match the spec (`tests/TestSpecs/ALIGN-LOCAL-001.md` §1.2, §7) exactly.

### Edge-case semantics
- **No similarity:** every cell = max(0, neg, neg, neg) = 0 ⇒ max score 0, empty alignment. (Sourced from the zero-floor property.)
- **Identical sequences (match M>0, M>|W_1|):** diagonal dominates, H(i,i)=i·M, max = n·M, full alignment, no gaps. (Derivable from the recurrence.)
- **Empty input:** a 0-length sequence leaves only the all-zero initialized row/column ⇒ max score 0. (Sourced from scoring-matrix definition.)

### Independent cross-check (hand computation)
Hand-recomputed the **full DP matrix** for the Wikipedia example A=`TGTTACGG`, B=`GGTTGACTA`, match +3, mismatch −3, linear gap W_1 = 2, recurrence `max(0, diag+s, up−2, left−2)`:

```
       -   G   G   T   T   G   A   C   T   A
  -  [ 0   0   0   0   0   0   0   0   0   0 ]
  T  [ 0   0   0   3   3   1   0   0   3   1 ]
  G  [ 0   3   3   1   1   6   4   2   1   0 ]
  T  [ 0   1   1   6   4   4   3   1   5   3 ]
  T  [ 0   0   0   4   9   7   5   3   4   2 ]
  A  [ 0   0   0   2   7   6  10   8   6   7 ]
  C  [ 0   0   0   0   5   4   8  13  11   9 ]
  G  [ 0   3   3   1   3   8   6  11  10   8 ]
  G  [ 0   3   6   4   2   6   5   9   8   7 ]
```

**Maximum = 13 at H[6,7]** (row C, col C). Traceback (diag → up → left priority):
H[6,7]=13 →diag H[5,6]=10 →diag H[4,5]=7 →left H[4,4]=9 →diag H[3,3]=6 →diag H[2,2]=3 →diag H[1,1]=0 STOP.
**Aligned region:** `GTT-AC` / `GTTGAC`, score **13**; positions start1=1,end1=5, start2=1,end2=6.

This reproduces Wikipedia's stated alignment `GTT-AC / GTTGAC` and matches the spec's §7.1 matrix cell-for-cell. The swapped case (M2) was also hand-computed: max 13 at H[7,6], yielding `GTTGAC / GTT-AC` (gap moves to seq2 via the "up" branch) — matches spec §7.2.

### Findings / divergences
None. The description is biologically and mathematically correct and faithfully sourced.

---

## Stage B — Implementation

### Code path reviewed
`SequenceAligner.cs`: `LocalAlign` (DnaSequence) lines 287–296; `LocalAlign` (string) lines 301–313; `LocalAlignCore` lines 315–348; `TracebackLocal` lines 350–398. `ScoringMatrix` record in `Seqeron.Genomics.Infrastructure/AlignmentTypes.cs:8`.

### Formula realised correctly? (evidence)
- **Zero floor:** line 335 `score[i,j] = Math.Max(0, Math.Max(diag, Math.Max(up, left)))`. Present — this is what makes it local (vs the global path which omits the `0`).
- **Recurrence terms:** diag = `score[i-1,j-1] + matchScore`, up = `score[i-1,j] + GapExtend`, left = `score[i,j-1] + GapExtend` (lines 331–333). Matches the linear-gap S–W recurrence; `matchScore` = `Match` on equal chars else `Mismatch` (line 329).
- **Gap model:** **linear**, penalty = `scoring.GapExtend` per position. `GapOpen` is intentionally unused here (no affine term in S–W local). With the test's `GapExtend = -2` this equals Wikipedia's `W_1 = 2`. Correct.
- **Initialization:** `int[m+1, n+1]` defaults to 0 for all cells, so row 0 and column 0 are 0 with no explicit penalty loop. Correct.
- **Max-cell optimum:** lines 337–342 track `maxScore`/`maxI`/`maxJ` over every filled cell; traceback starts at `(maxI,maxJ)` (line 347). It does **not** take F(m,n). Correct.
- **Score reported:** `score[endI,endJ]` (line 392) = value of the max cell. Correct.
- **Traceback-to-zero:** loop condition `while (i>0 && j>0 && score[i,j] > 0)` (line 359) terminates at the first 0 cell. Correct. Branch priority diag → up → left (lines 366–383) reproduces the spec's expected gap placement.
- **Positions:** `StartPosition*/EndPosition*` are 0-indexed origin/terminus of the traced region (lines 394–397), within bounds (INV-5).

### Cross-verification table recomputed vs code
| Case | Input | Expected (source/hand) | Code result | Match |
|------|-------|------------------------|-------------|-------|
| M1 Wikipedia | TGTTACGG / GGTTGACTA, +3/−3/−2 | score 13, `GTT-AC`/`GTTGAC`, (1,5)/(1,6) | identical | ✅ |
| M2 swapped | GGTTGACTA / TGTTACGG | score 13, `GTTGAC`/`GTT-AC`, (1,6)/(1,5) | identical | ✅ |
| M3 string overload | same as M1 | score 13, `GTT-AC`/`GTTGAC` | identical | ✅ |
| M4 empty | "" / "ACGT" | `AlignmentResult.Empty` | identical | ✅ |
| M5 null | null DnaSequence | `ArgumentNullException` | thrown | ✅ |
| S1 identical | ACGTACGT ×2, +3 | score 24, full match | identical | ✅ |
| S2 dissimilar | AAAA / TTTT | score 0, empty alignment | identical | ✅ |

All M1/M2 values were independently hand-recomputed in Stage A and confirmed against the actual code via the passing test suite.

### Variant/delegate consistency
String overload (`LocalAlign(string,...)`) applies `ToUpperInvariant()` then calls the same `LocalAlignCore`; returns `AlignmentResult.Empty` on empty input. The DnaSequence overload null-guards then calls the identical core. M3 confirms parity.

### Edge cases in code
- No similarity (S2): `maxScore` stays 0, `maxI=maxJ=0`, the `i>0 && j>0` guard fails immediately ⇒ empty aligned strings and `score[0,0]=0`. Confirmed.
- Identical (S1): diagonal dominates, max at H[n,n]=24, full diagonal traceback, no gaps. Confirmed.
- Empty/null handled at the public overloads (M4/M5).

### Test quality audit
Tests assert exact sourced values (score 13, exact aligned strings, exact 0-indexed positions), the zero-floor and AlignmentType invariants, gap-removal-to-substring (INV-4), and both gap-direction traceback branches (left in M1, up in M2). Deterministic; no tautology-only assertions. Edge cases (empty, null, all-match, all-mismatch) covered.

### Findings / defects
None.

---

## Verdict & follow-ups
- **Stage A: PASS** — recurrence (with zero floor), initialization, max-cell optimum, linear gap model, and traceback-to-zero all match Wikipedia / Smith–Waterman (1981); the full DP matrix for the Wikipedia example was hand-recomputed to max score 13 with alignment `GTT-AC`/`GTTGAC`.
- **Stage B: PASS** — the implementation faithfully realises the validated recurrence: zero floor present (line 335), optimum taken from the max cell not F(m,n) (lines 337–347), traceback terminates at 0 (line 359), linear gap via `GapExtend`. Hand-computed example reproduced by the code.
- **End-state: ✅ CLEAN.** No defects found; no code changes required. Build green; LocalAlign filter 12/12 pass; full `Seqeron.Genomics.Tests` suite 4461 passed, 0 failed.
