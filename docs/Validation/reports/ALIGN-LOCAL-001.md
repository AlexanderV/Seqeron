# Validation Report: ALIGN-LOCAL-001 — Local Alignment (Smith–Waterman)

- **Validated:** 2026-06-24   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` (canonical) and `LocalAlign(string, string, ScoringMatrix?)` (delegate), both routing to `LocalAlignCore` → `TracebackLocal` in `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

This is a re-validation in fresh context. The source file was last functionally touched for local alignment by the e19a8a02 ("O(n) traceback") perf work; the local core/traceback path reviewed below is logically identical to the previously-validated version. No defects found; no code changes made.

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "Smith–Waterman algorithm"** (https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm). Confirms: the local-alignment recurrence with the **zero floor**; first row/first column initialized to 0 (no end-gap penalty); the optimal score is the **maximum cell anywhere in the matrix** (not the corner); traceback starts at that max cell and **ends at the first cell with score 0**; linear (`W_k = k·W_1`) and affine (Gotoh) gap models; and the worked `TGTTACGG` / `GGTTGACTA` example (match +3, mismatch −3, gap 2 → score 13, alignment `GTT-AC` / `GTTGAC`). The page cites Smith, T.F. & Waterman, M.S. (1981), *J. Mol. Biol.* 147:195–197 as the primary reference.

### Formula check
Linear-gap S–W recurrence (Wikipedia §Linear gap penalty), matching TestSpec §1.2 point 5:

**H(i,j) = max( 0, H(i−1,j−1)+s(a_i,b_j), H(i−1,j) − W_1, H(i,j−1) − W_1 )**

- Zero floor (the `0` term) is the distinguishing feature vs Needleman–Wunsch and makes the alignment local.
- Initialization H(k,0)=H(0,l)=0.
- Optimum = max cell anywhere; traceback to first 0.

All match the spec (§1.2, §7) exactly.

### Edge-case semantics
- **No similarity:** every cell = max(0, neg, neg, neg) = 0 ⇒ max score 0, empty alignment (zero-floor property — sourced).
- **Identical sequences (match M>0, M>|W_1|):** diagonal dominates, H(i,i)=i·M, max = n·M, full alignment, no gaps (derivable).
- **Empty input:** a 0-length sequence leaves only the all-zero initialized row/column ⇒ max score 0 (sourced from scoring-matrix definition).

### Independent cross-check (hand computation)
Hand-recomputed the **full DP matrix** for the Wikipedia example A=`TGTTACGG`, B=`GGTTGACTA`, match +3, mismatch −3, linear gap W_1 = 2, recurrence `max(0, diag+s, up−2, left−2)` — reproduces the spec §7.1 matrix cell-for-cell:

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

**Maximum = 13 at H[6,7]** (row C = A[5], col C = B[6]; match → diag = H[5,6]+3 = 10+3 = 13). Traceback (diag → up → left priority):
H[6,7]=13 →diag H[5,6]=10 →diag H[4,5]=7 →**left** H[4,4]=9 →diag H[3,3]=6 →diag H[2,2]=3 →diag H[1,1]=0 STOP.
**Aligned region:** A[1..5]=`GTTAC` / B[1..6]=`GTTGAC` → `GTT-AC` / `GTTGAC`, score **13**; positions start1=1,end1=5, start2=1,end2=6. Matches Wikipedia's stated alignment exactly.

Swapped case (M2) hand-checked against spec §7.2: transposed matrix has max 13 at H^T[7,6], yielding `GTTGAC` / `GTT-AC` (gap moves to seq2 via the "up" branch).

### Findings / divergences
None. The description is biologically and mathematically correct and faithfully sourced.

---

## Stage B — Implementation

### Code path reviewed
`SequenceAligner.cs`: `LocalAlign(DnaSequence)` lines 287–296; `LocalAlign(string)` lines 301–313; `LocalAlignCore` lines 315–348; `TracebackLocal` lines 350–398. `ScoringMatrix` default `SimpleDna` at lines 22–26; gap penalty supplied via `GapExtend` (stored as a negative value, added to the predecessor).

### Formula realised correctly? (evidence)
- **Zero floor:** line 335 `score[i,j] = Math.Max(0, Math.Max(diag, Math.Max(up, left)))`. Present.
- **Recurrence terms:** diag = `score[i-1,j-1] + matchScore`, up = `score[i-1,j] + GapExtend`, left = `score[i,j-1] + GapExtend` (lines 331–333); `matchScore` = `Match` on equal chars else `Mismatch` (line 329). Linear gap = `GapExtend` per position; with the test's `GapExtend = -2` this equals Wikipedia's `W_1 = 2`. `GapOpen` is intentionally unused (no affine term). Correct.
- **Initialization:** `int[m+1, n+1]` defaults to 0 ⇒ row 0 / column 0 are 0 with no penalty loop. Correct.
- **Max-cell optimum:** lines 337–342 track `maxScore`/`maxI`/`maxJ` over every filled cell; traceback starts at `(maxI,maxJ)` (line 347). Does **not** take F(m,n). Correct.
- **Score reported:** `score[endI,endJ]` (line 392) = the max cell value. Correct.
- **Traceback-to-zero:** loop condition `while (i>0 && j>0 && score[i,j] > 0)` (line 359) terminates at the first 0 cell. Branch priority diag → up → left (lines 366–383) reproduces the spec's gap placement. Correct.
- **Positions:** 0-indexed origin/terminus of the traced region (lines 394–397), within bounds (INV-5).

### Cross-verification table recomputed vs code (tests run)
| Case | Input | Expected (source/hand) | Code result | Match |
|------|-------|------------------------|-------------|-------|
| M1 Wikipedia | TGTTACGG / GGTTGACTA, +3/−3/−2 | score 13, `GTT-AC`/`GTTGAC`, (1,5)/(1,6) | identical | ✅ |
| M2 swapped | GGTTGACTA / TGTTACGG | score 13, `GTTGAC`/`GTT-AC`, (1,6)/(1,5) | identical | ✅ |
| M3 string overload | same as M1 | score 13, `GTT-AC`/`GTTGAC` | identical | ✅ |
| M4 empty | "" / "ACGT" | `AlignmentResult.Empty` | identical | ✅ |
| M5 null | null DnaSequence | `ArgumentNullException` | thrown | ✅ |
| S1 identical | ACGTACGT ×2, +3 | score 24, full match | identical | ✅ |
| S2 dissimilar | AAAA / TTTT | score 0, empty alignment | identical | ✅ |

M1/M2 hand-recomputed in Stage A; all 7 confirmed by running `SequenceAligner_LocalAlign_Tests` (7/7 passed).

### Variant/delegate consistency
String overload applies `ToUpperInvariant()` then calls the same `LocalAlignCore`; returns `AlignmentResult.Empty` on empty/null-or-empty input (lines 306–312). DnaSequence overload null-guards then calls the identical core. M3 confirms parity.

### Test quality audit
Tests assert exact sourced values (score 13, exact aligned strings, exact 0-indexed positions), the zero-floor and AlignmentType invariants, gap-removal-to-substring (INV-4), and both gap-direction traceback branches (left in M1, up in M2). Deterministic; no tautology-only assertions. Edge cases (empty, null, all-match, all-mismatch) covered.

### Findings / defects
None.

---

## Verdict & follow-ups
- **Stage A: PASS** — recurrence (with zero floor), initialization, max-cell optimum, linear gap model, and traceback-to-zero all match Wikipedia / Smith–Waterman (1981); the full Wikipedia DP matrix was hand-recomputed to max score 13 with alignment `GTT-AC`/`GTTGAC`.
- **Stage B: PASS** — the implementation faithfully realises the validated recurrence: zero floor (line 335), optimum from the max cell not F(m,n) (lines 337–347), traceback terminates at 0 (line 359), linear gap via `GapExtend`. Hand-computed example reproduced by the code.
- **End-state: ✅ CLEAN.** No defects found; no code changes required. Build green; `SequenceAligner_LocalAlign_Tests` 7/7 pass.
