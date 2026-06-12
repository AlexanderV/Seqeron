# Validation Report: ALIGN-GLOBAL-001 — Global Alignment (Needleman–Wunsch)

- **Validated:** 2026-06-12   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` (canonical, routes to `GlobalAlignCore`); string overload; `(string, string, ScoringMatrix?, CancellationToken, IProgress?)` cancellation overload; `(DnaSequence, …, CancellationToken, …)` overload.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_GlobalAlign_Tests.cs`, `PerformanceExtensionsTests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: "Needleman–Wunsch algorithm"** (fetched 2026-06-12). Confirms verbatim:
  - Recurrence: `F(i,j) = max( F(i−1,j−1)+S(Aᵢ,Bⱼ), F(i,j−1)+d, F(i−1,j)+d )`.
  - Initialization (linear gap penalty `d`): `F(0,j) = d·j`, `F(i,0) = d·i`.
  - Gap model is **LINEAR** (a single uniform penalty `d`); affine (gap-open + gap-extend, Gotoh) is described only as an *extension*, not the base algorithm.
  - Optimal global score is the bottom-right cell `F(n,m)`.
  - Traceback runs from `F(n,m)` back to `F(0,0)`; diagonal = match/mismatch, vertical/horizontal = gap.
  - Worked example: GCATGCG vs GATTACA, match +1 / mismatch −1 / gap −1 → **optimal score 0**.
- **Spec** (`tests/TestSpecs/ALIGN-GLOBAL-001.md`) and **Evidence** (`docs/Evidence/ALIGN-GLOBAL-001-Evidence.md`) cite the same recurrence, initialization, linear gap model, and example. The spec correctly uses **GCATGCG vs GATTACA** (not the variant "GCATGCU vs GATTACA"). The Evidence border row/column (`0,−1,…,−7`) match.

### Formula check
Spec §1.2(3) and Evidence §3 state the recurrence and init exactly as the source. Gap model: linear, with `ScoringMatrix.GapExtend` serving as `d` and `GapOpen` unused by `GlobalAlign`. Confirmed correct.

### Edge-case semantics (sourced)
Two empty → score 0; one empty → all-gap column with linear penalty; identical → diagonal, score n·match; completely different equal length → all mismatches. All consistent with the recurrence and the cited corner cases.

### Independent cross-check — full hand-computed DP matrix
GCATGCG (rows) vs GATTACA (cols), match +1, mismatch −1, gap −1, init `F(i,0)=−i`, `F(0,j)=−j`:

```
       -   G   A   T   T   A   C   A
   -   0  -1  -2  -3  -4  -5  -6  -7
   G  -1   1   0  -1  -2  -3  -4  -5
   C  -2   0   0  -1  -2  -3  -2  -3
   A  -3  -1   1   0  -1  -1  -2  -1
   T  -4  -2   0   2   1   0  -1  -2
   G  -5  -3  -1   1   1   0  -1  -2
   C  -6  -4  -2   0   0   0   1   0
   G  -7  -5  -3  -1  -1  -1   0   0
```

**F(7,7) = 0** — matches spec M1, the Evidence dataset, and Wikipedia. Independently recomputed (Python) and confirmed.

### Findings / divergences
None. Description is mathematically correct and faithfully sourced. Stage A: **PASS**.

---

## Stage B — Implementation

### Code path reviewed
- `GlobalAlign(DnaSequence,…)` → `GlobalAlignCore` (`SequenceAligner.cs:58-67, 220-273`): pooled flat buffer DP.
  - Init: `buf[i*cols]=i*GapExtend`, `buf[j]=j*GapExtend` (`:239-242`) ⇒ `F(i,0)=d·i`, `F(0,j)=d·j`. ✔
  - Fill: `diag=F[i-1,j-1]+match`, `up=F[i-1,j]+d`, `left=F[i,j-1]+d`, `max` of the three (`:255-259`). ✔ Exactly the validated recurrence.
  - Score returned = `score[m,n]` via `Traceback` `endI=m, endJ=n` (`:528`). ✔ Correct optimal cell.
  - Traceback (`:468-534`): diagonal-first tie-break, then up (gap in seq2), then left; deterministic. Matches sourced rules.
- Cancellation overload (`:96-202`): separate inline matrix + traceback, **identical recurrence/init**; verified to produce the same score and alignment as the canonical overload.
- `GapExtend` = linear `d`; `GapOpen` unused by global alignment. ✔ (Matches spec §6.)

### Formula realised correctly? (evidence — actual library output)
Built a consumer against the real assembly:

| Input | Code score | Code alignment | Expected |
|-------|-----------|----------------|----------|
| GCATGCG / GATTACA | **0** | `GCA-TGCG` / `G-ATTACA` (=0 by hand) | 0 ✔ |
| ACGTACGT / ACGTACGT | **8** | no gaps | 8 ✔ |
| ACGT / AGT | **2** | `ACGT` / `A-GT` | 2 ✔ |
| A / A | **1** | `A`/`A` | 1 ✔ |
| A / T | **−1** | `A`/`T` | −1 ✔ |
| "" / GATTACA (string) | **0** | `AlignmentResult.Empty` | Empty ✔ (API guard) |
| cancellation overload, GCATGCG/GATTACA | **0** | `GCA-TGCG`/`G-ATTACA` | identical to canonical ✔ |

The returned alignment `GCA-TGCG`/`G-ATTACA` is a valid optimal alignment (gap-removal yields the originals; per-position sum = +1−1+1−1+1−1+1−1 = 0). It differs from the *particular* alignment shown on Wikipedia, which is expected — multiple optimal tracebacks exist and the source explicitly permits returning any one.

### Cross-verification table recomputed vs code
All MUST cases (M1–M11) and API cases (A1, A2) execute against the built code with the exact sourced values; hand-computed matrix score (0) matches `matrix[m,n]`.

### Variant/delegate consistency
String overload == DnaSequence overload (M9, passing). Cancellation overload == canonical (verified: same score + same aligned strings). Two-empty via DnaSequence overload → score 0, empty alignment (correct NW); string overload returns `Empty` by API guard — both consistent with sourced semantics.

### Test quality audit
Tests assert exact sourced scores (0, −2, −4, 8, 2, 20), exact statistics (M10: matches=3, gaps=1, identity 75%), border-init correctness (M3/M4 distinguish correct linear border from a faulty one), score symmetry (M11), and structural invariants INV-1/2/3 recomputed from the aligned output. Not tautological. Deterministic.

### Findings / defects
None.

---

## Verdict & follow-ups
- **Stage A: PASS** — recurrence, linear gap model, initialization, optimal cell, traceback, and the worked example (score 0) all match Wikipedia / Needleman–Wunsch (1970); full DP matrix recomputed by hand.
- **Stage B: PASS** — implementation realises the validated recurrence exactly; canonical, string, and cancellation overloads agree; all edge cases behave as sourced.
- **End state: CLEAN** — no defect found; no code changes required. `Seqeron.Genomics.Tests`: GlobalAlign filter 42/42 passing; full suite 4461 passed, 0 failed.
