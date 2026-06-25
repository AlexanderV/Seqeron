# Validation Report: ALIGN-GLOBAL-001 — Global Alignment (Needleman–Wunsch)

- **Validated:** 2026-06-24   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` → `GlobalAlignCore` + `Traceback`; delegates: `GlobalAlign(string,string,ScoringMatrix?)`, `GlobalAlign(string|DnaSequence, …, CancellationToken, IProgress?)`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Needleman–Wunsch algorithm** (fetched live this session). Confirms verbatim:
  - Initialization: `F(0,j) = d·j`, `F(i,0) = d·i`.
  - Recurrence: `F(i,j) = max(F(i−1,j−1)+S(Aᵢ,Bⱼ), F(i,j−1)+d, F(i−1,j)+d)`.
  - Optimal score location: "The entry F(n,m) gives the maximum score among all possible alignments."
  - Worked example: GCATGCG vs GATTACA, match +1 / mismatch −1 / gap −1 → **optimal score 0**; one optimal alignment `GCAT-GCG / G-ATTACA`.
- Evidence doc (`docs/Evidence/ALIGN-GLOBAL-001-Evidence.md`) and TestSpec reproduce these formulas accurately, including the border row `0,−1,−2,…,−7`.

### Formula check
Single linear gap penalty `d` (true global, end gaps penalized). The TestSpec/Evidence map `d = ScoringMatrix.GapExtend`; `GapOpen` is unused by NW. Matches the Wikipedia "Advanced presentation of algorithm" equations exactly. NW (1970) and the standard linear-gap formulation are correctly described.

### Edge-case semantics
- Identical sequences → diagonal, score `n·Match` (sourced, M6).
- Completely different equal-length → all mismatches, negative score (sourced corner case, M5).
- Unequal lengths → border `d·i`/`d·j` drives the optimal (M3/M4).
- Multiple optimal alignments → one returned deterministically; explicitly allowed by source. Note: Wikipedia shows `GCAT-GCG/G-ATTACA` while TestSpec quotes `GCATG-CG/G-ATTACA` — both score 0, both valid optima.

### Independent cross-check (numbers)
Hand-recomputed the DP matrix for GCATGCG/GATTACA. Row 0 = `0 −1 −2 −3 −4 −5 −6 −7`; Row 1 (G) = `−1 1 0 −1 −2 −3 −4 −5`, reproducing the canonical Wikipedia table from the recurrence. Source-stated optimal F(7,7) = 0, confirmed by code (test M1 passes with Score == 0).

### Findings / divergences
None. Description is mathematically faithful to the primary/reference sources.

## Stage B — Implementation

### Code path reviewed
`SequenceAligner.cs`: `GlobalAlignCore` (L220–273) + `Traceback` (L468–538) — canonical; pooled flat DP buffer copied to 2D for traceback. Cancellation overload (L96–202) carries its own inline fill+traceback. `GlobalAlign` DnaSequence/string delegates L58–84, 207–218.

### Formula realised correctly?
Yes. Init `buf[i*cols]=i·GapExtend`, `buf[j]=j·GapExtend` (L239–242); recurrence `max(diag+match/mismatch, up+GapExtend, left+GapExtend)` (L253–259); returned `Score = score[m,n]` (`Traceback` `score[endI,endJ]` with endI=m, endJ=n). The cancellation overload (L119–202) computes the identical matrix and returns `matrix[m,n]`. Both paths are the validated NW linear-gap model with end gaps penalized.

### Cross-verification table recomputed vs code (test run, 13/13 pass)
| ID | Input | Expected | Result |
|----|-------|----------|--------|
| M1 | GCATGCG/GATTACA (+1/−1/−1) | 0 | PASS |
| M3 | T/ACGT, gap −1 | −2 | PASS |
| M4 | ACGT/T, gap −1 | −2 | PASS |
| M5 | AAAA/TTTT | −4 | PASS |
| M6 | ACGTACGT (identical) | 8, no gaps | PASS |
| M7 | ACGT/AGT | 2 | PASS |
| M8 | ACGT/ACGT match +5 | 20 | PASS |
| M10 | ACGT/AGT stats | M3,Mm0,G1,Id75%,Gap25% | PASS |
| M11 | symmetry score(A,B)=score(B,A) | equal | PASS |

### Variant/delegate consistency
- String overload == DnaSequence overload (M9, PASS).
- Cancellation overload == standard overload (`AlignmentProperties.GlobalAlign_CancellationOverload_SameResultAsStandard`, PASS), and same matrix/score logic by inspection.
- Property tests (FsCheck): symmetry, equal aligned length, aligned length ≥ max(input), identity→max score, determinism — all present in `Properties/AlignmentProperties.cs`.

### Numerical robustness
Integer DP; no division in the core. Scores bounded by `len·max(|Match|,|GapExtend|)`; no overflow on stated DNA-length ranges.

### Test quality audit
Assertions check exact sourced scores/alignments and statistics (not just "no throw"); INV-1/2/3 verified by reconstructing originals and recomputing the per-position score. Deterministic. Edge cases (unequal lengths, all-mismatch, identical, single deletion, empty/null) covered.

### Findings / defects
- **Minor (documented, not a defect):** the `DnaSequence` empty-input case is not guarded — `GlobalAlign(new DnaSequence(""), …)` runs the core and returns a score-0 empty-coordinate alignment, whereas the **string** overload returns `AlignmentResult.Empty`. Both are defensible (NW of empty vs empty = 0). The TestSpec frames empty/null handling as API-contract, not algorithm spec. No correctness impact on the NW model.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.**
- No code changes. Full GlobalAlign test class 13/13 green; corroborating property/cancellation tests green.
- No defects logged. Optional future tidy (non-blocking): make the `DnaSequence` empty overload mirror the string overload's `AlignmentResult.Empty` guard for symmetry.
