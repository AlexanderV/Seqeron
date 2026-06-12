# Validation Report: PAT-APPROX-002 — Approximate Matching (Edit / Levenshtein Distance)

- **Validated:** 2026-06-12   **Area:** Pattern Matching
- **Canonical method(s):** `ApproximateMatcher.EditDistance(string, string)`, `ApproximateMatcher.FindWithEdits(string, string, int)`, `FindWithEdits(DnaSequence, string, int)` (wrapper)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Levenshtein distance"** (fetched): confirms the exact recurrence and base cases (see below); confirms each of insertion/deletion/substitution costs **1**, match costs **0**; confirms worked examples `kitten→sitting = 3` (k→s sub, e→i sub, +g insert) and `flaw→lawn = 2` (delete f, insert n); confirms transpositions are **NOT** free in standard Levenshtein (that is Damerau–Levenshtein, which the implementation does not claim).
- **Rosetta Code — "Levenshtein distance"**: page returned HTTP 403; its canonical vectors were instead independently reproduced with a from-scratch reference DP (below), which is a stronger check than trusting the citation label.
- **Navarro (2001), "A Guided Tour to Approximate String Matching"**: standard edit-distance foundation; consistent with unit-cost Levenshtein. (Not load-bearing beyond the recurrence already confirmed by Wikipedia.)

### Formula check (exact recurrence verified)
Base cases: `lev(i,0) = i`, `lev(0,j) = j`.
Recurrence (i,j ≥ 1):
```
lev(i,j) = min(
    lev(i-1, j)   + 1,                       // deletion
    lev(i,   j-1) + 1,                        // insertion
    lev(i-1, j-1) + (a[i-1]==b[j-1] ? 0 : 1) // substitution / match
)
```
Substitution cost = **1** (Levenshtein), not 2 (LCS-style). Match cost = 0. Confirmed verbatim against Wikipedia.

### Edge-case semantics check
- empty vs empty = 0 ✔ (both base cases give 0)
- empty vs length-n = n ✔ (`lev(0,j)=j`)
- identical = 0 ✔
- transposition is NOT free (e.g. `stop→tops = 2`, not 1) ✔
- case-sensitivity: standard definition treats characters as distinct symbols (`s[i]=t[j]`, no normalisation), so `A` vs `a` = 1 ✔

### Independent cross-check (numbers)
A from-scratch reference Wagner–Fischer DP (cost 1 ins/del/sub) reproduced **every** spec vector exactly:

| s1 | s2 | reference | spec |
|----|----|-----------|------|
| kitten | sitting | 3 | 3 |
| rosettacode | raisethysword | 8 | 8 |
| saturday | sunday | 3 | 3 |
| stop | tops | 2 | 2 |
| sleep | fleeting | 5 | 5 |
| flaw | lawn | 2 | 2 |
| "" | abc | 3 | 3 |
| abc | "" | 3 | 3 |
| "" | "" | 0 | 0 |

### Findings / divergences
None. The spec's recurrence, costs, edge conventions, invariants (symmetry, identity, empty-string, triangle inequality, `|len a − len b| ≤ d ≤ max(len)`) and all eight canonical vectors are correct and sourced.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs`
- `EditDistance` — lines 190–231
- `FindWithEdits(string,…)` — lines 113–150; `FindWithEdits(DnaSequence,…)` wrapper — lines 155–159

### Formula realised correctly? (evidence)
`EditDistance` is a two-row (rolling) Wagner–Fischer DP:
- base row init `prev[j]=j` (line 207–208) and `curr[0]=i` (line 212) → matches `lev(0,j)=j`, `lev(i,0)=i`;
- `cost = s1[i-1]==s2[j-1] ? 0 : 1` (line 216) → match 0 / mismatch 1;
- `min(prev[j]+1 /*del*/, curr[j-1]+1 /*ins*/, prev[j-1]+cost /*sub*/)` (lines 217–223) → exactly the validated recurrence.
- empty-string fast paths `if (m==0) return n; if (n==0) return m;` (lines 199–200) match the base cases.
- `EditDistance` does **not** uppercase its inputs, so it is correctly case-sensitive per the standard definition (distinct from `FindWithEdits`, which uppercases the search text/pattern by design for DNA matching — a separate, intentional concern).
- Null inputs throw `ArgumentNullException` (lines 192–193); two-row buffer never overflows; `int` cannot overflow on any realistic length (distance ≤ max length).

### Cross-verification table recomputed vs code
The unit test suite asserts the exact eight sourced vectors and all five invariants; all match the Stage-A reference numbers and pass (see Test quality / Tests below).

### Variant/delegate consistency
`FindWithEdits(DnaSequence,…)` simply forwards `sequence.Sequence` to the string overload — verified by test C01 (`FindWithEdits_DnaSequenceOverload_DelegatesToStringVersion`). `FindWithEdits` uses `EditDistance` as its window scorer (line 134) and reports `maxEdits ≥ 0` guard (lines 119–120).

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/ApproximateMatcher_EditDistance_Tests.cs` — 25 tests (M01–M13, S01–S10, C01). Assertions check **exact** sourced values via `Is.EqualTo(<constant>)` with the actual computed value as the first argument — the correct NUnit2007 ordering, so the file builds with **0 warnings** (the NUnit2007 concern in the task brief does not apply to the current file state; no swap needed). Tests cover every Stage-A edge case (identical, both-empty, empty/non-empty, single ins/del/sub, case-sensitivity, transposition-like, bounds, triangle inequality) and are deterministic.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — recurrence, unit costs, base cases, edge semantics, and all eight canonical vectors independently confirmed.
- **Stage B: PASS** — code realises the validated recurrence exactly; case-sensitive `EditDistance`; wrapper consistent; tests assert exact sourced values.
- **End-state: CLEAN** — no defect found. No code or test changes required.
- Build: `Seqeron.Genomics.Tests` builds with 0 warnings / 0 errors.
- Tests: `--filter FullyQualifiedName~EditDistance` → 35 passed / 0 failed; full suite → **4461 passed, 0 failed** (matches baseline).
