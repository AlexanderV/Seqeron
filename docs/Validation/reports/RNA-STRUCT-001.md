# Validation Report: RNA-STRUCT-001 — RNA Secondary Structure Prediction (Nussinov base-pair maximization)

- **Validated:** 2026-06-12   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.PredictStructure(...)`,
  `RnaSecondaryStructure.CalculateMinimumFreeEnergy(...)` (Zuker MFE),
  `RnaSecondaryStructure.CalculateMinimumFreeEnergyClassic(...)` (Nussinov-style weighted-pair DP),
  helpers `CanPair` / `GetBasePairType` / `PairType` / `GetPairEnergy`,
  `ToDotBracket`/`FromDotBracket`(`GenerateFullDotBracket`/`ParseDotBracket`), `ValidateDotBracket`, `DetectPseudoknots`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Nussinov algorithm"** (fetched). Confirms the DP recurrence in its compact two-case
  form and the scoring/traceback:
  - `M(i,j) = max( M(i,j−1) ,  max_{i≤k<j} { M(i,k−1) + M(k+1,j−1) + Score(S_k,S_j) } )`
  - `Score(S_k,S_j) = 1` if complementary else `0`.
  - Traceback: stop when `j ≤ i`; if `M(i,j)=M(i,j−1)` decrement `j`; else find the `k` that pairs with
    `j` and recover pair `(k,j)`, then recurse on `(i,k−1)` and `(k+1,j−1)`.
  This two-case form is mathematically equivalent to the four-case textbook form
  `M(i,j)=max[ M(i+1,j), M(i,j−1), M(i+1,j−1)+δ(i,j), max_k M(i,k)+M(k+1,j) ]`:
  bifurcation is absorbed into the `M(i,k−1)+M(k+1,j−1)` term (when `k=i`, `M(i,i−1)=0`, giving the plain pair case).
- **Wikipedia — "Nucleic acid structure prediction"** (fetched). Confirms Nussinov is an early
  *maximum-base-pairing* DP, O(n³), cannot find pseudoknots; and confirms the biological hairpin
  minimum-loop constraint (≥3 unpaired nucleotides).
- **Nussinov & Jacobson (1980), PNAS 77(11):6309–6313** and **Zuker & Stiegler (1981), NAR 9(1):133–148**
  cited in the Evidence doc as the primary references for, respectively, the pair-maximization DP and
  the thermodynamic (MFE) DP. The repository correctly treats these as two different algorithm classes.

### Formula check
- **Allowed pairs:** Watson–Crick A-U / U-A / G-C / C-G **and** wobble G-U / U-G. The spec/evidence
  include G-U wobble (standard for RNA); confirmed correct (Wikipedia "Nucleic acid secondary structure").
- **Recurrence:** matches the canonical Nussinov recurrence (two-case form above) exactly.
- **Minimum loop:** spec uses minimum hairpin loop size = 3 unpaired nucleotides (a pair `(i,j)` requires
  `j − i ≥ 4`). Confirmed steric constraint per NNDB Turner 2004 and Wikipedia "Stem-loop".
- **Nussinov ≠ MFE:** the description explicitly states the weighted-pair DP "maximizes weighted pair
  count (not thermodynamic MFE)" and keeps a *separate* Zuker DP (`CalculateMinimumFreeEnergy`) for
  physical kcal/mol. This distinction is correct and important (Evidence §3.1, D5).

### Edge-case semantics
Empty/null → empty structure, value 0; poly-A (no complement) → 0 pairs / all dots; sequence too short to
admit any pair (`len < minLoop+2`) → 0; minimum loop enforced (no sharp hairpins). All defined and sourced.

### Independent cross-check (hand-computed Nussinov DP)
Sequence **GGGAAAUCC** (0–8), pairs A-U/G-C/G-U, `j−i ≥ 4`:
- Feasible pairs include G0-C7/C8, G1-C7/C8, G2-C7/C8, G·-U6, A3-U6, A4-U6.
- Maximum nested non-crossing set: **(0,8) G-C, (1,7) G-C, (2,6) G-U → 3 base pairs**; the innermost
  window (3,5) has only 3 bases so no further pair is possible (would need `j−i ≥ 4`).
- **Max base pairs = 3**, valid dot-bracket **`(((...)))`**.
Weighted score (WC −2, wobble −1): −2 −2 −1 = **−5.0**.

### Findings / divergences
None. Description is biologically and mathematically correct.

---

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`
- Pair set: `BuildPairLookup` (l.412) and `PairType` (l.1381) — A-U/U-A/G-C/C-G = 1 (WC), G-U/U-G = 2 (wobble). Correct, includes wobble.
- Nussinov-style weighted-pair DP: `CalculateMinimumFreeEnergyClassic` (l.1394):
  `dp[i,j]=dp[i,j−1]` (j unpaired); inner `for k=i; k < j − minLoopSize; k++` with
  `dp[i,k−1] + dp[k+1,j−1] + pairEnergy(k,j)` when `CanPair(k,j)`. This is exactly the canonical Nussinov
  recurrence; the `k < j − minLoopSize` bound enforces `j − k ≥ 4` (≥3 unpaired) — minimum loop correct.
  Bifurcation is realized by the `dp[i,k−1]` left term (general two-case Nussinov). `GetPairEnergy` = −2 WC / −1 wobble.
- Zuker MFE: `CalculateMinimumFreeEnergy` (l.1021) — separate V/W/WM recurrence with Turner 2004 NNDB
  parameters; correctly returns physical kcal/mol and may be 0 when no thermodynamically stable fold exists.
- `PredictStructure` (l.1452): builds the dot-bracket by `FindStemLoops` + greedy non-overlapping
  selection (`SelectNonOverlapping`), **not** a Nussinov DP traceback (deviation D5).
- Dot-bracket: `GenerateFullDotBracket`/`ParseDotBracket`/`ValidateDotBracket`/`DetectPseudoknots` — standard, balanced.

### Cross-verification table recomputed vs code (verified by scratch run, then removed)
| Sequence | Nussinov (classic) | Max pairs | Pred dot-bracket | Valid | Zuker MFE |
|----------|-------------------|-----------|------------------|-------|-----------|
| GGGAAAUCC | −5.0 | 3 (0,8)(1,7)(2,6) | `(((...)))` | yes | 0 (no stable fold) |
| GGGGAAAACCCC | −8.0 | 4 | `((((....))))` | yes | −5.28 |
| GGGAAACCC | −6.0 | 3 | `(((...)))` | yes | −1.12 |
| AAAAAAA | 0 | 0 | `.......` | yes | 0 |
| GC (too short) | 0 | 0 | `..` | yes | 0 |
| "" / null | 0 | 0 | `` | yes | 0 |
| GGGAAACCCCCCAAAGGG (bifurcation) | −12.0 | 6 | `(((......)))......` | yes | −3.72 |

All match the hand computation. Zuker values GGGAAACCC=−1.12 and GGGGAAAACCCC=−5.28 match the spec's
cited Turner-2004 manual references. The bifurcation case (−12, two hairpins) confirms the bifurcation
term works. Edge cases (no-pair → 0/all-dots, min-loop enforced so no sharp hairpins, single/short
sequence → 0) all hold.

### Variant/delegate consistency
`CalculateMinimumFreeEnergyClassic` (Nussinov) and `CalculateMinimumFreeEnergy` (Zuker) are deliberately
distinct value functions and agree on which sequences are structurable. `PairType`/`PairLookup` agree.

### Test quality audit
137 RnaSecondaryStructure tests pass; they assert exact NNDB Turner values, validity invariants
(`ValidateDotBracket` true, length = sequence length), WC/wobble-only pairs, non-overlap, empty/poly-A
edge cases, and crossing/non-crossing pseudoknots. Assertions are real (exact values, not "no-throw").

### Findings / defects
No correctness defect in the Nussinov DP, the pair set, the minimum-loop constraint, the bifurcation term,
or the dot-bracket traceback for the maximization variants. One **documented limitation (D5)**:
`PredictStructure` derives its dot-bracket from greedy non-overlapping stem-loop selection rather than a
Nussinov/Zuker DP traceback, so its predicted structure is a valid (always balanced, non-crossing)
approximation that may not be the global optimum. This is sourced and intended, not a bug; the optimal
*pair count* (Nussinov) and optimal *energy* (Zuker) values are themselves computed correctly by the DPs.

---

## Verdict & follow-ups
- **Stage A: PASS** — Nussinov recurrence, pair set (incl. G-U wobble), minimum loop = 3, and the
  Nussinov-vs-MFE distinction are all correct against Nussinov & Jacobson (1980), Zuker & Stiegler (1981),
  and Wikipedia.
- **Stage B: PASS-WITH-NOTES** — DP recurrence, pair set, min-loop, bifurcation, and dot-bracket all
  validated; hand-computed example (GGGAAAUCC → 3 pairs, `(((...)))`, score −5) reproduced exactly by the
  code. Note D5: `PredictStructure` uses greedy selection, not DP traceback (documented approximation).
- **State: CLEAN** — no defect; no code change required. Full suite green (4486 passed).
- No new defects logged.
