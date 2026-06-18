# Validation Report: RNA-PARTITION-001 — RNA Partition Function (McCaskill) & Boltzmann Structure Probability

- **Validated:** 2026-06-16   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.CalculatePartitionFunction(string, double, double)`, `CalculateStructureProbability(double, double, double)`, `GenerateRandomRna(int[, Random], double)`
- **Stage A verdict:** PASS-WITH-NOTES (one description defect found and fixed)
- **Stage B verdict:** FAIL → FIXED (probability defect corrected this session; full suite green)

## Stage A — Description

### Sources opened & what they confirm (retrieved/derived this session)
- **McCaskill JS (1990)** *Biopolymers* 29:1105-1119 (PMID 1695107): partition function `Z = Σ_S exp(−E(S)/RT)` over pseudoknot-free structures, O(n³) recursive scheme, base-pair binding probabilities over the ensemble. Confirms the model.
- **Will S., MIT 18.417 (mccaskill.pdf / mccaskill2.pdf):** inside recursion `Q_ij = Q_{i,j-1} + Σ_{i≤k<j−m} Q_{i,k-1}·Q^b_{kj}`, base case `Q=1`; `Z=Q_{1n}`; `Pr[(i,j)|S]=Z⁻¹ Σ_{P∋(i,j)} exp(−βE(P))`. The base-pair probability is the **outside** recursion: external term **plus** enclosing-pair terms. The Evidence's earlier reading ("only the external term contributes in the flat model") is **not** what the slides say.
- **Freiburg RNA McCaskill teaching tool:** simplified fixed-per-pair energy `E_bp`, WC + GU pairing, min loop `l`. Matches the implementation's energy model.
- **ViennaRNA pf_fold reference:** `p(s)=exp(−βE)/Z`, `β=1/kT`, `k≈1.987e-3 kcal/(mol·K)`, default 37 °C = 310.15 K. Confirms `CalculateStructureProbability` and constants.

### Formula check
- `Z`, inside recursion, base case, `Z=Q_{1n}`: **correct**, matches Source 2 verbatim.
- Boltzmann structure probability `p=exp(−βE)/Z`: **correct** (Source 1, 5). RT = 1.987·310.15/1000 = 0.61626805 kcal/mol; exp(−1/RT) = 0.197370910785… (independently recomputed).
- **DEFECT (fixed):** the base-pair-probability formula in the Evidence + algorithm docs stated that in the flat fixed-per-pair model `p_kl = p^E_kl` (external decomposition only). This is **mathematically false**: a pair that can be nested inside another pair receives strictly more probability from the enclosing structures. The correct formula is the McCaskill **outside recursion** `P[i,j] = Q^b(i,j)·O(i,j)/Z` with `O(i,j) = Q(0,i-1)·Q(j+1,n-1) + Σ_{k<i,l>j,CanPair,l-k>m} w·Q(k+1,i-1)·Q(j+1,l-1)·O(k,l)`, `w = exp(−βE_bp)`. Docs corrected this session.

### Edge-case semantics
`null`→`ArgumentNullException`; empty→`Z=1` no pairs; no admissible pair (`AAAA`)→`Z=1`; span ≤ m (`GC`)→`Z=1`; `temperature ≤ 0`→`ArgumentOutOfRangeException`. All defined and sourced (base case `Q=1`; min-loop `j-i>m`). **Correct.**

### Independent cross-check (numbers, this session)
Two independent re-implementations, both separate from the library:
1. **Exhaustive brute-force enumeration** of all non-crossing base-disjoint pair subsets (Boltzmann-weighted) — `/tmp/rna_bruteforce.py`, `/tmp/bf_weighted.py`.
2. **Standalone McCaskill Q/Q^b + outside recursion** — `/tmp/mccaskill_recur.py`, `/tmp/outside3.py`.

| Sequence | E_bp | Z | Key probabilities (true, by enumeration) |
|----------|------|---|------------------------------------------|
| AAAA | 0 | 1 | — |
| GC | 0 | 1 | — |
| GAAAAC | 0 | 2 | P[0,5]=0.5 |
| GGGGCCCC | 0 | 16 | P[1,5]=P[2,6]=3/16; P[0,7]=4/16 |
| GGGAAACCC | 0 | 20 | P[0,8]=P[2,6]=6/20; P[1,7]=4/20 |
| GGGGCCCC | −1 | 180.0183448 | P[1,5]=P[2,6]=0.31334323; P[0,7]=0.45594238 |

Brute force vs outside recursion agree to **3.3e-16** over 300 random sequences (len 0–12, E_bp ∈ {0,−1,−2,0.5}); per-base pairing-probability sum ≤ 1 always (max 0.983).

### Findings / divergences
- **D1 (Stage A, fixed):** Evidence §"mccaskill2" point 2 and algorithm-doc §2.2/§5.2/§5.3 claimed external decomposition suffices. Corrected to the outside recursion with the exact datasets above.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`:
- `CalculatePartitionFunction` (≈L1922) — inside DP + probabilities.
- `CalculateStructureProbability` (≈L2010), `GenerateRandomRna` (≈L2024/2032).
- Constants: `GasConstant_CalPerMolK=1.987`, `DefaultTemperatureKelvin=310.15`, `MinHairpinLoop=3` (`j-i>MinHairpinLoop`).

### Formula realised correctly?
- **Z / inside recursion:** correct. Verified the library reproduces 1, 1, 2, 16, 20 (E_bp=0) and 180.0183448 (E_bp=−1) exactly (ran the library directly).
- **DEFECT (found & fixed):** `BasePairProbabilities` originally used only `GetQ(0,i-1)·qb[i,j]·GetQ(j+1,n-1)/z` (external decomposition). Confirmed wrong by running the library: `GGGAAACCC` P[2,6] returned **0.05** (true 0.30), P[1,7] returned 0.10 (true 0.20); `GGGGCCCC` P[1,5]/P[2,6] returned 0.0625 (true 0.1875). Replaced with the outside recursion (external + enclosing terms, pairs processed outermost-first). After the fix the library matches brute force to machine precision for both E_bp=0 and E_bp=−1 and all pairs.
- `CalculateStructureProbability`: returns `exp(−E_s/RT)/exp(−E_ens/RT)`; equals 1 when energies equal, exp(−1/RT)=0.19737… for (−5,−6). Correct.
- `GenerateRandomRna`: seeded determinism, length, alphabet {A,C,G,U}. Correct.

### Cross-verification table recomputed vs code (post-fix)
All values in the Stage-A table above now reproduced exactly by the library (`/tmp/pfverify`). Z unchanged by the fix (only probabilities were affected).

### Variant/delegate consistency
`GenerateRandomRna(int,double)` forwards to the seeded overload — consistent. `CalculateStructureProbability` default temperature 310.15 K consistent with `CalculatePartitionFunction`.

### Test quality audit (HARD gate)
- **Before:** M5/M6 only asserted probabilities for **outermost (never-enclosed)** pairs — (0,5) in GAAAAC; (0,6),(0,7),(0,8) in GGGAAACCC. These pass even with the buggy external-only formula, so the tests were a **code-echoing blind spot** that hid the defect. This is a Stage-B test defect.
- **Fix:** rewrote M6 to assert **all 9** GGGAAACCC pairs (incl. nested (1,6),(1,7),(1,8),(2,6),(2,7)); added **M6b** (all 10 GGGGCCCC pairs incl. nested (1,5),(2,6)); **M6c** weighted E_bp=−1 probabilities (locks the general w≠1 case where the naive `/Q^b(k,l)` shortcut also fails); **M6d** per-base pairing-probability ≤ 1 invariant. Every expected value traces to the independent brute-force enumeration, not to library output. No assertions weakened, no tolerances widened (only one Z literal made full-precision), nothing skipped.
- **Honest green:** full unfiltered suite `dotnet test` → **Passed 6596, Failed 0, Skipped 0**; `dotnet build` 0 errors (4 pre-existing warnings in unrelated files). **Gate: PASS.**

### Findings / defects
- **D2 (Stage B, FIXED):** external-only base-pair probabilities — see register **B5**. Code corrected to the outside recursion; strict nested-pair + weighted + per-base-≤1 tests added.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES — model/formulae correct; one false probability claim in the docs corrected (D1).
- **Stage B:** FAIL → **FIXED** — probability defect (D2) corrected; Z and the Boltzmann probability were already correct.
- **End-state:** ✅ **CLEAN** — code fixed, docs corrected, strict sourced tests added, full suite green (6596/0).
- **Logged:** Findings register B5; ledger row 79.
