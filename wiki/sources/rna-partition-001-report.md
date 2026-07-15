---
type: source
title: "Validation report: RNA-PARTITION-001 (RNA partition function — McCaskill & Boltzmann structure probability)"
tags: [validation, rna, governance]
doc_path: docs/Validation/reports/RNA-PARTITION-001.md
sources:
  - docs/Validation/reports/RNA-PARTITION-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: RNA-PARTITION-001

The two-stage **validation write-up** for test unit **RNA-PARTITION-001** (RNA Partition
Function — McCaskill 1990 — and Boltzmann Structure Probability), validated 2026-06-16 in the
RnaStructure area. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's **verdict** on both the algorithm
description and the shipped code. The algorithm itself is summarized in
[[rna-partition-function-mccaskill]]; the two-stage methodology is the
[[validation-protocol]]. Distinct from the pre-implementation [[rna-partition-001-evidence]]
artifact.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: FAIL → FIXED · State: ✅ CLEAN.** Full unfiltered suite
**6596 passed / 0 failed / 0 skipped**, `dotnet build` 0 errors (4 pre-existing warnings in
unrelated files). One description defect (D1) was found and the docs corrected; one code
defect (D2, wrong base-pair probabilities) was found and fixed. `Z` and the Boltzmann
structure probability were already correct; only the base-pair probabilities were affected.

Canonical methods: `RnaSecondaryStructure.CalculatePartitionFunction(string, double, double)`,
`CalculateStructureProbability(double, double, double)`, `GenerateRandomRna(int[, Random],
double)`. See report at `docs/Validation/reports/RNA-PARTITION-001.md`. Logged: findings
register B5, ledger row 79.

## Stage A — description (algorithm faithfulness)

- **Model source-confirmed.** McCaskill (1990) *Biopolymers* 29:1105-1119 (PMID 1695107)
  confirms `Z = Σ_S exp(−E(S)/RT)` over pseudoknot-free structures, O(n³) scheme, and
  ensemble base-pair probabilities. Will S., MIT 18.417 slides confirm the inside recursion
  `Q_ij = Q_{i,j-1} + Σ Q_{i,k-1}·Q^b_{kj}`, base case `Q=1`, `Z=Q_{1n}`. Freiburg RNA teaching
  tool matches the simplified fixed-per-pair `E_bp` energy model; ViennaRNA `pf_fold` confirms
  `p(s)=exp(−βE)/Z`, `k≈1.987e-3 kcal/(mol·K)`, default 37 °C = 310.15 K.
- **Formula check:** `Z`, inside recursion, base case, `Z=Q_{1n}` all correct; Boltzmann
  probability `p=exp(−βE)/Z` correct. RT = 1.987·310.15/1000 = 0.61626805 kcal/mol;
  exp(−1/RT) = 0.19737… (independently recomputed).
- **DEFECT D1 (Stage A, fixed):** the base-pair-probability formula in the Evidence + algorithm
  docs claimed that in the flat fixed-per-pair model `p_kl = p^E_kl` (external decomposition
  only). This is mathematically false — a pair nestable inside another gets strictly more
  probability. The correct formula is the McCaskill **outside recursion**
  `P[i,j] = Q^b(i,j)·O(i,j)/Z`. Docs corrected this session.
- **Edge cases** sourced and correct: `null`→`ArgumentNullException`; empty/`AAAA`/`GC`
  (span ≤ m)→`Z=1`; `temperature ≤ 0`→`ArgumentOutOfRangeException`.
- **Independent cross-check:** two separate re-implementations (exhaustive brute-force
  Boltzmann enumeration; standalone McCaskill Q/Q^b + outside recursion) agree to **3.3e-16**
  over 300 random sequences; per-base pairing sum ≤ 1 always (max 0.983). Oracles: `Z` = 1, 1,
  2, 16, 20 (E_bp=0) and 180.0183448 (E_bp=−1); `GGGAAACCC` P[2,6]=6/20, P[1,7]=4/20.

## Stage B — implementation (code review + cross-check)

- **Code path:** `RnaSecondaryStructure.cs` — `CalculatePartitionFunction` (≈L1922, inside DP
  + probabilities), `CalculateStructureProbability` (≈L2010), `GenerateRandomRna` (≈L2024/2032).
  Constants `GasConstant_CalPerMolK=1.987`, `DefaultTemperatureKelvin=310.15`,
  `MinHairpinLoop=3`.
- **Z / inside recursion correct:** library reproduces 1, 1, 2, 16, 20 (E_bp=0) and
  180.0183448 (E_bp=−1) exactly.
- **DEFECT D2 (Stage B, found & FIXED):** `BasePairProbabilities` originally used only the
  external decomposition `GetQ(0,i-1)·qb[i,j]·GetQ(j+1,n-1)/z`. Confirmed wrong by running the
  library: `GGGAAACCC` P[2,6] returned 0.05 (true 0.30), P[1,7] 0.10 (true 0.20); `GGGGCCCC`
  P[1,5]/P[2,6] returned 0.0625 (true 0.1875). Replaced with the outside recursion (external +
  enclosing terms, pairs processed outermost-first). After the fix the library matches brute
  force to machine precision for E_bp=0 and E_bp=−1 and all pairs.
- **`CalculateStructureProbability`** returns `exp(−E_s/RT)/exp(−E_ens/RT)`; equals 1 when
  energies equal, exp(−1/RT)=0.19737… for (−5,−6). Correct. **`GenerateRandomRna`** seeded
  determinism, length, alphabet {A,C,G,U}, forwards to the seeded overload. Correct.
- **Test-quality audit (HARD gate) PASS:** before the fix, M5/M6 only asserted probabilities
  for outermost (never-enclosed) pairs, which pass even with the buggy external-only formula —
  a code-echoing blind spot that hid the defect. Fix rewrote M6 to assert all 9 `GGGAAACCC`
  pairs incl. nested; added M6b (all 10 `GGGGCCCC` pairs incl. nested), M6c (weighted E_bp=−1),
  M6d (per-base pairing ≤ 1 invariant). Every expected value traces to the independent
  brute-force enumeration, not to library output. No assertions weakened, nothing skipped.

## Findings

- **D1 (Stage A, FIXED)** — false external-only base-pair-probability claim in Evidence
  §"mccaskill2" and algorithm-doc §2.2/§5.2/§5.3; corrected to the outside recursion with exact
  datasets.
- **D2 (Stage B, FIXED)** — external-only base-pair probabilities in code (register B5);
  corrected to the outside recursion; strict nested-pair + weighted + per-base-≤1 tests added.
- **End-state ✅ CLEAN** — code fixed, docs corrected, strict sourced tests added, full suite
  green (6596/0).
