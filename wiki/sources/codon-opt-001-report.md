---
type: source
title: "Validation report: CODON-OPT-001 (codon/sequence optimization — synonymous recoding to a host's preferred codons, CodonOptimizer.OptimizeSequence)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/CODON-OPT-001.md
sources:
  - docs/Validation/reports/CODON-OPT-001.md
source_commit: 9dfa56d90f2bf6d5d33c73c8ba12dd34ad267513
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CODON-OPT-001

The two-stage **validation write-up** for test unit **CODON-OPT-001** (sequence/codon optimization —
recoding a CDS by synonymous substitution toward a target organism's preferred codons, with CAI and
GC reported before/after), validated **2026-06-24**. This is the *report* artifact that feeds one row
of the [[validation-ledger]]; it records the validator's independent **verdict** on both the algorithm
description (Stage A) and the shipped code (Stage B). The wider campaign is [[validation-and-testing]];
[[test-unit-registry]] defines the unit. The algorithm, its five strategies, its invariants and edge
behaviour are synthesized in the concept [[codon-optimization]]. Distinct from
[[codon-opt-001-evidence]] — the pre-implementation evidence artifact sourced from `docs/Evidence/` —
this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** No code defect, no logged defect, no code change.
The codon-optimizer filter (`CodonOptimizer_OptimizeSequence_Tests` + `CodonOptimizer_CAI_Tests`) ran
**59 passed, 0 failed**; every asserted optimized sequence, change tuple, CAI and GC value reproduced.

## Canonical method & source under test

`CodonOptimizer.OptimizeSequence(string, CodonUsageTable, OptimizationStrategy, double gcTargetMin,
double gcTargetMax, double rareCodonThreshold)` in
`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs:239`. Code path reviewed:
`OptimizeSequence` (`:239–319`), `SelectOptimalCodon` (`:321–359`), `BalanceGcContent` (`:378–414`),
`CalculateCAI` (`:423–450`), `CalculateRelativeAdaptiveness` (`:452–468`), `SplitIntoCodons`
(`:687–695`), `TranslateCodon` (`:697–700`).

## Stage A — description (algorithm faithfulness)

Confirmed against **Wikipedia "Codon Adaptation Index"** (relative adaptiveness `w_i = f_i / max(f_j)`
over synonymous codons, `CAI = (∏ w_i)^(1/L)` = geometric mean; Sharp & Li 1987, *NAR* 15(3):1281–1295),
**Wikipedia "Codon usage bias"**, **Plotkin & Kudla 2011**, and the **Kazusa Codon Usage Database**
(species 316407 = E. coli W3110 K-12, fetched 2026-06-24). The `CalculateCAI` /
`CalculateRelativeAdaptiveness` code reproduces Sharp & Li exactly (`CAI = exp((1/L)·Σ ln w_i)`), and
canonical optimization (MaximizeCAI) = the max-frequency synonymous codon per residue. The defining
invariant is **`translate(optimize(seq)) == translate(seq)`** — synonymous substitution never alters
the encoded protein.

**Independent numeric cross-check.** The built-in E. coli table was **re-derived by hand from Kazusa
316407** (per-AA fraction = codon per-thousand ÷ Σ synonymous per-thousand): every tested value
(Leu CUG 0.50 / CUA 0.04, Arg CGC 0.40 / CGU 0.38, Ala GCG 0.36, Pro CCG 0.53, Thr ACC 0.44, …) matches
the stored 2-decimal precision. Hand-computed CAIs matched the test assertions exactly: `CUAAGACGA`
(E. coli) = **0.1063**; `CUGCCGACC` (Human) = **0.7005**; `AUGGCUUAA` original (E. coli) = **0.6667**.

**Edge-case semantics (all sourced):** empty → empty, CAI 0; single-codon AAs Met (AUG) / Trp (UGG)
unchanged; stop codons (UAA/UAG/UGA) preserved verbatim; length not a multiple of 3 trimmed to complete
codons; DNA `T → U` conversion; case-insensitive (uppercased).

**Sole documented divergence (non-blocking):** the CAI clamp to `1e-6` for codons *absent from the
table* diverges from Sharp & Li's `0.5/N` convention, but it fires only on **incomplete custom tables**
— all three built-in organism tables are complete, so it never triggers for the canonical organisms.

## Stage B — implementation

**Protein-preservation invariant holds by construction.** `SelectOptimalCodon` only ever returns a
member of `AminoAcidToCodons[aminoAcid]` (a synonym of the same AA) or the current codon; stop codons
(`aminoAcid == "*"`) bypass selection and are copied verbatim (`:273–278`); the BalancedOptimization GC
phase (`BalanceGcContent`) also substitutes only within the synonym set and skips stops; unknown codons
(no synonym set, translate to "X") are returned unchanged. Independently locked by
`OptimizeSequence_Invariant_ProteinPreserved` (re-translates original vs optimized via a *separate*
genetic-code map, 5 inputs) and `OptimizeSequence_PreservesProtein_AllStrategies` (all five strategies).

**Max-weight synonym rule.** `MaximizeCAI` = `synonymousCodons.OrderByDescending(freq).First()`.
`CUAAGACGA` → `CUGCGCCGC` (CUA→CUG, AGA→CGC, CGA→CGC), CAI 1.0. Worked examples recomputed vs code:
`AUGGCUUAA` BalancedOpt → `AUGGCGUAA`, 1 change (index 3, GCU→GCG), OrigCAI 0.667 / OptCAI 1.0, GC
3/9→4/9; E. coli vs Yeast on `CUGAGA` → `CUGCGC` vs `UUGAGA` (Leu best CUG 0.50 / UUG 0.29; Arg best
CGC 0.40 / AGA 0.48); Human vs E. coli CAI of `CUGCCGACC` = 1.0 vs 0.700.

**Variant/delegate consistency.** All 5 strategies preserve protein. `HarmonizeExpression` uses
weighted-random synonym selection (`new Random()` per call → **intentionally non-deterministic**), so
its test asserts only protein + CAI range. `AvoidRareCodons` replaces only sub-threshold codons.
`MinimizeSecondary` falls through to BalancedOptimization in codon *selection*; the dedicated
`ReduceSecondaryStructure` method is separate. **Numerical robustness:** geometric mean in log space
(no underflow); `maxFreq <= 0` → NaN guard skips no-data AAs; empty/zero-codon → 0; `SplitIntoCodons`
(`i+2 < len`) keeps all complete codons and drops only the incomplete tail.

**Test-quality audit.** 59 assertions check exact sourced values (specific optimized sequences per
Kazusa, hand-computed CAI/GC, exact change tuples), independently re-translate to lock protein
preservation, and cover every Stage-A edge case — real assertions, not no-throw tautologies. The single
stochastic strategy asserts only order-independent invariants.

## Findings

- **No code defect, no test change, no code change (State CLEAN).** Every worked example and cross-check
  reproduced exactly.
- **Non-defect notes carried (documented, none affects the canonical model or any invariant):**
  (1) `HarmonizeExpression` is intentionally stochastic; (2) `MinimizeSecondary` selects codons via
  BalancedOptimization while secondary-structure reduction lives in the separate `ReduceSecondaryStructure`;
  (3) the CAI `1e-6` zero-frequency clamp applies only to incomplete custom tables, never to the built-in
  E. coli / yeast / human organisms.
