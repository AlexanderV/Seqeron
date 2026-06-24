# Validation Report: PHYLO-DIST-001 — Phylogenetic Distance Matrix (p-distance, Jukes–Cantor JC69, Kimura 2-parameter K2P)

- **Validated:** 2026-06-24   **Area:** Phylogenetic Analysis
- **Canonical method(s):** `PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs, method)`, `PhylogeneticAnalyzer.CalculatePairwiseDistance(s1, s2, method)`
  - Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs`
  - Tests: `tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_DistanceMatrix_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

This is an independent re-validation in a fresh context. The source method is unchanged since the
prior validation (last distance-relevant commit `61a4c923`); formulas re-derived by hand and
re-confirmed against Wikipedia "Models of DNA evolution" (fetched 2026-06-24).

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "Models of DNA evolution"** (fetched 2026-06-24) confirms verbatim:
  - JC69: `d̂ = −(3/4) ln(1 − (4/3)p)`, **natural log**, constants 3/4 and 4/3; `p` = proportion of differing sites.
  - K80/K2P: `K = −(1/2) ln((1 − 2p − q)·√(1 − 2q))`, **natural log**, constant 1/2; `p` = proportion of transitional differences, `q` = proportion of transversional differences.
- Primary refs cited in the spec/evidence: Jukes & Cantor (1969); Kimura M (1980) *J. Mol. Evol.* 16:111–120.
- Transition classification (purine↔purine A↔G, pyrimidine↔pyrimidine C↔T; all other changes are transversions) is the standard genetics convention.

### Formula check
- **p-distance** = (differing sites)/(compared sites). ✓
- **JC69** `d = −(3/4)·ln(1 − (4/3)·p)`, domain `1 − 4p/3 > 0` ⇔ `p < 3/4` (else undefined → +∞). ✓
- **K2P** `K = −(1/2)·ln((1 − 2P − Q)·√(1 − 2Q))`. ✓
  - The prompt's alternate form `½·ln(1/(1−2P−Q)) + ¼·ln(1/(1−2Q))` is algebraically identical:
    `−½ln((1−2P−Q)√(1−2Q)) = −½ln(1−2P−Q) − ¼ln(1−2Q) = ½ln(1/(1−2P−Q)) + ¼ln(1/(1−2Q))`. ✓
  - Domain guards: `1 − 2P − Q > 0` and `1 − 2Q > 0`. ✓
- **Log base:** natural log (ln) throughout. ✓

### Edge-case semantics
- Gap/missing handling = **pairwise deletion**: a column is dropped if either sequence has `-` or a non-ACGT (ambiguous IUPAC) base. ✓ (MEGA-standard.)
- Identical sequences → 0; all-gap / no comparable sites → 0 (limit). ✓
- JC saturation at `p ≥ 3/4` → +∞; K2P saturation when `1−2P−Q ≤ 0` or `1−2Q ≤ 0` → +∞. ✓

### Independent cross-check (hand computed this session, natural log)
| Case | Inputs | Formula | Hand value | Spec/Test |
|------|--------|---------|-----------|-----------|
| JC, p=1/8 | 1 diff / 8 | −0.75·ln(0.833333) | **0.1367412** | 0.13674 ✓ |
| JC, p=1/4 | 2 diffs / 8 | −0.75·ln(0.666667) | **0.3040988** | ≈0.304 ✓ |
| K2P transition | S=1/4, V=0 | −0.5·ln(0.5·√1) | **0.3465736** | 0.34657 ✓ |
| K2P transversion | S=0, V=1/4 | −0.5·ln(0.75·√0.5) | **0.3171607** | 0.31713 ✓ |
| K2P mixed | S=1/8, V=1/8 | −0.5·ln(0.625·√0.75) | **0.3070114** | matches formula ✓ |

### Invariants
Symmetry, zero diagonal, non-negativity, JC≥p, K2P≥p, and the two saturation conditions are genuine mathematical properties of the validated formulas. ✓

**Stage A: PASS** — every formula, constant, log base, and the transition/transversion classification match Wikipedia "Models of DNA evolution" and the cited primary literature; hand computations confirm spec/test values to required precision.

## Stage B — Implementation

### Code path reviewed
`CalculatePairwiseDistance` (`PhylogeneticAnalyzer.cs:223–270`), helpers `IsStandardBase` (272–273), `IsTransition` (275–281), `JukesCantorDistance` (283–289), `Kimura2ParameterDistance` (291–298); matrix wrapper `CalculateDistanceMatrix` (199–218).

### Formula realised correctly?
- p-distance: `p = differences / comparableSites` (line 258). ✓
- JC69 (286–288): `arg = 1 − 4p/3; if (arg <= 0) +∞; else −0.75·Math.Log(arg)`. `Math.Log` is natural log. Constants and domain guard correct. ✓
- K2P (294–297): `arg1 = 1 − 2s − v; arg2 = 1 − 2v; if (arg1<=0 || arg2<=0) +∞; else −0.5·Math.Log(arg1·Math.Sqrt(arg2))`. √ applied to (1−2v) term only; both domain guards present. ✓
- Transition classification (278–280): both purines (A/G) or both pyrimidines (C/T) → transition; else transversion. Not swapped. ✓
- Gap/ambiguous handling (242–243): skips `-` and any non-ACGT char (pairwise deletion); `comparableSites==0 → return 0` (256). ✓
- Case-insensitive via `char.ToUpperInvariant` (238–239). ✓
- Pre-conditions: null → `ArgumentNullException` (226–227); unequal length → `ArgumentException` (228–229). ✓
- `CalculateDistanceMatrix` (199–218) delegates to `CalculatePairwiseDistance` for every i<j and mirrors to [j,i] → symmetric, zero-diagonal matrix; no divergent reimplementation. ✓

### Cross-verification vs code
Ran `--filter FullyQualifiedName~DistanceMatrix` → **32 passed / 0 failed**. M08 (JC=0.13674), M09 (K2P transition 0.34657 / transversion 0.31713), S06 (mixed K2P matches formula) all pass and match the Stage A hand computations above.

### Test quality audit
26 spec-mapped tests + extras (BothGaps, JC_TwoDifferences p=0.25, pairwise symmetry). Assertions check exact sourced numeric values (hardcoded 0.13674 / 0.34657 / 0.31713) plus formula recomputation, not "no throw" tautologies. Edge cases covered: identical→0 (M01), JC p≥0.75→+∞ (M13), K2P V≥0.5→+∞ (S07), all-gap→0 (S04), ambiguous bases skipped (M14), null/unequal-length throws (M10/M15). Deterministic.

### Findings / defects
None. The implementation faithfully realises the validated formulas with correct constants, natural log, correct transition/transversion classification, and both saturation guards.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.**
- No code changes. `DistanceMatrix` filter → 32 passed / 0 failed. Build succeeded, 0 warnings.
