# Validation Report: PHYLO-DIST-001 — Phylogenetic Distance Matrix (p-distance, Jukes–Cantor, Kimura 2-parameter)

- **Validated:** 2026-06-12   **Area:** Phylogenetic Analysis
- **Canonical method(s):** `PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs, method)`, `PhylogeneticAnalyzer.CalculatePairwiseDistance(s1, s2, method)`
  - Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs`
  - Tests: `tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_DistanceMatrix_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "Models of DNA evolution"** (fetched 2026-06-12). Confirms verbatim:
  - JC69: `d̂ = −(3/4) ln(1 − (4/3)p)`, natural log, constants 3/4 and 4/3.
  - K80/K2P: `K = −(1/2) ln((1 − 2p − q)·√(1 − 2q))`, natural log, constant 1/2, with `p` = proportion of transitional differences and `q` = proportion of transversional differences.
- **Wikipedia, "Transition (genetics)"** (fetched 2026-06-12). Confirms transition classification:
  - Transitions: A↔G (purine↔purine), C↔T (pyrimidine↔pyrimidine).
  - Transversions: A↔C, A↔T, G↔C, G↔T (purine↔pyrimidine).
- Primary refs cited in the spec/evidence: Jukes & Cantor (1969); Kimura M (1980) *J. Mol. Evol.* 16:111–120.

### Formula check
- **p-distance** = (differing sites)/(compared sites). ✓ matches spec.
- **Jukes–Cantor** `d = −(3/4)·ln(1 − (4/3)·p)`, domain `p < 3/4` (else `1 − 4p/3 ≤ 0`, undefined → +∞). ✓
- **Kimura 2-parameter** `K = −(1/2)·ln((1 − 2P − Q)·√(1 − 2Q))`. ✓
  - The prompt's alternate form `½·ln(1/(1−2P−Q)) + ¼·ln(1/(1−2Q))` is algebraically identical:
    `−½ln((1−2P−Q)√(1−2Q)) = −½ln(1−2P−Q) − ¼ln(1−2Q) = ½ln(1/(1−2P−Q)) + ¼ln(1/(1−2Q))`. ✓
  - Domain guards: `1 − 2P − Q > 0` and `1 − 2Q > 0`. ✓
- **Log base:** natural log (ln) for all corrected distances. ✓ (no log10 mistake).

### Edge-case semantics
- Gap/missing-data handling = **pairwise deletion**: a column is dropped if either sequence has `-` or a non-ACGT (ambiguous IUPAC) base. Sourced (Distance matrices in phylogeny; pairwise deletion is standard, e.g. MEGA). ✓
- Identical sequences → 0; all-gap / no comparable sites → 0 (0 differences / limit). ✓
- JC saturation at p ≥ 3/4 → +∞; K2P saturation when `1−2P−Q ≤ 0` or `1−2Q ≤ 0` → +∞. ✓

### Independent cross-check (hand computed, natural log)
| Case | Inputs | Formula | Hand value | Spec/Evidence |
|------|--------|---------|-----------|---------------|
| JC, p=1/8 | ACGTACGT vs TCGTACGT | −0.75·ln(1−4·0.125/3)=−0.75·ln(0.833333) | **0.136741** | 0.13674 ✓ |
| JC, p=1/4 | 2 diffs / 8 | −0.75·ln(0.666667) | **0.304099** | ≈0.304 ✓ |
| K2P transition | S=1/4, V=0 | −0.5·ln(0.5·√1) | **0.346574** | 0.34657 ✓ |
| K2P transversion | S=0, V=1/4 | −0.5·ln(0.75·√0.5) | **0.317162** | 0.31713 ✓ |
| K2P mixed | S=1/8, V=1/8 | −0.5·ln(0.625·√0.75) | **0.306955** | ≈0.30696 ✓ |

(Note: evidence §5.3 lists pure-transversion as ≈0.31726; the precise value is 0.317162. The test M09 asserts 0.31713 within 1e-4, which 0.317162 satisfies. Minor cosmetic rounding in the evidence prose only; no code/test impact.)

### Invariants
Symmetry, zero diagonal, non-negativity, JC≥p, K2P≥p, and the two saturation conditions are all genuine mathematical properties of the validated formulas. ✓

**Stage A: PASS** — every formula, constant, log base, and the transition/transversion classification match authoritative sources; hand computations confirm spec values to required precision.

## Stage B — Implementation

### Code path reviewed
`CalculatePairwiseDistance` (`PhylogeneticAnalyzer.cs:165–212`), helpers `IsStandardBase` (214), `IsTransition` (217–223), `JukesCantorDistance` (225–231), `Kimura2ParameterDistance` (233–240).

### Formula realised correctly?
- p-distance: `p = differences / comparableSites` (line 200). ✓
- JC69 (line 228–230): `arg = 1 − 4p/3; if (arg <= 0) +∞; else −0.75·Math.Log(arg)`. `Math.Log` is natural log. ✓ Constants 4/3 and 3/4 correct, domain guard correct.
- K2P (line 236–239): `arg1 = 1 − 2s − v; arg2 = 1 − 2v; if (arg1<=0 || arg2<=0) +∞; else −0.5·Math.Log(arg1·Math.Sqrt(arg2))`. ✓ Constant 1/2 correct, √ on the (1−2v) term only, both domain guards present.
- Transition classification (line 220–222): both purines (A/G) or both pyrimidines (C/T) → transition; else transversion. Matches Wikipedia. ✓ Not swapped.
- Gap/ambiguous handling (line 184–185): skips `-` and any non-ACGT char (pairwise deletion). `comparableSites==0 → return 0` (line 198). ✓
- Case-insensitive via `ToUpperInvariant` (line 180–181). ✓
- Pre-conditions: null → `ArgumentNullException` (168–169); unequal length → `ArgumentException` (170–171). ✓

### Cross-verification recomputed vs code (via tests)
M08 asserts JC = 0.13674; M09 asserts K2P transition 0.34657 / transversion 0.31713; S06 asserts mixed K2P matches formula — all pass, matching the Stage A hand computations.

### Variant/delegate consistency
`CalculateDistanceMatrix` (141–158) delegates to `CalculatePairwiseDistance` for every i<j and mirrors to [j,i] giving a symmetric, zero-diagonal matrix; no divergent reimplementation. ✓

### Test quality audit
26 spec-mapped tests + 6 extra; assertions check exact sourced numeric values (hardcoded references 0.13674 / 0.34657 / 0.31713) plus formula recomputation, not just "no throw". Edge cases covered: identical→0 (M01), JC p≥0.75→+∞ (M13), K2P V≥0.5→+∞ (S07), all-gap→0 (S04), ambiguous bases skipped (M14), null/unequal-length throws (M10/M15). Deterministic.

### Findings / defects
None. Implementation faithfully realises the validated formulas with correct constants, natural log, correct transition/transversion classification, and both saturation guards.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.**
- Tests: `--filter FullyQualifiedName~DistanceMatrix` → 32 passed / 0 failed. Full suite `Seqeron.Genomics.Tests` → 4461 passed / 0 failed.
- No code changes required. One cosmetic note: evidence §5.3 rounds the pure-transversion K2P to 0.31726 in prose; exact value is 0.317162 (test asserts 0.31713 ±1e-4 — fine). Not a defect.
