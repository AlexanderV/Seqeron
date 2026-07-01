# Validation Report: MHC-MATRIX-001 — SMM / BIMAS Matrix pMHC Prediction

- **Validated:** 2026-06-25   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.PredictIc50Smm`, `PredictBindingHalfLifeBimas`, `PredictAndClassifySmm`, `LoadScoringMatrix` (+ `PmhcScoringMatrix` record, `SmmIc50Base`)
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Scope note

No redistributable trained HLA coefficient matrix ships with the library (BIMAS CGI defunct; Parker 1994
table paywalled; IEDB SMM matrices non-commercial/no-redistribution). The matrix is therefore
**caller-supplied** and only the published *scoring rules* + a loader are implemented. This is an
**acceptable, documented boundary** (ONCO-MHC-001), not a limitation: the validation targets the
SMM/BIMAS math and the predict→classify chain against the published transforms, not a bundled matrix.

## Stage A — Description

### Sources opened (this session)
- **Peters & Sette (2005), *BMC Bioinformatics* 6:132** — "Generating quantitative models … stabilized
  matrix method (SMM)". Via PMC (PMC1173087): the SMM prediction is the linear model **H·w = y_pred**, where
  the scoring matrix `w` "quantifies the contribution of each residue at each position" and **"the first
  column of each row is set to 1, which serves as a constant offset added to each prediction"** — i.e. the
  score is the *additive sum of position-specific contributions plus a constant intercept*. Confirms the
  code's `score = intercept + Σ_i contribution_i`.
- **IEDB MHC-I tool description / SMM linearisation** — IC50 (nM) is log-transformed via
  **`log50k = 1 − log(IC50)/log(50000)`**; peptides with affinity > 50000 nM map to log-transformed value 0.
  Inverting: **`IC50 = 50000^(1 − score)`**. Confirms `PredictIc50Smm` and `SmmIc50Base = 50000`.
- **Parker, Bednarek & Coligan (1994), *J. Immunol.* 152(1):163–175 (PMID 8254189) + BIMAS docs** — table of
  **180 coefficients (20 aa × 9 positions)**; the **"theoretical binding stability was calculated by
  multiplying together the corresponding coefficients"**. BIMAS scoring: running score starts at 1.0, is
  multiplied by each position's coefficient (unlisted/ambiguous residue = 1.00, leaving the score unchanged),
  then multiplied by a final constant to give the predicted half-time of dissociation. Confirms
  `PredictBindingHalfLifeBimas = FinalConstant · ∏_i coeff_i`, neutral default 1.0.
  (The live BIMAS CGI is unreachable — `ECONNREFUSED` — consistent with the code's "defunct" claim.)
- **IEDB / Sette 1994 thresholds** — IC50 < 50 nM high (strong), < 500 nM intermediate (weak), else
  non-binder; strict inequalities. Confirms `ClassifyBindingAffinity` cutoffs used by the chain.

### Formula check
- **SMM IC50:** `IC50 = 50000^(1 − (intercept + Σ contribution_i))`, missing residue → +0. ✅ matches source.
- **BIMAS T½:** `FinalConstant · ∏ coeff_i`, missing residue → ×1.0, running score init 1.0. ✅ matches source.
- **Chain:** `PredictAndClassifySmm` = `ClassifyBindingAffinity(PredictIc50Smm(...))`, strong<50/weak<500. ✅.

### Edge-case semantics (sourced/defined)
- Peptide length ≠ matrix rows → ArgumentException (defined contract). ✅
- Empty matrix (no rows) → ArgumentException. ✅
- Null peptide / null lines → ArgumentNullException. ✅
- Missing residue → neutral (SMM +0 additive identity / BIMAS ×1.0 multiplicative identity) — *sourced*. ✅
- Non-AA character → simply unlisted ⇒ neutral default (consistent with the "ambiguous residue" rule). ✅
- Malformed / non-numeric / multi-char token in loader → FormatException (defined). ✅

### Independent cross-check (hand-computed, this session — NOT code echoes)
| Method | Input | Transform | Expected |
|---|---|---|---|
| SMM | score 0 | 50000^(1−0) | **50000.0** |
| SMM | score 1 | 50000^(1−1) | **1.0** |
| SMM | score 0.5 | √50000 | **223.60679774997897** |
| SMM | intercept 0.3 + contrib 0.2 = 0.5 | √50000 | **223.60679774997897** |
| SMM (multi-pos) | "KAY", intercept 0.05, K=0.30, V=0.20(unlisted A⇒0), Y=0.25 ⇒ score 0.60 | 50000^0.40 | **75.78582832551992** ⇒ Weak |
| SMM chain | "LV", intercept 0.1, L=0.3, V=0.2 ⇒ score 0.6 | 50000^0.40 | **75.79 nM ⇒ Weak** |
| SMM chain | "GILGFVFTL", Σ contrib = 1.0 | 50000^0 | **1.0 nM ⇒ Strong** |
| SMM chain | "WWWWWWWWW" (no match) | 50000^1 | **50000 nM ⇒ NonBinder** |
| BIMAS | "LMV", const 10, 2·3·1.5 | product×const | **90.0** |
| BIMAS | "AAA" (unlisted) | 10·1·1·1 | **10.0** |
| BIMAS | "LV" 5·4, const 1 | product×const | **20.0** |
| BIMAS | "AA" 0.1·0.2, const 1 | product×const | **0.02** |
| BIMAS | "LM" 2·3, const 10 (loader) | product×const | **60.0** |
| BIMAS | "LM" 2·3, no CONST ⇒ const 1.0 | product×const | **6.0** |

All values reproduced by an independent Python computation of the *published transforms*.

**Stage A verdict: ✅ PASS** — all formulas, constants and edge conventions match the cited primary sources.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `PmhcScoringMatrix` record (8341), `SmmIc50Base = 50000` (8350), neutral defaults (8357/8363).
- `LoadScoringMatrix` (8382): parses `RESIDUE=VALUE` / `CONST=VALUE`, comments, blanks, uppercases; default
  FinalConstant = 1.0 identity.
- `PredictBindingHalfLifeBimas` (8470): `score=1.0`; `score *= coeff` (missing → 1.0); `return score·FinalConstant`. ✅
- `PredictIc50Smm` (8499): `score=FinalConstant` (intercept); `score += value` (missing → 0); `return 50000^(1−score)`. ✅
- `PredictAndClassifySmm` (8526): chains `ClassifyBindingAffinity(PredictIc50Smm(...))`. ✅
- `ValidateMatrixAgainstPeptide` (8533): empty-matrix + length-mismatch guards. ✅

### Formula realised correctly?
Yes — the additive-sum-with-intercept (SMM) and multiplicative-product-with-final-constant (BIMAS) transforms
are realised exactly as in the sources; the IC50 transform uses `Math.Pow(50000, 1−score)`. Verified by tracing
and by the cross-check table above (recomputed against the actual code via the test fixture, all match).

### Variant / delegate consistency
`PredictAndClassifySmm` calls the same `PredictIc50Smm` and the established `ClassifyBindingAffinity`; no divergent
duplicate path. Loader round-trips into both predictors consistently (P14/P14b).

### Numerical robustness
`Math.Pow` on finite scores is stable across the validated range (score 0→50000, score 1→1). Product/sum do not
overflow on the tested coefficient ranges. Missing-residue defaults avoid any KeyNotFound/NaN.

### Test-quality audit (HARD gate)
Existing 43-test fixture covers every public method/overload with exact, source-traced values (not code echoes).
Audit found three genuine Stage-A-path gaps; **added 4 tests** locking hand-computed values:
- `PredictIc50Smm_MultiPositionAdditiveSumWithMissingResidue_HandComputed` — multi-position additive sum +
  intercept + missing residue together (prior IC50 tests were single-position only). IC50 = 75.78582832551992.
- `PredictIc50Smm_NonAminoAcidCharacter_ContributesZero` — non-AA char → neutral (Stage-A "non-AA" edge).
- `PredictAndClassifySmm_WeakBand_ReturnsWeak` — the middle (Weak) classification band of the chain, untested
  by P9 (Strong) / P10 (NonBinder). IC50 = 75.79 nM ⇒ Weak.
- `LoadScoringMatrix_NoConst_DefaultsToIdentityOne` — default FinalConstant = 1.0 identity (T½ = 6.0).
All expected values trace to the published transform / hand-computation. Full unfiltered suite green.

**Stage B verdict: ✅ PASS** — implementation faithfully realises the validated description; suite is real,
exact, and now covers the additive multi-position path, the Weak band, non-AA input and the loader default.

## Verdict & follow-ups
- **State: ✅ CLEAN.** No defect found. The caller-supplied matrix is a documented, acceptable boundary; the
  SMM/BIMAS math and the predict→classify chain verify against the published transforms (Peters & Sette 2005;
  Parker 1994 / BIMAS; IEDB log50k). 4 tests added to harden coverage. No findings logged.
- Full suite: `dotnet test Seqeron.sln -c Debug` → Failed: 0, Passed: 18756 (Genomics), 0 warnings on changed file.
