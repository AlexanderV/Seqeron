# Test Specification: MHC-MATRIX-001

**Test Unit ID:** MHC-MATRIX-001
**Area:** Oncology
**Algorithm:** SMM / BIMAS Matrix pMHC Prediction
**Status:** ☑ Validated (Stage A ✅ / Stage B ✅ / CLEAN) — 2026-06-25
**Last Updated:** 2026-06-25

> Independently re-validated under the two-stage protocol (see
> `docs/Validation/reports/MHC-MATRIX-001.md`). The matrix is **caller-supplied** (no redistributable trained
> HLA matrix: BIMAS CGI defunct, Parker 1994 table paywalled, IEDB SMM non-commercial); validation targets the
> published scoring rules + the predict→classify chain, an acceptable documented boundary (ONCO-MHC-001).

---

## 1. Evidence Summary

| # | Source | Confirms |
|---|--------|----------|
| 1 | Peters & Sette (2005), *BMC Bioinformatics* 6:132 (PMC1173087) | SMM = additive sum of position-specific contributions **plus a constant offset/intercept** (`H·w`, first column = 1). |
| 2 | IEDB MHC-I SMM linearisation | `log50k = 1 − log(IC50)/log(50000)` ⇒ `IC50 = 50000^(1 − score)`; affinity > 50000 nM ⇒ score 0. |
| 3 | Parker, Bednarek & Coligan (1994), *J. Immunol.* 152(1):163–175 (PMID 8254189) + BIMAS docs | BIMAS T½ = final-constant × **product** of per-position coefficients (180 = 20 aa × 9 pos); unlisted residue = 1.00; running score init 1.0. |
| 4 | IEDB / Sette 1994 thresholds | IC50 < 50 nM strong, < 500 nM weak, else non-binder (strict). |

## 2. Canonical Method(s)

`PredictIc50Smm`, `PredictBindingHalfLifeBimas`, `PredictAndClassifySmm`, `LoadScoringMatrix`
(+ `PmhcScoringMatrix` record, `SmmIc50Base` const)

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMhcBinding_Tests.cs`

## 3. Contract / Invariants

- R: SMM IC50 = `50000^(1 − (intercept + Σ contribution_i))` > 0; missing residue ⇒ +0 (additive identity).
- R: BIMAS T½ = `FinalConstant · ∏ coeff_i` ≥ 0; missing residue ⇒ ×1.0 (multiplicative identity).
- R: peptide length ≠ matrix rows / empty matrix ⇒ ArgumentException; null ⇒ ArgumentNullException.
- M: anchor (favorable residue) match ⇒ stronger binding (higher BIMAS T½ / lower SMM IC50).
- Chain: `PredictAndClassifySmm` = `ClassifyBindingAffinity(PredictIc50Smm(...))` (strong<50/weak<500).

## 4. Cross-check / Differential Oracle

Hand-computed against the published transforms (Python, this session) — all reproduced by the code:
SMM score 0→50000, 0.5→√50000=223.60679774997897, 1→1; multi-pos "KAY" score 0.60→75.78582832551992 (Weak);
"GILGFVFTL" Σ=1.0→1 nM (Strong); BIMAS "LMV"×10→90, "AAA"×10→10, "LV"→20, "AA"→0.02, loader "LM"×10→60,
no-CONST "LM"→6.0.

## 5. Validation Checklist (✅ all met)

- [x] Stage A: every source retrieved this session; SMM additive-sum-plus-offset + IC50 transform + BIMAS product confirmed against the publications.
- [x] Stage B: implementation traced; cross-check table recomputed vs code, all match.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Genomics 18756 passed), 0 warnings on changed file.
- [x] 4 tests added (multi-position additive sum + missing residue; non-AA; Weak band of chain; loader default constant).
- [x] Flipped `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the `docs/checklists/*.md`.
