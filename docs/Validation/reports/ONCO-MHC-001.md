# Validation Report: ONCO-MHC-001 — MHC-Peptide Binding Classification

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.ClassifyBindingAffinity(double)`, `ClassifyBindingRank(double, MhcClass)`, `IsValidPeptideLength(int, MhcClass)`, `ClassifyMhcBinding(int, double, MhcClass)` (`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:5306-5396`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Scope

This unit is the **classification half** of MHC-binding prediction: it maps a *caller-supplied* predicted binding affinity (IC50, nM) or %Rank percentile into Strong / Weak / NonBinder using community-standard cutoffs, and validates peptide length per MHC class. It does **not** predict the IC50/%Rank — that needs a trained model (NetMHCpan/MHCflurry) and is explicitly out of scope. There are therefore **no learned weights, PSSM, or anchor-residue computations** to validate; the only externally-sourced numbers are the threshold constants and the length ranges. The session validated each constant against the primary literature, not against the code.

## Stage A — Description

### Sources opened this session (independently re-retrieved)

1. **NetMHCpan-4.1 (Reynisson et al. 2020, *Nucleic Acids Res* 48(W1):W449–W454, PMC7319546)** — WebFetched full text. Verbatim:
   - Class I: "%Rank < 0.5% and %Rank < 2% thresholds are considered for detecting SBs and WBs for class I".
   - Class II: "%Rank < 2% and %Rank < 10%, for SBs and WBs for class II".
   - Class I length: "for class I, the length range goes from 8 to 14 amino acids, default is 8–11".
   - Class II length (NetMHCIIpan): "for class II only one length is admitted with 15 being the default value".
2. **IEDB threshold guidance** (help.iedb.org, articles 114094152371 / 114094151811) — via WebSearch result body. Verbatim:
   - "peptides with IC50 values <50 nM are considered high affinity, <500 nM intermediate affinity and <5000 nM low affinity".
   - "an absolute binding affinity (IC50) threshold of 500 nM identifies strong binders" (binder/non-binder demarcation).
3. **Sette et al. (1994, *J Immunol* 153(12):5586–92, PMID 7527444)** — confirmed via WebSearch snippet (PubMed page is reCAPTCHA-gated to WebFetch): "an affinity threshold of approximately 500 nM (preferably 50 nM or less) apparently determines the capacity" of a peptide to elicit a CTL response. This is the biological basis for the 50 nM (strong) / 500 nM (binder) cutoffs.
4. **IEDB MHC class II peptide length** — WebSearch confirmed: class II peptides "typically range between 13 and 25 amino acids long"; binding core "usually nine amino acids long (9-mer)".

### Formula / cutoff check (each constant traced to a source retrieved this session)

| Constant (code) | Value | Source | Match |
|---|---|---|---|
| `StrongBinderIc50Nm` | 50 | IEDB "<50 nM high affinity"; Sette "preferably 50 nM" | ✅ |
| `WeakBinderIc50Nm` | 500 | IEDB "<500 nM intermediate"; Sette "≈500 nM"; Roomp 2010 500 nM demarcation | ✅ |
| `ClassIStrongBinderRankPercent` | 0.5 | Reynisson "%Rank < 0.5% … SBs … class I" | ✅ |
| `ClassIWeakBinderRankPercent` | 2.0 | Reynisson "%Rank < 2% … WBs … class I" | ✅ |
| `ClassIIStrongBinderRankPercent` | 2.0 | Reynisson "%Rank < 2% … SBs … class II" | ✅ |
| `ClassIIWeakBinderRankPercent` | 10.0 | Reynisson "%Rank < 10% … WBs … class II" | ✅ |
| Class I length 8–11 | 8–11 | Reynisson default 8–11 (full 8–14) | ✅ (default; see Note 1) |
| Class II length 13–25 | 13–25 | IEDB class II tool description | ✅ |

The classification model — Strong iff `v < s`, Weak iff `s ≤ v < w`, NonBinder iff `v ≥ w`, all strict `<` — matches the verbatim "<" inequalities in every source. Boundary semantics (50 nM → Weak; 500 nM → NonBinder; %Rank 0.5 → Weak; 2.0 → NonBinder, class I) follow directly from the strict inequalities.

### Edge-case semantics
- IC50 must be a finite positive concentration → invalid (≤0, NaN, ±∞) throws. Sound (concentration semantics).
- %Rank is a percentile ∈ [0,100] → out-of-range / NaN throws. Sound (Reynisson defines %Rank as the top X% of random-peptide scores).
- Length outside the class range → not a valid candidate. Sound.
- INV-3 monotonicity and INV-4 strict-< are genuine properties of the nested-threshold model.

### Findings / divergences (Stage A notes)
- **Note 1 (assumption, documented):** Class I length range is the **8–11 default**, not the full **8–14** that Reynisson also admits. The code/doc disclose this explicitly and align it to the ONCO-NEO-001 windowing constants. Consequence: lengths 12–14 are reported invalid for class I. This is a documented default choice, not a sourcing error — hence **PASS-WITH-NOTES** rather than PASS.
- **Note 2:** The class II length range 13–25 is from IEDB, not from NetMHCpan-4.1 itself (NetMHCIIpan admits a configurable peptide length, default core handling 15). The Evidence already attributes 13–25 to IEDB correctly. No defect.

## Stage B — Implementation

### Code path reviewed
- `ClassifyBindingAffinity` (`OncologyAnalyzer.cs:5306`): validates finite & >0, then `< 50 → Strong`, else `< 500 → Weak : NonBinder`. Exact strict-`<` realisation. ✅
- `ClassifyBindingRank` (`:5337`): validates NaN / [0,100], selects per-class cutoffs, `< strong → Strong`, else `< weak → Weak : NonBinder`. ✅
- `IsValidPeptideLength` (`:5369`): inclusive `[min,max]` per class; total (non-positive lengths fall outside, return false). ✅
- `ClassifyMhcBinding` (`:5388`): length gate → NonBinder, else delegates to `ClassifyBindingAffinity` (so IC50 validation is propagated on valid length). ✅

### Cross-verification table recomputed vs code (values from sourced cutoffs, not from the code)

| Input | Sourced expected | Code result | Match |
|---|---|---|---|
| IC50 10 | Strong | Strong | ✅ |
| IC50 50 (boundary) | Weak | Weak | ✅ |
| IC50 200 | Weak | Weak | ✅ |
| IC50 500 (boundary) | NonBinder | NonBinder | ✅ |
| IC50 1000 | NonBinder | NonBinder | ✅ |
| %Rank 0.4 / I | Strong | Strong | ✅ |
| %Rank 0.5 / I (boundary) | Weak | Weak | ✅ |
| %Rank 2.0 / I (boundary) | NonBinder | NonBinder | ✅ |
| %Rank 1.5 / II | Strong | Strong | ✅ |
| %Rank 2.0 / II (boundary) | Weak | Weak | ✅ |
| %Rank 10.0 / II (boundary) | NonBinder | NonBinder | ✅ |
| length 8/11 / I | valid | valid | ✅ |
| length 12 / I | invalid | invalid | ✅ |
| length 13/25 / II | valid | valid | ✅ |
| length 26 / II | invalid | invalid | ✅ |

### Variant/delegate consistency
`ClassifyMhcBinding` is a thin wrapper that reuses `IsValidPeptideLength` + `ClassifyBindingAffinity`; results agree with the canonical methods by construction (verified by M20/M21/M21b).

### Test quality audit (HARD gate)
- **Sourced, not code-echoed:** every `Is.EqualTo` asserts the exact category derived from the published cutoffs. A deliberately-wrong implementation fails them: flipping any `<` to `<=` breaks M2/M4/M7/M9/M11b/M12; swapping class I/II cutoffs breaks M11/M11b/M13; widening a length bound breaks the length cases.
- **No green-washing:** all assertions are exact `Is.EqualTo` / `Is.True` / `Is.False` / `Throws<ArgumentOutOfRangeException>`; no Greater/AtLeast/ranges/tolerances, no skips, no commented-out tests.
- **Coverage gaps found and fixed (test-only, 0 code change):**
  1. **Class II strong-cutoff boundary** (`%Rank = 2.0`, class II → Weak) was unasserted — the class-II analogue of the class-I 0.5 boundary, and the only direct test of INV-4 on the class-II strong cutoff. Added `ClassifyBindingRank_ClassIITwoBoundary_ReturnsWeak` (sourced to Reynisson "%Rank < 2% … class II", 2.0<10).
  2. **`ClassifyMhcBinding` IC50-validation propagation** on a valid length was unasserted. Added `ClassifyMhcBinding_ValidLengthInvalidIc50_Throws` (len 9 valid, IC50 0 → throws per INV-01).
- **Remaining (non-defects):** S2 asserts `+∞` but not `−∞` (same `IsInfinity` branch); default `minLength`/`MhcClassI*` defaults are exercised transitively. These do not leave any sourced value untested.

### Findings / defects
None. All cutoffs and length ranges in the implementation match the externally-retrieved sources exactly. The two test gaps were closed in-session.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — every threshold/length constant traces to a source retrieved this session; the only divergence is the *documented* class I 8–11 default (vs full 8–14).
- **Stage B: PASS** — code realises the strict-`<` nested-threshold model exactly; cross-check table matches; tests now cover all four methods, both classes, all boundaries, and all validation branches.
- **Test-quality gate: PASS** — sourced exact assertions, no green-washing, full logic covered after adding two cases (25 → 27 tests for this unit; original count met).
- **End-state: CLEAN** — no code defect; two test-coverage gaps fixed; `dotnet build` 0 errors; full unfiltered suite **6663 passed, 0 failed**.
- No defect logged in FINDINGS_REGISTER beyond the test-coverage note.

## 2026-06-25 — Opt-in matrix-based predictor added (data-block partially lifted)

The earlier block ("no redistributable, cross-verifiable trained model") was revisited with web tools. An opt-in **matrix-based pMHC binding predictor** is now implemented, with the existing classification and all defaults unchanged.

### What was added
- `OncologyAnalyzer.PredictBindingHalfLifeBimas(peptide, matrix)` — BIMAS / Parker (1994) product rule: `T½ = FinalConstant · ∏ coefficients` (unlisted residue coefficient 1.0).
- `OncologyAnalyzer.PredictIc50Smm(peptide, matrix)` — SMM / Peters & Sette (2005): score = intercept + Σ contributions; `IC50 = 50000^(1 − score)` (algebraic inverse of the IEDB log50k linearisation).
- `OncologyAnalyzer.PredictAndClassifySmm(peptide, matrix)` — chains `PredictIc50Smm` → `ClassifyBindingAffinity` (predict→classify end-to-end).
- `OncologyAnalyzer.LoadScoringMatrix(lines)` + `PmhcScoringMatrix` / `PmhcScoringMethod` — caller-supplied matrix loader and types.

### Sources retrieved this session (verbatim scoring rules)
- **BIMAS** scoring documentation (NIH/CIT; restating Parker 1994), retrieved via Internet Archive (the live `www-bimas.cit.nih.gov` no longer resolves): "The initial (running) score is set to 1.0 … multiplied by the coefficient … multiplied by a final constant to yield an estimate of the half time of disassociation"; unlisted residue coefficient 1.00. The HLA-A2 coefficient table "is published … (except for HLA-A2)" in Parker 1994.
- **Parker, Bednarek & Coligan (1994)** PubMed abstract: "180 coefficients (20 amino acids x 9 positions)"; binding stability "calculated by multiplying together the corresponding coefficients" (≈ factor-of-5 accuracy).
- **SMM / IEDB log50k:** `log50k = 1 − log(IC50)/log(50000)` ⇒ `IC50 = 50000^(1 − score)` (Peters & Sette 2005, BMC Bioinformatics 6:132).

### Reproduced reference predictions (exact, independent of the code)
- SMM transform anchors: score 0 → **50000 nM**; score 1 → **1 nM**; score 0.5 → **√50000 = 223.6067977499790 nM** — asserted with `.Within` tolerances ≤ 1e-6.
- Strong-vs-non-binder ranking: an SMM matrix whose contributions for `GILGFVFTL` (influenza M1 58–66, the paradigm HLA-A*02:01 binder) sum to 1.0 → IC50 = 1 nM → **Strong**; a poly-W 9-mer → score 0 → IC50 = 50000 nM → **NonBinder**; the binder's IC50 (1 nM) is far below the non-binder's (50000 nM).
- BIMAS product: const 10 · (2.0·3.0·1.5) = **90.0**; `AAA` (unlisted) → 10 · 1·1·1 = **10.0**.

### Residual (honest)
No redistributable, cross-verifiable trained HLA coefficient matrix was obtainable this session: the public BIMAS coefficient files are served only by a now-defunct dynamic CGI (the Internet Archive captured the input form, not the generated value table), the Parker 1994 180-value table is paywalled, and the IEDB SMM matrices are non-commercial / no-redistribution. The library therefore bundles the published **scoring rules** (fully sourced and cross-verifiable) and a **caller-supplied matrix loader**, mirroring ONCO-IMMUNE-001 (CIBERSORT LM22). The pan-allele NetMHCpan/MHCflurry **neural** model remains out of scope.

### End-state
`dotnet build` 0 warnings/0 errors in changed files; this unit's fixture 43 passed / 0 failed (27 classification + 16 prediction); full unfiltered suite green (see commit). No code defect; classification and defaults unchanged. Validation Status confirmed **☐** (unchanged); Quick-Reference counts unchanged.
