# Validation Report: ONCO-MHC-001 — MHC-Peptide Binding Classification

- **Validated:** 2026-06-26 (fresh re-validation after the 8–11 → 8–14 length-window fix, commit 66c24491; supersedes 2026-06-25)   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.MhcClassIMaxPeptideLength` (now 14), `IsValidPeptideLength(int, MhcClass)`, `GenerateNeoantigenPeptides(...)` default window, `ClassifyBindingAffinity(double)`, `ClassifyBindingRank(double, MhcClass)`, `ClassifyMhcBinding(int, double, MhcClass)`; and the MHCflurry pan-allele predictor `MhcflurryAffinityPredictor.EncodePeptide` / `PredictIc50` / `PredictAndClassify` (`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:7954,8268,8025`; `MhcflurryAffinityPredictor.cs:250,408,640`).
- **Stage A verdict:** PASS (length window 8–14 = full NetMHCpan-4.1 class I window, now matching source exactly — the prior 8–11 PASS-WITH-NOTES divergence is resolved)
- **Stage B verdict:** PASS

---

## 2026-06-26 — Fresh re-validation after the class I length-window widening (8–11 → 8–14)

Commit **66c24491** widened the class I peptide-length window from 8–11 to **8–14** (`MhcClassIMaxPeptideLength` 11 → 14), propagating to `IsValidPeptideLength` (class I gate) and the `GenerateNeoantigenPeptides` default window. ONCO-MHC-001 was reset to ⬜ pending; this is its fresh re-validation, focusing on the length change and confirming the rest still holds. No prior repo TestSpec/Evidence/test was trusted; sources and the oracle were retrieved THIS session.

### Stage A — sources retrieved this session
- **NetMHCpan-4.1 (Reynisson et al. 2020, *Nucleic Acids Res.* 48(W1):W449–W454, PMC7319546)** — WebFetched the primary article. Verbatim: *"for class I, the length range goes from 8 to 14 amino acids, default is 8–11"*; *"a %Rank < 0.5% and %Rank < 2% thresholds are considered for detecting SBs and WBs for class I and %Rank < 2% and %Rank < 10%, for SBs and WBs for class II"*. The widened **8–14** window is therefore the FULL class I range NetMHCpan-4.1 accepts (the previous 8–11 was the *default* sub-range). The %Rank thresholds (0.5/2 class I, 2/10 class II) and the IC50 cutoffs (50/500 nM, IEDB/Sette, confirmed in prior sessions) are unchanged by this fix.
- **MHCflurry encoder (`encodable_sequences.py`)** — read the installed `mhcflurry` 2.1.5 source (`Library/Python/3.9/.../mhcflurry/encodable_sequences.py`). The `left_pad_centered_right_pad` branch: `min_length = 5` ("We arbitrarily set a minimum length of 5"), `max_length = 15` (default), encoded width `3 * max_length`; left edge `result[:length]`, right edge `result[-length:]`, center `center_left_padding = floor((max_length - length)/2)`, `center_left_offset = max_length + center_left_padding`. This is reproduced verbatim by the C# `EncodePeptide` (`MhcflurryAffinityPredictor.cs:277–294`), with `PeptideMinLength=5`, `PeptideMaxLength=15`. The 8–14 class I window sits entirely inside the encoder's supported 5–15, so lengths 12/13/14 encode and score with no special-casing.

**Stage A: PASS.** The length window now equals the source's stated class I range exactly — the only prior divergence (the documented 8–11 default vs full 8–14) is *resolved* by this fix, so the note is retired and Stage A is a clean PASS.

### Stage B — independent oracle cross-check (lengths 12/13/14)

Re-ran `mhcflurry` 2.1.5 in-session, loading the **exact embedded fixture network** `PAN-CLASS1-1-3ed9fb2d2dcc9803` (the smallest `models_class1_pan` member, release 20200610, present locally as the full 10-network pack) via `Class1NeuralNetwork.from_config` + its `weights_*.npz`, and predicting through the official `predict()` path with the `Class1AffinityPredictor` allele-pseudosequence table — an oracle wholly independent of the C# code. Comparison to the test goldens:

| Peptide / allele | len | Test golden (nM) | This-session oracle (nM) | Rel. diff |
|---|---|---|---|---|
| SIINFEKL / HLA-A\*02:01 | 8 | 11483.195201 | 11483.199313 | 3.6e-7 |
| GILGFVFTL / HLA-A\*02:01 | 9 | 19.123150 | 19.123144 | 3.1e-7 |
| NLVPMVATV / HLA-A\*02:01 | 9 | 17.542640 | 17.542638 | 1.1e-7 |
| ELAGIGILTV / HLA-A\*02:01 | 10 | 119.054961 | 119.054974 | 1.1e-7 |
| AAAWYLWEV / HLA-A\*02:01 | 9 | 16.559303 | 16.559300 | 1.8e-7 |
| SIINFEKL / HLA-B\*07:02 | 8 | 28830.796646 | 28830.799575 | 1.0e-7 |
| SLYNTVATL / HLA-A\*02:01 | 9 | 28.972028 | 28.972053 | 8.6e-7 |
| CINGVCWTV / HLA-A\*02:01 | 9 | 92.105940 | 92.105926 | 1.5e-7 |
| **GILGFVFTLAAA** / HLA-A\*02:01 | **12** | **25274.910033** | **25274.911516** | **5.9e-8** |
| **GILGFVFTLAAAA** / HLA-A\*02:01 | **13** | **32389.125801** | **32389.127574** | **5.5e-8** |
| **GILGFVFTLAAAAA** / HLA-A\*02:01 | **14** | **32972.178346** | **32972.173041** | **1.6e-7** |

All within the unit's claimed RelTol = 1e-3 (<0.1%); observed gap ≤ 9e-7 (≪ the asserted <0.03% in the report). The **len-12/13/14 goldens are independently confirmed**, and the **8/9/10-mer goldens are byte-identical to before** (no regression). `MhcflurryAffinityPredictor.cs` was *not* touched by commit 66c24491 (only 3 `[TestCase]` rows were added to the test file), so the forward-pass/encoding code is bit-for-bit unchanged for the 8–11 lengths — confirmed both by the unchanged source and by the oracle reproducing the old 8/9/10-mer values.

### Stage B — code path & test-quality gate
- **Length gate:** `IsValidPeptideLength` (`OncologyAnalyzer.cs:8268`) is inclusive `[8,14]` for class I; the gate test `IsValidPeptideLength_ClassI_RespectsEightToFourteen` asserts 8/9/11/14 valid, 7 too short, **15** above the max — exactly the validated 8–14 window with 15 rejected.
- **Generator:** `GenerateNeoantigenPeptides` default `maxLength = MhcClassIMaxPeptideLength = 14`; the generator test asserts 7 lengths × 5 windows = **35** candidates (was 20 for 8–11), with a per-length count loop covering 12/13/14, and last peptide = a 14-mer.
- **Encoder boundaries:** `EncodePeptide_OutOfRangeLength_Throws` covers length 5 (min) and 15 (max) encode, length 4 and 16 throw — the MHCflurry encoder's 5–15 range.
- **Classification chain (unchanged):** the strict-`<` nested-threshold model (Strong<50 / Weak<500 nM; %Rank class I 0.5/2, class II 2/10) and the predict→classify chain all re-confirmed against the IEDB/Sette/NetMHCpan cutoffs; goldens are exact oracle values with citation comments, not code echoes. A broken encoder/forward-pass would fail the new 12/13/14 cases (they assert exact oracle IC50s within 1e-3); a `<`→`<=` threshold flip still fails the classification boundary cases. No green-washing, no skips.

**Stage B: PASS.** Full unfiltered suite `dotnet test Seqeron.sln -c Debug` = **18861 passed / 0 failed** (Genomics project); no files changed this session.

### End-state
**✅ CLEAN.** The length-window widening 8–11 → 8–14 is correct and fully sourced (NetMHCpan-4.1 class I window, Reynisson 2020, retrieved this session); lengths 12/13/14 are oracle-verified against `mhcflurry` 2.1.5's exact embedded network to ≤ 9e-7 relative error; the 8/9/10-mer predictions are byte-identical (no regression). The classification thresholds and chain are unchanged and still correct. The prior 8–11 PASS-WITH-NOTES divergence is resolved → Stage A is now a clean PASS. No code or test defect found; no change made this session. **Residual (honest, by-design):** the full 10-network MHCflurry ensemble weights (~80 MB) are caller-loaded by `LoadWeightPack` (one ~4.6 MB member embedded only as the CI parity fixture); the proprietary NetMHCpan-4.1 ANN and the SMM/BIMAS matrices remain caller-supplied / out of scope — an accepted open boundary, not a defect.

---

## 2026-06-25 — Fresh re-validation (post-reset; predictors split out to MHC-MATRIX-001 / MHC-NN-001)

ONCO-MHC-001 was reset to ⬜ pending after the SMM/BIMAS matrix predictor and the MHCflurry NN were carved out into their own already-CLEAN units (MHC-MATRIX-001, MHC-NN-001). This unit's OWN canonical surface — the IC50/%Rank binding-affinity **classification** thresholds and the predict→classify **chain** — was re-validated fresh against externally retrieved sources THIS session. No prior repo TestSpec/Evidence/test was trusted.

### Stage A — sources re-retrieved this session
- **NetMHCpan-4.1 / NetMHCIIpan-4.0 (Reynisson et al. 2020, *Nucleic Acids Res* 48(W1):W449–W454)** — WebSearch confirmed verbatim: class I "strong binder if the %Rank is below … (by default, 0.5%), and a weak binder if … below … (by default, 2%)"; class II SB < 2% / WB < 10%; %Rank defined as the position of the score in a distribution from random natural peptides. Confirms the four %Rank cutoffs (0.5 / 2 / 2 / 10) and the strict-`<` "below" semantics.
- **IEDB threshold guidance** (help.iedb.org) — WebSearch confirmed verbatim: "peptides with IC50 values <50 nM are considered high affinity, <500 nM intermediate affinity and <5000 nM low affinity"; "the threshold of IC50 = 500 nM or 50 nM selects peptide binders or strong binders, respectively." Confirms the IC50 cutoffs 50 / 500 nM and that 500 nM is the binder/non-binder demarcation, all strict `<`.
- **Sette et al. (1994)** and class II length 13–25 (IEDB tool description) were confirmed in the prior session and are unchanged; the two primary conventions above were independently re-grounded this session.

The classification model — `Strong iff v < 50`, `Weak iff 50 ≤ v < 500`, `NonBinder iff v ≥ 500` (IC50); analogous per-class strict-`<` nesting for %Rank — matches the verbatim "<"/"below" inequalities in every source. **Stage A: PASS-WITH-NOTES** (the only divergence is the *documented* class I 8–11 default vs the full 8–14 NetMHCpan admits; see original §Note 1).

### Independent boundary hand-evaluation (recorded labels — derived from sources, NOT code output)

`ClassifyBindingAffinity(ic50)` — validate finite & > 0, then `< 50 → Strong`, else `< 500 → Weak : NonBinder`:

| IC50 (nM) | Hand-derived label (IEDB/Sette strict `<`) | Code trace |
|---|---|---|
| 10 | Strong | Strong ✅ |
| 49.9 | Strong (49.9 < 50) | Strong ✅ |
| 50 (boundary) | Weak (not `< 50`; 50 < 500) | Weak ✅ |
| 500 (boundary) | NonBinder (not `< 500`) | NonBinder ✅ |
| 5000 | NonBinder | NonBinder ✅ |

`ClassifyBindingRank` boundaries (class I 0.5/2, class II 2/10), all strict `<`: 0.4/I→Strong, 0.5/I→Weak, 2.0/I→NonBinder; 1.5/II→Strong, 2.0/II→Weak, 10.0/II→NonBinder. All match the code.

**Predict→classify chain end-to-end** (`PredictAndClassifySmm` → `ClassifyBindingAffinity`, SMM transform `IC50 = 50000^(1−score)` from MHC-MATRIX-001):
- GILGFVFTL (influenza M1 58–66, paradigm HLA-A*02:01 binder), contributions summing to score 1.0 → IC50 = 50000⁰ = **1 nM → Strong**.
- score 0.6 → IC50 = 50000⁰·⁴ = 75.786 nM ∈ [50,500) → **Weak**.
- poly-W 9-mer, score 0 → IC50 = 50000¹ = **50000 nM → NonBinder**.
Chain composes correctly: predictor yields IC50, classifier applies the sourced thresholds. Binder IC50 (1 nM) ≪ non-binder IC50 (50000 nM), correct ranking.

### Stage B — implementation & test-quality gate
- Code path re-read (`OncologyAnalyzer.cs:8197–8287`, chain `:8526–8531`): all four classifiers realise the strict-`<` nested-threshold model exactly; the chain delegates the predicted IC50 into `ClassifyBindingAffinity`. **Stage B: PASS.**
- Test fixture `OncologyAnalyzer_ClassifyMhcBinding_Tests.cs` (47 tests) re-audited against the HARD gate: every classification assertion is an exact `Is.EqualTo`/`Throws<>` carrying the published cutoff value in its citation comment — sourced, not code-echoed. A `<`→`<=` flip fails M2/M4/M7/M9/M11b/M12; swapping class I/II cutoffs fails M11/M11b/M13. Coverage is complete on the unit's own surface: every IC50 boundary (50/500) and %Rank boundary (0.5/2 class I, 2/10 class II), all three bands (Strong/Weak/NonBinder), the length gate (class I 8/11/12, class II 13/25/26, non-positive), all invalid-input paths (IC50 ≤0/NaN/±∞, %Rank out-of-range/NaN, null peptide, empty matrix), and the predict→classify chain in all three bands (P9 Strong, P10b Weak, P10 NonBinder). No green-washing, no skips; the only inequality assertions are ranking sanity in the predictor (MHC-MATRIX-001) tests, paired with exact value checks.

### End-state
**✅ CLEAN.** No code or description defect on ONCO-MHC-001's own surface; thresholds, categories and strict-`<` boundary handling all trace to IEDB/NetMHCpan/Sette retrieved this session; the predict→classify chain is correct and tested end-to-end. **No code/test change made this session** (the surface was already complete and correctly sourced). Unit fixture **47 passed / 0 failed**; the prior full-suite baseline (18213 green) is unchanged since no files were touched. The trained predictor matrices being caller-supplied (MHC-MATRIX-001) and the MHCflurry weight pack being the separate MHC-NN-001 are the documented, accepted scope boundaries.

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

## 2026-06-25 — MHCflurry Class I pan-allele binding-AFFINITY neural predictor ported (Apache-2.0)

### What was added
A faithful C# port of MHCflurry's Class I pan-allele binding-affinity network in a new class `MhcflurryAffinityPredictor` (Oncology project):
- `EncodePeptide` — BLOSUM62 `left_pad_centered_right_pad` (3×15×21 = 945), per `encodable_sequences.py`.
- `EncodePseudosequence` / `GetPseudosequence` / `GetAllelePseudosequences` — the bundled 37-residue allele pseudosequence table (`Resources/mhcflurry.allele_sequences.csv`, Apache-2.0, 14 993 alleles).
- `ToIc50` — `IC50 = 50000^(1−x)` (`regression_target.to_ic50`).
- `LoadWeightPack` + `Network.ForwardRaw / PredictIc50` — the feed-forward pass (`tanh` hidden, `sigmoid` output) supporting both MHCflurry topologies (`feedforward` and `with-skip-connections`).
- `PredictIc50(networks, …)` — geometric-mean ensemble (`exp(mean(log(ic50)))`, `ensemble_centrality.py`).
- `PredictAndClassify(networks, …)` — chains the predicted IC50 into the existing `OncologyAnalyzer.ClassifyBindingAffinity`.

### Evidence retrieved this session (no citation from memory)
Installed `mhcflurry` 2.1.5 (`pip3 install --user mhcflurry`) and fetched `models_class1_pan` (release 20200610) via `mhcflurry-downloads`. Read the actual source modules (`amino_acid.py`, `encodable_sequences.py`, `allele_encoding.py`, `class1_neural_network.py`, `regression_target.py`, `ensemble_centrality.py`, `class1_affinity_predictor.py`) and the model `manifest.csv` (`network_json` Keras graphs) + `weights_*.npz` + `allele_sequences.csv`. The `LICENSE` file (`mhcflurry-2.1.5.dist-info/licenses/LICENSE`) was opened and confirmed to be the full **Apache License 2.0**; the NOTICE/attribution is preserved in `Resources/MHCFLURRY_NOTICE.txt`.

### Oracle cross-check (independent of the C# code)
Re-ran the MHCflurry Python API in-session. A from-scratch numpy reimplementation of the forward pass reproduced the model exactly; the C# port matches:
- **Single-network** (smallest member, embedded as the CI test fixture): 8 peptide/allele pairs (SIINFEKL 11483.20, GILGFVFTL 19.12, NLVPMVATV 17.54, ELAGIGILTV 119.05, AAAWYLWEV 16.56, SIINFEKL/B\*07:02 28830.80, SLYNTVATL 28.97, CINGVCWTV 92.11 nM) — C# within < 0.1% (observed < 0.03%).
- **Full 10-network ensemble** (geometric mean): SIINFEKL/HLA-A\*02:01 = 11927.16 nM, GILGFVFTL = 19.96 nM, etc. — numpy/C# engine reproduces to < 0.03% (verified in-session).
- Predict→classify chain: GILGFVFTL/HLA-A\*02:01 → Strong; SIINFEKL/HLA-A\*02:01 → NonBinder.

Note: the README's SIINFEKL ≈ 6029 nM refers to a different/older model snapshot (the default *presentation* model); the authoritative oracle here is the exact installed `models_class1_pan` (20200610), which gives 11927 nM, and the C# port matches it.

### Residual (honest)
The MHCflurry trained ensemble weights total ≈ 80 MB of near-incompressible float32 across the 10 networks. Embedding 80 MB would be the single largest data artifact in the repo, so it is **not embedded** for repo health: the pseudosequence table + forward-pass engine ship bundled, **one** ensemble member (~4.6 MB) is embedded only as the CI parity test fixture, and the full ensemble is loaded from a caller-supplied MHCflurry weight pack via `LoadWeightPack` (the user fetches the Apache-2.0 weights they already have via `mhcflurry-downloads`). The algorithm itself is exact and oracle-verified — only the weight payload is caller-supplied (analogous to ONCO-IMMUNE-001's LM22). The *latest* NetMHCpan-4.1 ANN and the MHCflurry **presentation/processing** heads remain out of scope.

### End-state
`dotnet build` 0 warnings/0 errors in changed files; the new fixture `MhcflurryAffinityPredictor_PredictIc50_Tests.cs` 19 passed / 0 failed; full unfiltered suite green (see commit). Classification and defaults unchanged; the neural predictor is opt-in. Validation Status confirmed **☐** (unchanged); Quick-Reference counts unchanged.

## Runtime enforcement (LimitationPolicy)

This unit's guarded branch — matrix (SMM/BIMAS) pMHC scoring from a caller-supplied matrix (no bundled validated matrix) — has **minimum access mode `Moderate`** (`Seqeron.Genomics.Core.LimitationCatalog`). Under the default `LimitationPolicy.DefaultMode = Moderate` it is **allowed** (this guarded branch throws only under `Strict`); see [LIMITATIONS.md](../LIMITATIONS.md) › Runtime enforcement. Additive policy layer; the validated contract and `✅ CLEAN` verdict are unchanged.
