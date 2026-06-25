# Test Specification: ONCO-MHC-001

**Test Unit ID:** ONCO-MHC-001
**Area:** Oncology
**Algorithm:** MHC-Peptide Binding (length/affinity/%rank classification + matrix-based prediction)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Reynisson et al. (2020) NetMHCpan-4.1, *Nucleic Acids Res* 48(W1):W449–W454 | 1 | https://doi.org/10.1093/nar/gkaa379 (PMC7319546) | 2026-06-14 |
| 2 | Sette et al. (1994), *J Immunol* 153(12):5586–92 | 1 | https://pubmed.ncbi.nlm.nih.gov/7527444/ | 2026-06-14 |
| 3 | Roomp, Antes & Lengauer (2010), *BMC Bioinformatics* 11:90 | 1 | https://doi.org/10.1186/1471-2105-11-90 (PMC2836306) | 2026-06-14 |
| 4 | IEDB threshold help article | 5 | https://help.iedb.org/hc/en-us/articles/114094152371 | 2026-06-14 |
| 5 | IEDB MHC class II tool description | 5 | https://help.iedb.org/hc/en-us/articles/114094151731 | 2026-06-14 |
| 6 | Parker, Bednarek & Coligan (1994), *J Immunol* 152(1):163–175 (BIMAS product rule) | 1 | https://pubmed.ncbi.nlm.nih.gov/8254189/ | 2026-06-25 |
| 7 | BIMAS HLA peptide-motif-search scoring docs (archived) | 5 | https://web.archive.org/web/20041016022153/http://www-bimas.cit.nih.gov/molbio/hla_bind/hla_motif_search_info.html | 2026-06-25 |
| 8 | Peters & Sette (2005), *BMC Bioinformatics* 6:132 (SMM) | 1 | https://doi.org/10.1186/1471-2105-6-132 | 2026-06-25 |
| 9 | IEDB log50k transform `log50k = 1 − log(IC50)/log(50000)` | 3 | https://dmnfarrell.github.io/bioinformatics/mhclearning | 2026-06-25 |

### 1.2 Key Evidence Points

1. Class I %Rank: strong binder < 0.5%, weak binder < 2% (default) — Reynisson 2020 (PMC7319546).
2. Class II %Rank: strong binder < 2%, weak binder < 10% (default) — Reynisson 2020.
3. IC50 tiers: high (strong) < 50 nM, intermediate (weak) < 500 nM, low < 5000 nM — IEDB (#4); the 50/500 nM cutoffs trace to Sette 1994 "≈500 nM (preferably 50 nM or less)".
4. 500 nM is the binder/non-binder demarcation — Roomp 2010 (PMC2836306), corroborating IEDB.
5. Class I peptide length: **8–14** (the full NetMHCpan-4.1 class I window; the service offers 8/9/10/11/12/13/14-mer options) — Reynisson 2020. Class II: 13–25 — IEDB (#5). (Within the MHCflurry encoding's 5–15 support, `max_length = 15`.)
6. The PAN-ALLELE prediction (NetMHCpan neural model) is out of scope. An opt-in MATRIX-BASED predictor is now provided: the BIMAS product rule and the SMM IC50 transform (below). The trained coefficient MATRIX is caller-supplied (no redistributable matrix was obtainable).
7. BIMAS scoring (verbatim): running score starts at 1.0, is multiplied by each position's coefficient, then by a final constant → estimated half-time of dissociation (HLA-A2); unlisted residue coefficient = 1.0 — source #7; Parker 1994 (#6) "calculated by multiplying together the corresponding coefficients".
8. SMM/IEDB transform (verbatim): `log50k = 1 − log(IC50)/log(50000)` ⇒ `IC50 = 50000^(1 − score)`; the score is the sum of position-specific additive contributions plus an intercept — sources #8, #9.

### 1.3 Documented Corner Cases

- Boundary semantics are strict `<`: IC50 = 50 nM is NOT strong (→ weak); IC50 = 500 nM is NOT a binder (→ non-binder). %Rank = 0.5 is NOT strong (→ weak); %Rank = 2.0 is NOT a weak binder (→ non-binder). Source: verbatim "<" inequalities in Reynisson 2020 and IEDB tiers.
- Peptide length outside the class range is not a valid binder candidate (Reynisson 2020; IEDB #5).

### 1.4 Known Failure Modes / Pitfalls

1. IC50 must be a positive concentration (invariant IC50 > 0); zero/negative/NaN/∞ are invalid — Registry invariant; concentration semantics.
2. %Rank is a percentile in [0, 100]; values < 0, > 100, or NaN are invalid — percentile definition (Reynisson 2020 §%Rank).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ClassifyBindingAffinity(double ic50Nm)` | OncologyAnalyzer | Canonical | IC50 → Strong/Weak/NonBinder (50/500 nM) |
| `ClassifyBindingRank(double percentRank, MhcClass mhcClass)` | OncologyAnalyzer | Canonical | %Rank → Strong/Weak/NonBinder (class I 0.5/2; class II 2/10) |
| `IsValidPeptideLength(int length, MhcClass mhcClass)` | OncologyAnalyzer | Canonical | class I 8–14, class II 13–25 |
| `ClassifyMhcBinding(int peptideLength, double ic50Nm, MhcClass mhcClass)` | OncologyAnalyzer | Delegate | length gate + `ClassifyBindingAffinity` |
| `PredictIc50Smm(string peptide, PmhcScoringMatrix matrix)` | OncologyAnalyzer | Canonical | SMM sum → `IC50 = 50000^(1−score)` |
| `PredictBindingHalfLifeBimas(string peptide, PmhcScoringMatrix matrix)` | OncologyAnalyzer | Canonical | BIMAS product → finalConstant·∏ coefficients |
| `PredictAndClassifySmm(string peptide, PmhcScoringMatrix matrix)` | OncologyAnalyzer | Delegate | `PredictIc50Smm` + `ClassifyBindingAffinity` |
| `LoadScoringMatrix(IEnumerable<string> lines)` | OncologyAnalyzer | Canonical | caller-supplied matrix loader |
| `EncodePeptide(string)` | MhcflurryAffinityPredictor | Canonical | BLOSUM62 `left_pad_centered_right_pad` (945) |
| `EncodePseudosequence(string)` / `GetPseudosequence(string)` | MhcflurryAffinityPredictor | Canonical | bundled 37-residue allele pseudosequence table |
| `ToIc50(double)` | MhcflurryAffinityPredictor | Canonical | `IC50 = 50000^(1−x)` |
| `Network.ForwardRaw / PredictIc50` + `LoadWeightPack(Stream)` | MhcflurryAffinityPredictor | Canonical | per-network feed-forward pass (both topologies) |
| `PredictIc50(networks, peptide, allele)` | MhcflurryAffinityPredictor | Canonical | geometric-mean ensemble IC50 |
| `PredictAndClassify(networks, peptide, allele)` | MhcflurryAffinityPredictor | Delegate | predict → `ClassifyBindingAffinity` chain |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | IC50 input must be > 0 (positive concentration); else exception | Yes | Registry invariant "IC50 > 0"; concentration semantics |
| INV-2 | %Rank input must be in [0, 100]; else exception | Yes | Reynisson 2020 (%Rank is a percentile) |
| INV-3 | Strong ⇒ Weak ⇒ NonBinder are mutually exclusive, monotone in the score (smaller IC50/%Rank ⇒ stronger or equal class) | Yes | Reynisson 2020; IEDB tiers |
| INV-4 | Classification cutoffs are strict `<` (boundary value excluded from the stronger class) | Yes | Verbatim "<" in Reynisson 2020 and IEDB tiers |
| INV-5 | SMM: `IC50 = 50000^(1−score)` — strictly decreasing in score; score 0→50000, 1→1, 0.5→√50000; IC50 always finite & > 0 | Yes | IEDB log50k (#9); Peters & Sette 2005 (#8) |
| INV-6 | BIMAS: `T½ = finalConstant · ∏ coefficients`; an unlisted residue contributes the neutral coefficient 1.0 | Yes | BIMAS scoring docs (#7); Parker 1994 (#6) |
| INV-7 | `Predict*` require `peptide.Length == matrix.Rows.Count` and a non-empty matrix; else `ArgumentException` (null peptide → `ArgumentNullException`) | Yes | implementation contract (one row per peptide position) |
| INV-8 | MHCflurry peptide encoding: length ∈ [5, 15]; out-of-range → `ArgumentOutOfRangeException`; output length 3·15·21 = 945 | Yes | `encodable_sequences.py` (#11) |
| INV-9 | MHCflurry `to_ic50(x) = 50000^(1−x)`: x=0→50000, x=1→1, x=0.5→√50000 | Yes | `regression_target.to_ic50` (#11) |
| INV-10 | MHCflurry ensemble = geometric mean of per-network IC50s (`exp(mean(log(ic50)))`); duplicating one network reproduces its IC50 | Yes | `ensemble_centrality.py` + `predict_to_dataframe` (#11) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | IC50 strong | ic50 = 10 nM | Strong | IEDB <50; Sette 1994 |
| M2 | IC50 strong boundary | ic50 = 50 nM (strict `<`) | Weak | IEDB "<50"; INV-4 |
| M3 | IC50 weak | ic50 = 200 nM | Weak | IEDB <500 |
| M4 | IC50 weak boundary | ic50 = 500 nM (strict `<`) | NonBinder | IEDB "<500"; Roomp 2010 |
| M5 | IC50 non-binder | ic50 = 1000 nM | NonBinder | IEDB tiers |
| M6 | %Rank class I strong | rank = 0.4, class I | Strong | Reynisson 2020 (<0.5) |
| M7 | %Rank class I strong boundary | rank = 0.5, class I | Weak | Reynisson 2020 ("<0.5"); INV-4 |
| M8 | %Rank class I weak | rank = 1.0, class I | Weak | Reynisson 2020 (<2) |
| M9 | %Rank class I weak boundary | rank = 2.0, class I | NonBinder | Reynisson 2020 ("<2"); INV-4 |
| M10 | %Rank class I non-binder | rank = 5.0, class I | NonBinder | Reynisson 2020 |
| M11 | %Rank class II strong | rank = 1.5, class II | Strong | Reynisson 2020 (<2) |
| M12 | %Rank class II weak boundary | rank = 10.0, class II | NonBinder | Reynisson 2020 ("<10"); INV-4 |
| M13 | %Rank class II weak | rank = 5.0, class II | Weak | Reynisson 2020 (<10) |
| M14 | Length class I valid | len = 8/9/11/14, class I | true | Reynisson 2020 (class I 8–14) |
| M15 | Length class I too short | len = 7, class I | false | Reynisson 2020 |
| M16 | Length class I above range | len = 15, class I | false | Reynisson 2020 (class I window tops at 14) |
| M17 | Length class II valid | len = 15, class II | true | IEDB #5 (13–25) |
| M18 | Length class II too short | len = 12, class II | false | IEDB #5 |
| M19 | Length class II too long | len = 26, class II | false | IEDB #5 |
| M20 | Combined gate fails on length | len = 7, ic50 = 10, class I | NonBinder (invalid length ⇒ not a candidate) | Reynisson 2020 length gate + IEDB affinity |
| M21 | Combined gate passes length, strong affinity | len = 9, ic50 = 10, class I | Strong | combined of M1 + M14 |
| M21b | Combined gate, invalid IC50 on valid length | len = 9, ic50 = 0, class I | ArgumentOutOfRangeException (INV-1 propagates) | INV-1 |
| P1 | SMM score 0 | 1-mer, score 0 | IC50 = 50000 nM | INV-5; IEDB log50k |
| P2 | SMM score 1 | 1-mer, contribution 1 | IC50 = 1 nM | INV-5 |
| P3 | SMM score 0.5 | 1-mer, contribution 0.5 | IC50 = √50000 = 223.6067977499790 nM | INV-5 |
| P4 | SMM intercept summed | intercept 0.3 + contribution 0.2 = 0.5 | IC50 = √50000 | INV-5 (intercept is part of the sum) |
| P5 | SMM unlisted residue | residue absent from row | contributes 0 ⇒ IC50 = 50000 | INV-6 analogue (additive identity) |
| P6 | SMM null peptide | null | ArgumentNullException | INV-7 |
| P7 | SMM length mismatch | len 2 vs 1 row | ArgumentException | INV-7 |
| P8 | SMM empty matrix | 0 rows | ArgumentException | INV-7 |
| P9 | Predict→classify strong | GILGFVFTL, score 1.0 | IC50 = 1 nM, Strong | INV-5 + classification |
| P10 | Strong-vs-non-binder ranking | GILGFVFTL vs poly-W | binder IC50 (1) ≪ non-binder IC50 (50000); Strong vs NonBinder | INV-5 + classification |
| P11 | BIMAS product | const 10 · (2·3·1.5) | T½ = 90.0 | INV-6 |
| P12 | BIMAS unlisted = 1.0 | AAA on same matrix | T½ = 10.0 | INV-6 |
| P13 | BIMAS ranking | favorable vs unfavorable anchors | 20.0 > 0.02 | INV-6 |
| P14 | Loader round-trip | CONST + RESIDUE=VALUE rows | parses 2 rows, const 10, T½(LM)=60 | source #7 (format) |
| P15 | Loader malformed token | "L 2.0" (no '=') | FormatException | loader contract |
| P16 | Loader non-numeric | "L=abc" | FormatException | loader contract |
| P17 | Loader multi-char residue | "LM=2.0" | FormatException | loader contract |
| P18 | Loader null input | null | ArgumentNullException | loader contract |
| MF1 | MHCflurry peptide encoding layout | EncodePeptide("SIINFEKL") | 945 values; positions match `left_pad_centered_right_pad` index layout; pos 0 = BLOSUM62 S-row, pad = X-row | INV-8; `encodable_sequences.py` (#11) |
| MF2 | MHCflurry peptide length bounds | len 5 & 15 ok; 4 & 16 throw; null throws | encode / `ArgumentOutOfRangeException` / `ArgumentNullException` | INV-8 |
| MF3 | Bundled pseudosequence lookup | GetPseudosequence("HLA-A\*02:01") / ("HLA-B\*07:02") | `YFAMYGEKVAHTHVDTLYGVRYDHYYTWAVLAYTWYA` / `YYSEYRNIYAQTDESNLYGLSYDDYYTWAERAYEWYA` (37 res) | `allele_sequences.csv` (#11) |
| MF4 | Unknown/null allele | "HLA-Z\*99:99" / null | KeyNotFoundException / ArgumentNullException | predictor contract |
| MF5 | Table scale | GetAllelePseudosequences() | > 5000 HLA- alleles | `allele_sequences.csv` (#11) |
| MF6 | `to_ic50` anchors | ToIc50(0/1/0.5) | 50000 / 1 / √50000 | INV-9; `regression_target.to_ic50` (#11) |
| MF7 | Single-network oracle parity | 11 peptide/allele pairs (8 of len 8–10 incl. SIINFEKL/A\*02:01, plus len-12/13/14 windows GILGFVFTLAAA / GILGFVFTLAAAA / GILGFVFTLAAAAA) | single-net IC50 within 0.1% of mhcflurry oracle (len 12/13/14 = 25274.910033 / 32389.125801 / 32972.178346 nM) | mhcflurry 2.1.5 oracle (#11) |
| MF8 | Strong vs non-binder ranking | GILGFVFTL vs SIINFEKL on HLA-A\*02:01 | strong < 50 nM; non-binder > 5000 nM; ratio > 100× | mhcflurry oracle (#11) |
| MF9 | Ensemble = geometric mean | duplicate one network ×3 | equals the single-network IC50 | INV-10; `ensemble_centrality.py` (#11) |
| MF10 | Empty ensemble | 0 networks | ArgumentException | predictor contract |
| MF11 | Predict→classify chain | GILGFVFTL (Strong) / SIINFEKL (NonBinder) on HLA-A\*02:01 | Strong / NonBinder via `ClassifyBindingAffinity` | predict→classify (#11 + classification cutoffs) |
| MF12 | Weight-pack loader validation | bad magic / null stream | InvalidDataException / ArgumentNullException | pack-format contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | IC50 ≤ 0 rejected | ic50 = 0 / −1 | ArgumentOutOfRangeException | INV-1 |
| S2 | IC50 non-finite rejected | ic50 = NaN / ∞ | ArgumentOutOfRangeException | INV-1 |
| S3 | %Rank out of [0,100] rejected | rank = −0.1 / 100.1 | ArgumentOutOfRangeException | INV-2 |
| S4 | %Rank NaN rejected | rank = NaN | ArgumentOutOfRangeException | INV-2 |
| S5 | Length ≤ 0 | len = 0 / −1 | false (not valid) | length is a count |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Monotonicity %Rank class I | rank 0.1 < 0.5 ≤ 1 < 2 ≤ 3 ⇒ classes non-increasing in strength | Strong, Weak, Weak, NonBinder, NonBinder | INV-3 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Classification half (M1–C1) implemented in the 2026-06-14 session (27 tests); the matrix-based predictor (P1–P18) in the 2026-06-25 session. The 2026-06-25 MHCflurry extension adds the ported pan-allele neural binding-affinity predictor (MF1–MF12) in a new canonical file; those cases start ❌ Missing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M21, M21b | ✅ Covered | classification half (prior session) |
| S1–S5 | ✅ Covered | classification half (prior session) |
| C1 | ✅ Covered | classification half (prior session) |
| P1–P18 (matrix predictor) | ✅ Covered | matrix-based predictor (prior 2026-06-25 step) |
| MF1–MF12 (MHCflurry neural) | ❌ Missing | ported pan-allele affinity predictor (this step) |

### 5.3 Consolidation Plan

- **Canonical files:** `OncologyAnalyzer_ClassifyMhcBinding_Tests.cs` (classification + matrix predictor) and `MhcflurryAffinityPredictor_PredictIc50_Tests.cs` (MHCflurry neural predictor). Two files because the MHCflurry predictor lives in a separate class (`MhcflurryAffinityPredictor`).
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_ClassifyMhcBinding_Tests.cs | canonical (classification + matrix) | 43 (27 classification + 16 matrix prediction) |
| MhcflurryAffinityPredictor_PredictIc50_Tests.cs | canonical (MHCflurry neural) | 19 (encoders, table, transform, 8 oracle-parity cases, ensemble, chain, loader) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | MF1–MF2 (peptide encoding + bounds) | ❌ Missing | implemented | ✅ Done |
| 2 | MF3–MF5 (pseudosequence table) | ❌ Missing | implemented | ✅ Done |
| 3 | MF6 (`to_ic50`) | ❌ Missing | implemented | ✅ Done |
| 4 | MF7–MF8 (single-net oracle parity + ranking) | ❌ Missing | implemented | ✅ Done |
| 5 | MF9–MF10 (ensemble geometric mean) | ❌ Missing | implemented | ✅ Done |
| 6 | MF11 (predict→classify chain) | ❌ Missing | implemented | ✅ Done |
| 7 | MF12 (weight-pack loader validation) | ❌ Missing | implemented | ✅ Done |

**Total items (this step):** 7 groups (19 MHCflurry tests)
**✅ Done:** 7 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M21, M21b | ✅ Covered | exact expected values from Evidence |
| S1–S5 | ✅ Covered | exception-type assertions |
| C1 | ✅ Covered | monotonicity property |
| P1–P5 (SMM transform) | ✅ Covered | exact IC50 = 50000^(1−score) anchors |
| P6–P8 (SMM validation) | ✅ Covered | null / length-mismatch / empty-matrix exceptions |
| P9–P10 (predict→classify) | ✅ Covered | strong→1 nM Strong; non-binder→50000 nM NonBinder; ranking |
| P11–P13 (BIMAS) | ✅ Covered | exact product·constant; neutral 1.0; ranking |
| P14–P18 (loader) | ✅ Covered | round-trip + malformed/non-numeric/multi-char/null |
| MF1–MF2 (peptide encoding) | ✅ Covered | exact 945-value layout vs `left_pad_centered_right_pad`; [5,15] bounds |
| MF3–MF5 (pseudosequence table) | ✅ Covered | exact 37-residue strings; unknown/null exceptions; table scale |
| MF6 (`to_ic50`) | ✅ Covered | exact 50000^(1−x) anchors |
| MF7–MF8 (oracle parity + ranking) | ✅ Covered | 8 pairs within 0.1% of mhcflurry; strong≪non-binder |
| MF9–MF10 (ensemble geometric mean) | ✅ Covered | duplicate-network identity; empty-ensemble exception |
| MF11 (predict→classify) | ✅ Covered | Strong / NonBinder through `ClassifyBindingAffinity` |
| MF12 (loader validation) | ✅ Covered | bad-magic / null-stream exceptions |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | The coefficient MATRIX is caller-supplied — no redistributable trained matrix was obtainable (BIMAS CGI dead/unarchived; Parker 1994 paywalled; IEDB SMM non-commercial). Only the published scoring RULES are embedded. | P1–P18, `Predict*`, `LoadScoringMatrix` |

> **Resolved 2026-06-26:** the earlier assumption "class I accepted length range = 8–11" is removed. The accepted class I range is now the full NetMHCpan-4.1 window **8–14** (`MhcClassIMaxPeptideLength` 11 → 14), within the MHCflurry encoding's 5–15 support; lengths 12/13/14 are oracle-verified (see MF7).

---

## 7. Open Questions / Decisions

1. The PAN-ALLELE NetMHCpan neural prediction remains **out of scope**. An opt-in matrix-based predictor (BIMAS product / SMM transform) is now provided; the trained coefficient matrix is caller-supplied because no redistributable, cross-verifiable matrix could be retrieved (the scoring rules ARE fully sourced and cross-verifiable; the weights are not embedded). No model weights are fabricated. The existing classification + defaults are unchanged.
