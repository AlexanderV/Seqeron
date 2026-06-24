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
5. Class I peptide length: 8–14, default 8–11 — Reynisson 2020. Class II: 13–25 — IEDB (#5).
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
| `IsValidPeptideLength(int length, MhcClass mhcClass)` | OncologyAnalyzer | Canonical | class I 8–11, class II 13–25 |
| `ClassifyMhcBinding(int peptideLength, double ic50Nm, MhcClass mhcClass)` | OncologyAnalyzer | Delegate | length gate + `ClassifyBindingAffinity` |
| `PredictIc50Smm(string peptide, PmhcScoringMatrix matrix)` | OncologyAnalyzer | Canonical | SMM sum → `IC50 = 50000^(1−score)` |
| `PredictBindingHalfLifeBimas(string peptide, PmhcScoringMatrix matrix)` | OncologyAnalyzer | Canonical | BIMAS product → finalConstant·∏ coefficients |
| `PredictAndClassifySmm(string peptide, PmhcScoringMatrix matrix)` | OncologyAnalyzer | Delegate | `PredictIc50Smm` + `ClassifyBindingAffinity` |
| `LoadScoringMatrix(IEnumerable<string> lines)` | OncologyAnalyzer | Canonical | caller-supplied matrix loader |

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
| M14 | Length class I valid | len = 9, class I | true | Reynisson 2020 (8–11) |
| M15 | Length class I too short | len = 7, class I | false | Reynisson 2020 |
| M16 | Length class I above range | len = 12, class I | false | Reynisson 2020 (default 8–11) |
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

- Classification half (M1–C1) implemented in the 2026-06-14 session (27 tests). The 2026-06-25 extension adds the matrix-based predictor (P1–P18); those cases start ❌ Missing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M21, M21b | ✅ Covered | classification half (prior session) |
| S1–S5 | ✅ Covered | classification half (prior session) |
| C1 | ✅ Covered | classification half (prior session) |
| P1–P18 (prediction) | ❌ Missing | matrix-based predictor (this session) |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMhcBinding_Tests.cs` — all cases for this unit.
- **Remove:** none (new unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_ClassifyMhcBinding_Tests.cs | canonical | 43 (27 classification + 16 prediction) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | P1–P5 (SMM transform) | ❌ Missing | implemented | ✅ Done |
| 2 | P6–P8 (SMM validation) | ❌ Missing | implemented | ✅ Done |
| 3 | P9–P10 (predict→classify, ranking) | ❌ Missing | implemented | ✅ Done |
| 4 | P11–P13 (BIMAS product) | ❌ Missing | implemented | ✅ Done |
| 5 | P14–P18 (loader) | ❌ Missing | implemented | ✅ Done |

**Total items (this session):** 5 groups (16 prediction tests; the 27 classification tests were done in the prior session)
**✅ Done:** 5 | **⛔ Blocked:** 0 | **Remaining:** 0

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

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Class I accepted length range = 8–11 (Reynisson default; matches ONCO-NEO-001 constants) rather than the full 8–14 | M16, `IsValidPeptideLength` class I |
| 2 | The coefficient MATRIX is caller-supplied — no redistributable trained matrix was obtainable (BIMAS CGI dead/unarchived; Parker 1994 paywalled; IEDB SMM non-commercial). Only the published scoring RULES are embedded. | P1–P18, `Predict*`, `LoadScoringMatrix` |

---

## 7. Open Questions / Decisions

1. The PAN-ALLELE NetMHCpan neural prediction remains **out of scope**. An opt-in matrix-based predictor (BIMAS product / SMM transform) is now provided; the trained coefficient matrix is caller-supplied because no redistributable, cross-verifiable matrix could be retrieved (the scoring rules ARE fully sourced and cross-verifiable; the weights are not embedded). No model weights are fabricated. The existing classification + defaults are unchanged.
