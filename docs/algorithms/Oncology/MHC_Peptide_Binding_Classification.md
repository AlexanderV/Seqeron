# MHC-Peptide Binding Classification

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Neoantigen prediction |
| Test Unit ID | ONCO-MHC-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-25 |

## 1. Overview

This unit (a) classifies a peptide–MHC pair into binding-strength categories (Strong / Weak / NonBinder) from a predicted binding affinity (IC50 in nM) or %Rank and validates the peptide length, and (b) provides an **opt-in matrix-based predictor** that turns a position-specific scoring matrix into a predicted affinity: the BIMAS / Parker (1994) product rule (predicted half-time of dissociation) [6][7] and the SMM / Peters & Sette (2005) transform `IC50 = 50000^(1−score)` [8][9], which chains into the same classifier. The threshold model is specification-driven with no learned parameters [1][2][3]. The trained coefficient **matrix is caller-supplied** — no redistributable, cross-verifiable trained HLA matrix was obtainable (BIMAS coefficient files are served by a now-defunct dynamic CGI; the Parker 1994 table is paywalled; IEDB SMM matrices are non-commercial) [6][7][8] — so only the published scoring *rules* are bundled (hence Implementation Status **Framework**). The pan-allele NetMHCpan neural model remains out of scope [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Neoantigen prediction (ONCO-NEO-001) enumerates candidate mutant 8–11-mer peptides. Whether such a peptide is presented depends on its binding affinity to a patient HLA allele. Affinity is predicted as an IC50 (half-maximal inhibitory concentration, lower = stronger) or as a %Rank (percentile of the score against random natural peptides, lower = stronger) [1]. Downstream selection then thresholds these predicted values into binder categories using community-standard cutoffs [1][2][4].

### 2.2 Core Model

Given a predicted value `v` and per-class strong/weak cutoffs `(s, w)` with `s < w`, the category is:

- Strong  iff `v < s`
- Weak    iff `s ≤ v < w`
- NonBinder iff `v ≥ w`

Cutoffs (all strict `<`):

- IC50 (nM): `s = 50`, `w = 500` — Sette et al. 1994 "≈500 nM (preferably 50 nM or less)" [2]; IEDB "<50 nM high affinity, <500 nM intermediate affinity" [4]; 500 nM binder demarcation corroborated by Roomp et al. 2010 [3].
- %Rank class I: `s = 0.5`, `w = 2` — Reynisson et al. 2020 "by default, %Rank < 0.5% and %Rank < 2% ... for SBs and WBs for class I" [1].
- %Rank class II: `s = 2`, `w = 10` — Reynisson et al. 2020 "%Rank < 2% and %Rank < 10%, for SBs and WBs for class II" [1].

Peptide-length validity (inclusive bounds):

- Class I: 8–11 — Reynisson et al. 2020 ("8 to 14 amino acids, default is 8–11"), adopting the 8–11 default to match the ONCO-NEO-001 windowing constants [1].
- Class II: 13–25 — IEDB MHC class II tool description ("typically range between 13 and 25 amino acids long") [5].

**Matrix-based prediction (opt-in).** A position-specific scoring matrix has one residue→value row per peptide position plus a scalar final constant. Two published scoring conventions are supported:

- **BIMAS / Parker (1994)** — the predicted half-time of dissociation is the **product** of the position-specific coefficients times a final constant: `T½ = FinalConstant · ∏_i coeff[i][peptide[i]]`. The running score starts at 1.0; an unlisted/ambiguous residue contributes the neutral coefficient 1.0 (no effect). Higher T½ = stronger binder. Source: BIMAS scoring documentation ("The initial (running) score is set to 1.0 … multiplied by the coefficient … multiplied by a final constant to yield an estimate of the half time of disassociation") [7]; Parker 1994 ("calculated by multiplying together the corresponding coefficients") [6].
- **SMM / Peters & Sette (2005)** — the log50k score is the **sum** of position-specific additive contributions plus an intercept; the IC50 is recovered by `IC50 = 50000^(1 − score)`, the algebraic inverse of the IEDB linearisation `log50k = 1 − log(IC50)/log(50000)` [8][9]. Lower IC50 = stronger binder. An unlisted residue contributes 0 (additive identity).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | IC50 input must be a finite value > 0 | IC50 is a positive concentration; Registry invariant "IC50 > 0"; the method throws otherwise |
| INV-02 | %Rank input must be in [0, 100] | %Rank is a percentile of scores from random peptides [1]; the method throws otherwise |
| INV-03 | Categories are monotone in the score: a smaller IC50/%Rank yields an equal-or-stronger category | `s < w` and nested strict thresholds |
| INV-04 | Cutoffs are strict `<`; a value exactly at a cutoff falls in the weaker category | Verbatim "<" in [1] and IEDB tiers [4] |
| INV-05 | SMM `IC50 = 50000^(1−score)` is strictly decreasing in score; score 0→50000, 1→1, 0.5→√50000; always finite & > 0 | Algebraic inverse of IEDB log50k [9]; Peters & Sette 2005 [8] |
| INV-06 | BIMAS `T½ = FinalConstant · ∏ coefficients`; an unlisted residue contributes the neutral coefficient 1.0 | BIMAS scoring docs [7]; Parker 1994 [6] |
| INV-07 | `Predict*` require a non-empty matrix and `peptide.Length == matrix.Rows.Count`; else `ArgumentException` (null peptide → `ArgumentNullException`) | Implementation contract (one row per position) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `ic50Nm` | `double` | required | Predicted IC50 in nanomolar | finite, > 0 |
| `percentRank` | `double` | required | Predicted %Rank percentile | [0, 100] |
| `mhcClass` | `MhcClass` | required | ClassI or ClassII | enum |
| `peptideLength` | `int` | required | Peptide residue count | any int (out-of-range ⇒ invalid) |
| `peptide` | `string` | required | Peptide for matrix scoring | length must equal `matrix.Rows.Count` |
| `matrix` | `PmhcScoringMatrix` | required | Caller-supplied scoring matrix (per-position residue→value rows + final constant) | ≥ 1 row |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `BindingStrength` | Strong / Weak / NonBinder |
| return | `bool` | (`IsValidPeptideLength`) whether length is in the class range |
| return | `double` | (`PredictBindingHalfLifeBimas`) predicted half-time of dissociation (BIMAS units); (`PredictIc50Smm`) predicted IC50 (nM) |
| return | `(double Ic50Nm, BindingStrength Strength)` | (`PredictAndClassifySmm`) predicted IC50 + its category |
| return | `PmhcScoringMatrix` | (`LoadScoringMatrix`) parsed caller-supplied matrix |

### 3.3 Preconditions and Validation

`ClassifyBindingAffinity` throws `ArgumentOutOfRangeException` for NaN, ∞, or `ic50Nm ≤ 0`. `ClassifyBindingRank` throws `ArgumentOutOfRangeException` for NaN or `percentRank ∉ [0,100]`. `IsValidPeptideLength` is total (returns `false` for any out-of-range or non-positive length). `ClassifyMhcBinding` returns `NonBinder` for an invalid length before evaluating affinity (so an invalid-length peptide is never reported as a binder), and otherwise propagates `ClassifyBindingAffinity` validation. The matrix predictors throw `ArgumentNullException` for a null peptide and `ArgumentException` when the matrix has no rows or the peptide length differs from the row count (INV-07); `LoadScoringMatrix` throws `ArgumentNullException` for null input and `FormatException` for a malformed token, a non-numeric value, or a multi-character residue key.

## 4. Algorithm

### 4.1 High-Level Steps

1. (Optional, `ClassifyMhcBinding`) Reject peptides whose length is outside the MHC class range → NonBinder.
2. Validate the supplied value (IC50 > 0, or %Rank ∈ [0,100]).
3. Compare against the strong cutoff, then the weak cutoff, using strict `<`.
4. Return Strong / Weak / NonBinder.

### 4.2 Decision Rules, Scoring, Reference Tables

| Class | Metric | Strong (`<`) | Weak (`<`) | Source |
|-------|--------|-------------|-----------|--------|
| I/II | IC50 (nM) | 50 | 500 | [2][4][3] |
| I | %Rank | 0.5 | 2 | [1] |
| II | %Rank | 2 | 10 | [1] |
| I | length | — | — | 8–11 inclusive [1] |
| II | length | — | — | 13–25 inclusive [5] |

Matrix-based prediction (opt-in):

| Convention | Combine rule | Predicted quantity | Source |
|------------|--------------|--------------------|--------|
| BIMAS | `T½ = FinalConstant · ∏ coeff` (missing residue = 1.0) | half-time of dissociation (higher = stronger) | [6][7] |
| SMM | `score = intercept + Σ contrib` (missing residue = 0); `IC50 = 50000^(1 − score)` | IC50 nM (lower = stronger) → `ClassifyBindingAffinity` | [8][9] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| each classification / validation | O(1) | O(1) | two comparisons |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ClassifyBindingAffinity(double)`: IC50 (nM) → BindingStrength.
- `OncologyAnalyzer.ClassifyBindingRank(double, MhcClass)`: %Rank → BindingStrength (per-class cutoffs).
- `OncologyAnalyzer.IsValidPeptideLength(int, MhcClass)`: length range check.
- `OncologyAnalyzer.ClassifyMhcBinding(int, double, MhcClass)`: length gate + affinity classification.
- `OncologyAnalyzer.PredictBindingHalfLifeBimas(string, PmhcScoringMatrix)`: BIMAS product rule → half-time of dissociation.
- `OncologyAnalyzer.PredictIc50Smm(string, PmhcScoringMatrix)`: SMM sum + `50000^(1−score)` → IC50 (nM).
- `OncologyAnalyzer.PredictAndClassifySmm(string, PmhcScoringMatrix)`: `PredictIc50Smm` → `ClassifyBindingAffinity` (predict→classify chain).
- `OncologyAnalyzer.LoadScoringMatrix(IEnumerable<string>)`: parse a caller-supplied matrix (`CONST=…` + `RESIDUE=VALUE` rows).
- `OncologyAnalyzer.PmhcScoringMatrix` / `PmhcScoringMethod`: the matrix record and scoring-convention enum.

### 5.2 Current Behavior

The threshold half is a pure O(1) comparison; no search/matching, so the repository suffix tree is not applicable. The matrix predictor is a single linear pass over the peptide (product for BIMAS, sum for SMM). No trained model weights are embedded: the coefficient matrix is caller-supplied via `LoadScoringMatrix` (or constructed directly), exactly as ONCO-IMMUNE-001 handles the CIBERSORT LM22 signature matrix. All cutoffs and the SMM base (`SmmIc50Base = 50000`) are named source-cited constants. Default behaviour is unchanged: the predictors are opt-in additions; existing callers of the classifiers see no change.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- IC50 cutoffs 50 nM (strong) / 500 nM (weak), strict `<` [2][4][3].
- %Rank cutoffs class I 0.5/2, class II 2/10, strict `<` [1].
- Length ranges class I 8–11, class II 13–25 [1][5].
- BIMAS product rule `T½ = FinalConstant · ∏ coeff`, neutral coefficient 1.0 [6][7].
- SMM transform `IC50 = 50000^(1 − score)`, score = intercept + Σ contributions [8][9].

**Intentionally simplified:**

- Class I length range uses the 8–11 default rather than the full 8–14 [1]; **consequence:** lengths 12–14 are reported invalid for class I, matching the ONCO-NEO-001 canonical neoantigen search.
- The coefficient **matrix is caller-supplied**, not embedded; **consequence:** the user must obtain a matrix (the public-domain BIMAS/Parker 1994 HLA-A2 table, or an IEDB SMM matrix under its non-commercial licence) under their own licence — no redistributable, cross-verifiable trained matrix was obtainable this session.

**Not implemented:**

- The pan-allele NetMHCpan / MHCflurry **neural** prediction; **users should rely on:** an external predictor (NetMHCpan / MHCflurry / IEDB) for affinities outside the matrix model, or supply a coefficient matrix for the BIMAS/SMM predictors here.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Class I length range 8–11 (vs 8–14) | Assumption | Lengths 12–14 reported invalid for class I | accepted | Reynisson default; matches ONCO-NEO-001 |
| 2 | Coefficient matrix caller-supplied (no embedded trained matrix) | Assumption | User must obtain a licensed matrix; only scoring rules are bundled | accepted | BIMAS CGI dead/unarchived; Parker 1994 paywalled; IEDB SMM non-commercial |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| IC50 = 50 nM | Weak | strict `<` strong cutoff [4] / INV-04 |
| IC50 = 500 nM | NonBinder | strict `<` weak cutoff [4] / INV-04 |
| %Rank = 0.5 (class I) | Weak | strict `<` [1] |
| %Rank = 2.0 (class I) | NonBinder | strict `<` [1] |
| IC50 ≤ 0 / NaN / ∞ | ArgumentOutOfRangeException | INV-01 |
| %Rank ∉ [0,100] / NaN | ArgumentOutOfRangeException | INV-02 |
| length out of class range | invalid / NonBinder | length gate [1][5] |

### 6.2 Limitations

The matrix-based predictor needs a caller-supplied coefficient matrix (no trained matrix is embedded; the pan-allele NetMHCpan neural model is out of scope). The BIMAS half-life and SMM IC50 are only as good as the supplied matrix and the independent-side-chain approximation (Parker 1994 reports accuracy "to within a factor of 5" [6]). Class I length default 8–11 excludes 12–14. %Rank cutoffs are the NetMHCpan-4.1 defaults; allele-specific cutoffs are not modeled.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// A caller-supplied NetMHCpan IC50 of 42 nM for a 9-mer on an HLA-A class I allele:
var s = OncologyAnalyzer.ClassifyMhcBinding(peptideLength: 9, ic50Nm: 42.0, OncologyAnalyzer.MhcClass.ClassI);
// s == BindingStrength.Strong  (length 9 valid; 42 < 50)

var r = OncologyAnalyzer.ClassifyBindingRank(0.5, OncologyAnalyzer.MhcClass.ClassI);
// r == BindingStrength.Weak  (0.5 is NOT < 0.5; 0.5 < 2)

// Opt-in matrix-based prediction (caller supplies the coefficient matrix under their own licence):
var matrix = OncologyAnalyzer.LoadScoringMatrix(new[] { "CONST=0.0", /* one RESIDUE=VALUE line per position */ });
var (ic50, strength) = OncologyAnalyzer.PredictAndClassifySmm("GILGFVFTL", matrix);
// SMM: IC50 = 50000^(1 - score); chained into ClassifyBindingAffinity → Strong/Weak/NonBinder.
```

**Numerical walk-through:** an SMM matrix whose contributions for `GILGFVFTL` sum to `score = 1.0` (intercept 0) gives `IC50 = 50000^(1 − 1) = 1 nM` → Strong; a peptide matching none of the listed residues scores `0` → `IC50 = 50000^(1 − 0) = 50000 nM` → NonBinder. A BIMAS matrix with constant 10 and coefficients 2.0, 3.0, 1.5 for `LMV` gives `T½ = 10 · (2·3·1.5) = 90`.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_ClassifyMhcBinding_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMhcBinding_Tests.cs) — covers `INV-01`–`INV-07`
- Evidence: [ONCO-MHC-001-Evidence.md](../../../docs/Evidence/ONCO-MHC-001-Evidence.md)
- Related algorithms: [Neoantigen peptide windowing (ONCO-NEO-001)](../../Evidence/ONCO-NEO-001-Evidence.md)

## 8. References

1. Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. 2020. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research* 48(W1):W449–W454. https://doi.org/10.1093/nar/gkaa379
2. Sette A, Vitiello A, Reherman B, et al. 1994. The relationship between class I binding affinity and immunogenicity of potential cytotoxic T cell epitopes. *Journal of Immunology* 153(12):5586–5592. https://pubmed.ncbi.nlm.nih.gov/7527444/
3. Roomp K, Antes I, Lengauer T. 2010. Predicting MHC class I epitopes in large datasets. *BMC Bioinformatics* 11:90. https://doi.org/10.1186/1471-2105-11-90
4. IEDB. What thresholds (cut-offs) should I use for MHC class I and II binding predictions. Accessed 2026-06-14. https://help.iedb.org/hc/en-us/articles/114094152371-What-thresholds-cut-offs-should-I-use-for-MHC-class-I-and-II-binding-predictions
5. IEDB. T Cell Epitopes - MHC Class II Binding Prediction Tools Description. Accessed 2026-06-14. https://help.iedb.org/hc/en-us/articles/114094151731-T-Cell-Epitopes-MHC-Class-II-Binding-Prediction-Tools-Description
6. Parker KC, Bednarek MA, Coligan JE. 1994. Scheme for ranking potential HLA-A2 binding peptides based on independent binding of individual peptide side-chains. *Journal of Immunology* 152(1):163–175. https://pubmed.ncbi.nlm.nih.gov/8254189/
7. BIMAS — Information & background on the HLA peptide motif searches (NIH/CIT/CBEL; R. Taylor & K. Parker). Accessed 2026-06-25 (server retired; via Internet Archive). https://web.archive.org/web/20041016022153/http://www-bimas.cit.nih.gov/molbio/hla_bind/hla_motif_search_info.html
8. Peters B, Sette A. 2005. Generating quantitative models describing the sequence specificity of biological processes with the stabilized matrix method. *BMC Bioinformatics* 6:132. https://doi.org/10.1186/1471-2105-6-132
9. IEDB MHC class I log50k linearisation `log50k = 1 − log(IC50)/log(50000)` (restated in: D. Farrell, "Create an MHC-Class I binding predictor in Python"). Accessed 2026-06-25. https://dmnfarrell.github.io/bioinformatics/mhclearning
