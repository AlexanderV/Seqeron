# MHC-Peptide Binding Classification

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Neoantigen prediction |
| Test Unit ID | ONCO-MHC-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-14 |

## 1. Overview

This unit classifies a peptide–MHC pair into binding-strength categories (Strong / Weak / NonBinder) from a **caller-supplied** predicted binding affinity (IC50 in nM) or %Rank, and validates the peptide length for the chosen MHC class. It is the classification half of MHC-binding prediction: the actual peptide–MHC affinity prediction requires a trained model (e.g. NetMHCpan) and is out of scope / caller-supplied [1]. The algorithm is specification-driven — it applies the standard IEDB / NetMHCpan thresholds exactly, with no learned parameters [1][2][3].

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

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | IC50 input must be a finite value > 0 | IC50 is a positive concentration; Registry invariant "IC50 > 0"; the method throws otherwise |
| INV-02 | %Rank input must be in [0, 100] | %Rank is a percentile of scores from random peptides [1]; the method throws otherwise |
| INV-03 | Categories are monotone in the score: a smaller IC50/%Rank yields an equal-or-stronger category | `s < w` and nested strict thresholds |
| INV-04 | Cutoffs are strict `<`; a value exactly at a cutoff falls in the weaker category | Verbatim "<" in [1] and IEDB tiers [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `ic50Nm` | `double` | required | Predicted IC50 in nanomolar | finite, > 0 |
| `percentRank` | `double` | required | Predicted %Rank percentile | [0, 100] |
| `mhcClass` | `MhcClass` | required | ClassI or ClassII | enum |
| `peptideLength` | `int` | required | Peptide residue count | any int (out-of-range ⇒ invalid) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `BindingStrength` | Strong / Weak / NonBinder |
| return | `bool` | (`IsValidPeptideLength`) whether length is in the class range |

### 3.3 Preconditions and Validation

`ClassifyBindingAffinity` throws `ArgumentOutOfRangeException` for NaN, ∞, or `ic50Nm ≤ 0`. `ClassifyBindingRank` throws `ArgumentOutOfRangeException` for NaN or `percentRank ∉ [0,100]`. `IsValidPeptideLength` is total (returns `false` for any out-of-range or non-positive length). `ClassifyMhcBinding` returns `NonBinder` for an invalid length before evaluating affinity (so an invalid-length peptide is never reported as a binder), and otherwise propagates `ClassifyBindingAffinity` validation.

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

### 5.2 Current Behavior

Pure threshold comparison; no search/matching, so the repository suffix tree is not applicable. No trained model, PSSM, or learned weights are present: prediction is the caller's responsibility. All cutoffs are named source-cited constants (`StrongBinderIc50Nm`, `ClassIStrongBinderRankPercent`, etc.).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- IC50 cutoffs 50 nM (strong) / 500 nM (weak), strict `<` [2][4][3].
- %Rank cutoffs class I 0.5/2, class II 2/10, strict `<` [1].
- Length ranges class I 8–11, class II 13–25 [1][5].

**Intentionally simplified:**

- Class I length range uses the 8–11 default rather than the full 8–14 [1]; **consequence:** lengths 12–14 are reported invalid for class I, matching the ONCO-NEO-001 canonical neoantigen search.

**Not implemented:**

- Peptide–MHC affinity / %Rank **prediction** (the trained NetMHCpan model and PSSM weights); **users should rely on:** an external predictor (NetMHCpan / MHCflurry / IEDB) supplying the IC50 or %Rank that these methods classify.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Class I length range 8–11 (vs 8–14) | Assumption | Lengths 12–14 reported invalid for class I | accepted | Reynisson default; matches ONCO-NEO-001 |

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

Does not predict affinity (no model). Class I length default 8–11 excludes 12–14. %Rank cutoffs are the NetMHCpan-4.1 defaults; allele-specific cutoffs are not modeled.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// A caller-supplied NetMHCpan IC50 of 42 nM for a 9-mer on an HLA-A class I allele:
var s = OncologyAnalyzer.ClassifyMhcBinding(peptideLength: 9, ic50Nm: 42.0, OncologyAnalyzer.MhcClass.ClassI);
// s == BindingStrength.Strong  (length 9 valid; 42 < 50)

var r = OncologyAnalyzer.ClassifyBindingRank(0.5, OncologyAnalyzer.MhcClass.ClassI);
// r == BindingStrength.Weak  (0.5 is NOT < 0.5; 0.5 < 2)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_ClassifyMhcBinding_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMhcBinding_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-MHC-001-Evidence.md](../../../docs/Evidence/ONCO-MHC-001-Evidence.md)
- Related algorithms: [Neoantigen peptide windowing (ONCO-NEO-001)](../../Evidence/ONCO-NEO-001-Evidence.md)

## 8. References

1. Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. 2020. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research* 48(W1):W449–W454. https://doi.org/10.1093/nar/gkaa379
2. Sette A, Vitiello A, Reherman B, et al. 1994. The relationship between class I binding affinity and immunogenicity of potential cytotoxic T cell epitopes. *Journal of Immunology* 153(12):5586–5592. https://pubmed.ncbi.nlm.nih.gov/7527444/
3. Roomp K, Antes I, Lengauer T. 2010. Predicting MHC class I epitopes in large datasets. *BMC Bioinformatics* 11:90. https://doi.org/10.1186/1471-2105-11-90
4. IEDB. What thresholds (cut-offs) should I use for MHC class I and II binding predictions. Accessed 2026-06-14. https://help.iedb.org/hc/en-us/articles/114094152371-What-thresholds-cut-offs-should-I-use-for-MHC-class-I-and-II-binding-predictions
5. IEDB. T Cell Epitopes - MHC Class II Binding Prediction Tools Description. Accessed 2026-06-14. https://help.iedb.org/hc/en-us/articles/114094151731-T-Cell-Epitopes-MHC-Class-II-Binding-Prediction-Tools-Description
