# HLA Allele Nomenclature Parsing and Allele-Specific HLA LOH

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-HLA-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-15 |

## 1. Overview

This unit provides two retrievable, formally specified pieces of HLA analysis: (1) parsing and validation of HLA allele names against the WHO HLA Nomenclature (`HLA-A*02:01` etc.), and (2) classification of allele-specific HLA loss of heterozygosity (HLA LOH) from caller-supplied per-allele copy number and allelic-imbalance evidence, per the LOHHLA method [1][2][3]. Both are specification-/decision-rule-driven and exact: parsing is a deterministic validator of the nomenclature grammar, and LOH classification is a deterministic threshold rule. The unit deliberately does **not** perform HLA genotyping from reads (no trained model / reference allele database); genotype and copy number are taken as caller-supplied inputs.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The HLA (human leukocyte antigen) genes are the most polymorphic in the human genome. The WHO Nomenclature Committee for Factors of the HLA System assigns each allele a name of the form `HLA-<Gene>*Field1:Field2[:Field3[:Field4]]` with an optional expression-status suffix [1][2]. In cancer, tumors can lose one HLA haplotype (allele-specific HLA LOH), removing the ability to present a subset of neoantigens; LOHHLA quantifies this from sequencing by estimating the copy number of each homologous HLA allele [3].

### 2.2 Core Model

**Nomenclature grammar** [1][2]: an allele name is `HLA-` + gene + `*` + 2–4 colon-separated numeric fields + optional suffix. Field 1 is the type/allele group ("often corresponds to the serological antigen") [1]; Field 2 the specific HLA protein/subtype [1]; Field 3 synonymous coding-region substitutions [1]; Field 4 non-coding (intron / UTR) differences [1]. "All alleles receive at least a four digit name … the first two sets of digits" → minimum two fields, maximum four [1][2]. The expression suffix ∈ {N, L, S, C, A, Q} [1].

**HLA LOH decision rule** [3]: for each HLA gene, LOHHLA infers the allele-specific copy number of both homologous alleles. "A copy number < 0.5, is classified as subject to loss, and thereby indicative of LOH." [3]. To avoid over-calling, "Allelic imbalance is determined if p < 0.01 using the paired Student's t-Test." [3]. Thus HLA LOH is called iff one allele has CN < 0.5 **and** the allelic-imbalance p value is < 0.01.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Both alleles CN < 0.5 with significant imbalance is reported as homozygous loss, not allele-specific LOH (the source defines a single lost allele as CN < 0.5) [3]. | A genuine homozygous deletion would otherwise be mislabeled as allele-specific LOH. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A parsed allele has 2 ≤ field count ≤ 4 | Validator rejects <2 or >4 fields per [1][2]. |
| INV-02 | The normalized `Name` round-trips the parsed gene/fields/suffix | `Name` is reconstructed from the parsed components. |
| INV-03 | HLA LOH called ⇔ (exactly one allele CN < 0.5) ∧ (imbalance p < 0.01) | Direct encoding of the two LOHHLA thresholds [3]. |
| INV-04 | When LOH is called, the lost allele is the one with CN < 0.5 | The branch selects the sub-threshold allele [3]. |
| INV-05 | CN = 0.5 is retained, p = 0.01 is not significant | Both comparisons are strict (`< 0.5`, `< 0.01`) per the verbatim source wording [3]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| alleleName | string | required | HLA allele name to parse | `HLA-` prefix; `*` gene/field separator; 2–4 digit fields; optional N/L/S/C/A/Q suffix |
| alleleCopyNumber | HlaAlleleCopyNumber | required | Per-allele CN + imbalance p | CN ≥ 0; p ∈ [0, 1] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| HlaAllele | record | Gene (upper-cased), Fields (2–4 digit strings), Suffix; `Name` = normalized string |
| HlaLohResult | record | IsLoh, LostAllele ∈ {None, Allele1, Allele2, Both}, AllelicImbalanceSignificant |

### 3.3 Preconditions and Validation

`ParseHlaAllele` throws `ArgumentNullException` for null, `ArgumentException` for empty/whitespace, and `FormatException` for any grammar violation (missing prefix/`*`, wrong field count, non-numeric field, invalid suffix). `TryParseHlaAllele` returns false instead of throwing for null/format/argument errors. The gene is upper-cased; the `HLA-` prefix and a single suffix letter are matched case-insensitively. `DetectHlaLoh` throws `ArgumentException` for negative copy number or a p value outside [0, 1]. Copy-number and p-value comparisons are strict.

## 4. Algorithm

### 4.1 High-Level Steps

1. **Parse:** require `HLA-` prefix; split gene from fields at `*`; strip an optional trailing expression letter; split remaining block on `:`; validate 2–4 numeric fields.
2. **LOH:** validate inputs; test allelic-imbalance significance (p < 0.01); test each allele's CN against 0.5; emit LOH (one allele lost + significant), homozygous loss (both lost), or no LOH.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Field-count bounds: min 2, max 4 [1][2].
- Expression suffix map: N→Null, L→Low, S→Secreted, C→Cytoplasm, A→Aberrant, Q→Questionable [1].
- LOH thresholds: allele CN < 0.5 (strict) [3]; allelic-imbalance paired t-test p < 0.01 (strict) [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ParseHlaAllele | O(n) | O(f) | n = name length, f = field count (≤ 4). |
| DetectHlaLoh | O(1) | O(1) | Constant-time threshold rule. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ParseHlaAllele(string)`: parse + validate a WHO HLA nomenclature name.
- `OncologyAnalyzer.TryParseHlaAllele(string?, out HlaAllele)`: non-throwing wrapper.
- `OncologyAnalyzer.DetectHlaLoh(HlaAlleleCopyNumber)`: LOHHLA allele-specific LOH classification.

### 5.2 Current Behavior

The expression suffix is detected only when the final character of the field block is a letter, so purely numeric names parse normally. The `HLA-` prefix and the suffix letter are matched case-insensitively and the gene is upper-cased; field digit strings are preserved verbatim (leading zeros kept, e.g. `02`). No substring/pattern search is performed, so the repository suffix tree is not applicable (see §7.3 note).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Nomenclature grammar `HLA-<Gene>*F1:F2[:F3[:F4]][suffix]`, 2–4 fields, suffix set {N,L,S,C,A,Q} [1][2].
- HLA LOH rule: allele CN < 0.5 AND allelic-imbalance paired-t-test p < 0.01 [3].

**Intentionally simplified:**

- (none)

**Not implemented:**

- HLA genotyping / allele calling from reads (LOHHLA's coverage/BAF copy-number estimation and the IPD-IMGT/HLA reference alignment); **users should rely on:** an external HLA caller (e.g. LOHHLA, OptiType, Polysolver) to supply genotype and per-allele copy number, which this unit then validates / classifies.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Both alleles CN < 0.5 → homozygous loss label | Assumption | Distinguishes homozygous deletion from allele-specific LOH | accepted | ASM-01; thresholds unchanged |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `HLA-A*02` (single field) | FormatException | Two-field minimum [1][2] |
| `HLA-A*02:01:01:01:01` (five fields) | FormatException | Four-field maximum [1][2] |
| `HLA-A*02:01X` | FormatException | X ∉ {N,L,S,C,A,Q} [1] |
| `A*02:01` | FormatException | Missing `HLA-` prefix [1] |
| CN exactly 0.5 | Not lost | Strict `< 0.5` [3] |
| Imbalance p exactly 0.01 | Not significant → no LOH | Strict `< 0.01` [3] |
| Low CN but p ≥ 0.01 | No LOH | Over-calling guard [3] |

### 6.2 Limitations

Does not genotype HLA alleles or estimate copy number from sequencing reads; both are caller-supplied. The nomenclature validator checks structural grammar, not membership in the IPD-IMGT/HLA database (a syntactically valid but non-existent allele is accepted). Homozygous loci (two identical alleles) cannot be assessed for allele-specific loss because no polymorphic sites distinguish the homologs [3].

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var allele = OncologyAnalyzer.ParseHlaAllele("HLA-A*24:02:01:02L");
// allele.Gene == "A"; allele.Fields == ["24","02","01","02"]; allele.Suffix == HlaExpressionSuffix.Low

var loh = OncologyAnalyzer.DetectHlaLoh(
    new OncologyAnalyzer.HlaAlleleCopyNumber("HLA-A*02:01", 1.8, "HLA-A*11:01", 0.30, 0.001));
// loh.IsLoh == true; loh.LostAllele == HlaLostAllele.Allele2
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_HlaAnalysis_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_HlaAnalysis_Tests.cs) — covers `INV-01`…`INV-05`
- Evidence: [ONCO-HLA-001-Evidence.md](../../../docs/Evidence/ONCO-HLA-001-Evidence.md)
- Search reuse: the repository suffix tree was evaluated and is **not used** — this unit performs no substring/pattern search (nomenclature parsing is a single linear validation; LOH is a constant-time threshold rule).

## 8. References

1. WHO Nomenclature Committee for Factors of the HLA System. Naming Alleles. IPD-IMGT/HLA. https://hla.alleles.org/pages/nomenclature/naming_alleles/
2. Marsh SGE, Albert ED, Bodmer WF, et al. 2010. Nomenclature for factors of the HLA system, 2010. Tissue Antigens 75(4):291–455. https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1399-0039.2010.01466.x
3. McGranahan N, Rosenthal R, Hiley CT, et al. 2017. Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution. Cell 171(6):1259–1271. https://pmc.ncbi.nlm.nih.gov/articles/PMC5720478/
4. mskcc/lohhla. LOHHLAscript.R (reference implementation). https://raw.githubusercontent.com/mskcc/lohhla/master/LOHHLAscript.R
