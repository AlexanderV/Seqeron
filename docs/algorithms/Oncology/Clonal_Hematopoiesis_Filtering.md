# Clonal Hematopoiesis (CHIP) Filtering for cfDNA Liquid Biopsy

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-CHIP-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-15 |

## 1. Overview

In plasma cell-free DNA (cfDNA) liquid biopsy, somatic variants arise not only from the tumour but also from clonal hematopoiesis (CH) — expanded blood-cell clones carrying driver mutations. CH is the dominant non-tumour confounder of cfDNA variant calls [3]. This algorithm flags candidate clonal-hematopoiesis-of-indeterminate-potential (CHIP) variants by a gene + variant-allele-fraction (VAF) heuristic [1] and removes CH confounders from a cfDNA call set by matched white-blood-cell (WBC) subtraction [3], leaving candidate tumour-derived variants. It is a deterministic, rule-based filter; the CHIP driver-gene panel is caller-supplied or a labelled canonical default [1][2].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

CHIP is "the presence of detectable somatic clonal mutations in genes recurrently mutated in hematologic malignancies" in individuals who "lack a known hematologic malignancy or other clonal disorder", with "the mutant allele fraction must be ≥2% in the peripheral blood" [1]. CH variants account for 81.6% of cfDNA mutations in controls and 53.2% in cancer patients [3], so distinguishing them from tumour signal is essential for liquid-biopsy interpretation.

### 2.2 Core Model

Let a cfDNA variant be `v = (locus, gene, VAF)` where `locus = (chromosome, position, ref, alt)`.

**CHIP candidate flag (gene + VAF heuristic):**

`isCHIP(v) ⟺ gene(v) ∈ G  ∧  VAF(v) ≥ τ`

where `G` is the CHIP driver-gene panel and `τ = 0.02` (≥ 2%) [1]. The threshold is inclusive [1].

**Matched-WBC subtraction (definitive origin test) [3]:** a cfDNA variant whose locus is also present in the matched WBC sample (≥ 1 alt read [5]) is WBC/CH-derived and removed, *regardless of gene*. `FilterCHIP` removes `v` when it is matched in WBC OR when `isCHIP(v)` holds with no WBC evidence; the complement is retained as candidate tumour.

Canonical driver genes `G` = {DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D}: "Four genes (DNMT3A, TET2, ASXL1, and PPM1D) had disproportionately high numbers of somatic mutations" with recurrent JAK2 V617F and SF3B1 K700E [2]; the three most prevalent are DNMT3A, TET2, ASXL1 [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `isCHIP(v)` ⟺ gene ∈ panel ∧ VAF ≥ τ | direct from the CHIP definition [1] |
| INV-02 | VAF threshold inclusive at 0.02 (≥) | "≥2%" [1] |
| INV-03 | `FilterCHIP` output ⊆ input, input order preserved | filter iterates once, keeps subset |
| INV-04 | A cfDNA variant present in matched WBC is removed regardless of gene | matched-WBC origin test [3] |
| INV-05 | A cfDNA variant absent from matched WBC and not meeting the heuristic is retained | only confounders are removed [3] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Matched-WBC subtraction (this) | Gene+VAF heuristic alone |
|--------|-------------------------------|--------------------------|
| Origin test | definitive (per-locus WBC evidence) [3] | candidate flag only [4] |
| Requires matched WBC | yes (rule a) | no (rule b) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| variants | `IEnumerable<ChipVariant>` | required | cfDNA variants to screen / filter | non-null |
| whiteBloodCellVariants | `IEnumerable<ChipVariant>` | required (FilterCHIP) | matched WBC variants | non-null |
| chipGenes | `IReadOnlyCollection<string>?` | `DefaultChipGenes` | CHIP driver panel (HGNC symbols) | case-insensitive |
| minVaf | `double` | `0.02` | CHIP VAF threshold | (0, 1] |
| minWbcAltReads | `int` | `1` | WBC alt-read cutoff for locus presence | ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| IdentifyCHIPVariants | `IReadOnlyList<ChipVariant>` | candidate CHIP variants, input order |
| FilterCHIP | `IReadOnlyList<ChipVariant>` | retained candidate tumour variants, input order |
| IsCanonicalChipGene | `bool` | gene ∈ panel (case-insensitive) |

### 3.3 Preconditions and Validation

Null `variants` / `whiteBloodCellVariants` → `ArgumentNullException`. `minVaf ∉ (0, 1]` → `ArgumentOutOfRangeException`. `minWbcAltReads < 1` → `ArgumentOutOfRangeException`. Gene comparison is case-insensitive (HGNC symbols are upper-case); a null/empty gene is not a CHIP gene. Locus identity is exact on (chromosome, 1-based position, ref, alt). Empty inputs are valid (empty result / heuristic-only filtering).

## 4. Algorithm

### 4.1 High-Level Steps

1. `IdentifyCHIPVariants`: for each variant, flag when gene ∈ panel and VAF ≥ minVaf.
2. `FilterCHIP`: build a set of WBC loci with ≥ minWbcAltReads alt reads; for each cfDNA variant, drop it if its locus is in that set OR it meets the gene+VAF heuristic; keep the rest.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

- CHIP VAF threshold `τ = 0.02` (≥) [1].
- Default panel {DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D} [1][2].
- WBC presence cutoff = 1 alt read (configurable) [5].
- `HashSet<(string,int,string,string)>` of WBC loci for O(1) membership.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| IdentifyCHIPVariants | O(n·g) | O(n) | n variants, g panel size (small constant) |
| FilterCHIP | O(n + w) | O(w) | w WBC loci hashed once; n cfDNA lookups O(1) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.IdentifyCHIPVariants(...)`: gene+VAF candidate-CHIP flag.
- `OncologyAnalyzer.FilterCHIP(...)`: matched-WBC subtraction + gene+VAF fallback.
- `OncologyAnalyzer.IsCanonicalChipGene(...)`: case-insensitive panel membership.

### 5.2 Current Behavior

The CHIP panel is a labelled canonical default that callers may override (Framework status) — the algorithm is panel-agnostic [3]. Matched-WBC presence is decided by alt-read evidence at an exact locus key, mirroring the repository's MRD `IsVariantDetected` convention [5]. No substring/pattern search is performed (locus equality + gene membership), so the repository suffix tree is **not used** (it does not fit set-membership / exact-key lookups).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- CHIP = driver-gene somatic mutation at VAF ≥ 0.02, inclusive [1].
- Canonical driver-gene panel from Steensma 2015 / Genovese 2014 [1][2].
- Matched-WBC subtraction to remove CH confounders regardless of gene [3].

**Intentionally simplified:**

- WBC presence uses a configurable alt-read cutoff (default ≥ 1); **consequence:** assay-specific error-model cutoffs are the caller's responsibility, the VAF–origin relationship being unclear in general [4].

**Not implemented:**

- Hematologic-malignancy exclusion (cytopenia / morphology / diagnostic criteria) [1]; **users should rely on:** upstream clinical phenotyping — this assay operates on blood/plasma variant evidence only.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Default gene panel is a labelled canonical set | Assumption | overridable; affects which variants flag | accepted | caller may pass `chipGenes` |
| 2 | WBC presence = ≥1 alt read at locus | Assumption | affects subtraction sensitivity | accepted | configurable `minWbcAltReads` [5] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| VAF exactly 0.02 | flagged CHIP | inclusive ≥ [1] |
| Non-CHIP gene, high VAF | not CHIP | driver-gene requirement [1] |
| Variant in matched WBC, non-CHIP gene | removed | matched-WBC origin [3] |
| Empty WBC set | only gene+VAF heuristic applies | nothing to subtract [3] |
| null/empty gene | not a CHIP gene | undefined symbol |

### 6.2 Limitations

Does not assess hematologic-malignancy diagnostic criteria, copy-number/structural CH, or assay-specific error models; the gene+VAF heuristic is a candidate flag, not a definitive origin call (matched WBC is definitive) [3][4].

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var cfDna = new[]
{
    new OncologyAnalyzer.ChipVariant("2", 25234374, "C", "T", "DNMT3A", 0.06),  // CH driver
    new OncologyAnalyzer.ChipVariant("7", 55259515, "T", "G", "EGFR", 0.30),    // tumour candidate
};
var wbc = Array.Empty<OncologyAnalyzer.ChipVariant>();
var tumourCandidates = OncologyAnalyzer.FilterCHIP(cfDna, wbc); // -> only the EGFR variant
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_FilterCHIP_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_FilterCHIP_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [ONCO-CHIP-001-Evidence.md](../../../docs/Evidence/ONCO-CHIP-001-Evidence.md)

## 8. References

1. Steensma DP, Bejar R, Jaiswal S, et al. 2015. Clonal hematopoiesis of indeterminate potential and its distinction from myelodysplastic syndromes. *Blood* 126(1):9–16. https://doi.org/10.1182/blood-2015-03-631747
2. Genovese G, Kähler AK, Handsaker RE, et al. 2014. Clonal Hematopoiesis and Blood-Cancer Risk Inferred from Blood DNA Sequence. *N Engl J Med* 371(26):2477–2487. https://doi.org/10.1056/NEJMoa1409405
3. Razavi P, Li BT, Brown DN, et al. 2019. High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants. *Nat Med* 25:1928–1937. https://doi.org/10.1038/s41591-019-0652-7
4. Arango-Argoty G, et al. 2025. An artificial intelligence-based model for prediction of clonal hematopoiesis variants in cell-free DNA samples. *NPJ Precis Oncol* 9:147. https://doi.org/10.1038/s41698-025-00921-w
5. Wan JCM, Heider K, Gale D, et al. 2020. ctDNA monitoring using patient-specific sequencing and integration of variant reads. *Sci Transl Med* 12(548):eaaz8084. https://doi.org/10.1126/scitranslmed.aaz8084
