# Cancer-Specific Variant Annotation (AMP/ASCO/CAP Tiers)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-ANNOT-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Classifies a somatic sequence variant into the four-tier clinical-significance system of the AMP/ASCO/CAP 2017 guideline [1]: Tier I (strong clinical significance), Tier II (potential clinical significance), Tier III (unknown clinical significance), and Tier IV (benign / likely benign). The classification is specification-driven: it applies the decision criteria of Figure 2 / Tables 3–7 of the guideline [1] to caller-supplied evidence (assigned clinical evidence level, population minor allele frequency, and whether the variant has a published cancer association). It is used to interpret and report somatic variants found in tumor sequencing. The library does not reproduce curated knowledgebases (population, somatic, or germline databases); those facts are supplied by the caller, which is why the implementation status is **Framework**.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Interpretation of somatic variants in cancer focuses on clinical impact (does the variant predict sensitivity, resistance, diagnosis, or prognosis?) rather than germline pathogenicity/causality [1]. The AMP/ASCO/CAP working group standardized this into four tiers ranked by the strength of clinical and experimental evidence, drawing on professional guidelines, population databases (to exclude polymorphisms), somatic databases (COSMIC, TCGA), germline databases, and the literature [1].

### 2.2 Core Model

The guideline defines four **levels of evidence** [1, Table 3] and maps them to tiers [1, Figure 2]:

- **Level A** — biomarkers predicting response/resistance to FDA-approved therapies for a specific tumor type, or in professional guidelines [1].
- **Level B** — based on well-powered studies with expert consensus [1].
- **Level C** — FDA/guideline therapies for a *different* tumor type (off-label), clinical-trial inclusion criteria, or diagnostic/prognostic from multiple small studies [1].
- **Level D** — plausible therapeutic significance from preclinical studies, or diagnostic/prognostic from small studies / few case reports without consensus [1].

Tier assignment [1, Figure 2]:

- **Tier I** — strong clinical significance: **Level A or B** evidence [1].
- **Tier II** — potential clinical significance: **Level C or D** evidence [1].
- **Tier III** — unknown clinical significance: not observed at a significant allele frequency, with no convincing published evidence of cancer association [1, Table 6].
- **Tier IV** — benign / likely benign: observed at a significant allele frequency (**MAF ≥ 1%** in the general population), or with no published evidence of cancer association [1, Table 7, Figure 2].

The benign population-frequency cutoff is the guideline's recommended primary cutoff: "In the absence of paired normal tissue, the work group recommends using 1% (0.01) as a primary cutoff" for eliminating polymorphic/benign variants [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every variant maps to exactly one of the four tiers. | The four tiers are exhaustive and the decision rule is a total, deterministic cascade [1, Figure 2]. |
| INV-02 | Level A/B ⇒ Tier I and Level C/D ⇒ Tier II, independent of MAF / cancer association. | The guideline categorizes by evidence level first [1, Figure 2]. |
| INV-03 | With no evidence level, MAF ≥ 0.01 OR no cancer association ⇒ Tier IV; otherwise Tier III. | Tier IV criteria (MAF ≥ 1% / no cancer association) [1, Table 7, Figure 2]; Tier III otherwise [1, Table 6]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| variant | `CancerVariantAnnotationInput` | required | Gene, protein change, evidence level, population MAF, cancer-association flag | `PopulationMaf ∈ [0, 1]` |
| variants | `IEnumerable<CancerVariantAnnotationInput>` | required | Batch of variant evidence records | non-null |
| cosmicCatalog | `IReadOnlyDictionary<(string Gene, string ProteinChange), string>` | required | Caller-supplied COSMIC records (e.g. COSMIC IDs) | non-null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (tier) | `VariantTier` | One of Tier I–IV (`ClassifyVariantTier`). |
| annotations | `IReadOnlyList<CancerVariantAnnotation>` | One `(variant, tier)` per input, in input order (`AnnotateCancerVariants`). |
| (cosmic id) | `string?` | Catalog value on a hit, `null` on a miss (`GetCOSMICAnnotation`). |

### 3.3 Preconditions and Validation

`PopulationMaf` must be a real number in `[0, 1]`; NaN, negative, or `> 1` throws `ArgumentOutOfRangeException`. Null `variants` or null `cosmicCatalog` throws `ArgumentNullException`. An empty `variants` batch returns an empty list. Gene / protein-change strings are compared by exact (ordinal) equality for the COSMIC lookup.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `PopulationMaf ∈ [0, 1]`.
2. If evidence level is A or B ⇒ Tier I.
3. Else if evidence level is C or D ⇒ Tier II.
4. Else (no evidence level): if `MAF ≥ 0.01` OR no cancer association ⇒ Tier IV.
5. Else ⇒ Tier III.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

| Evidence level | Population MAF | Cancer association | Tier |
|----------------|----------------|--------------------|------|
| A or B | any | any | Tier I [1, Fig 2] |
| C or D | any | any | Tier II [1, Fig 2] |
| None | ≥ 0.01 | any | Tier IV [1, Table 7] |
| None | < 0.01 | false | Tier IV [1, Fig 2] |
| None | < 0.01 | true | Tier III [1, Table 6] |

Constant: `BenignPopulationMafThreshold = 0.01` (1% primary cutoff [1]).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ClassifyVariantTier` | O(1) | O(1) | constant-time cascade. |
| `AnnotateCancerVariants` (n variants) | O(n) | O(n) | one pass; n annotations. |
| `GetCOSMICAnnotation` | O(1) avg | O(1) | dictionary lookup. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ClassifyVariantTier(CancerVariantAnnotationInput)`: per-variant tier decision.
- `OncologyAnalyzer.AnnotateCancerVariants(IEnumerable<CancerVariantAnnotationInput>)`: batch annotation, input order.
- `OncologyAnalyzer.GetCOSMICAnnotation(CancerVariantAnnotationInput, IReadOnlyDictionary<(string,string),string>)`: COSMIC catalog lookup.

### 5.2 Current Behavior

Clinical evidence (Tier I/II) is evaluated before the benign-frequency rule: a Level A/B biomarker remains Tier I even if it is also common in population databases, matching the guideline's evidence-level-first categorization [1, Figure 2]. The COSMIC lookup uses a caller-supplied dictionary (mirroring the existing `MatchCancerHotspots` caller-supplied-set pattern in this class) rather than an embedded database. **Search reuse:** the repository suffix tree was evaluated and is not applicable — this unit performs constant-time enum/threshold decisions and a hash-map lookup, not substring/occurrence search.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Four-tier system and Level→Tier mapping: Tier I = Level A/B, Tier II = Level C/D, Tier III = unknown, Tier IV = benign/likely benign [1, Figure 2].
- Benign population cutoff MAF ≥ 1% (0.01) [1, Population Databases / Table 7].
- Tier III vs Tier IV discrimination by population frequency and cancer association [1, Tables 6/7, Figure 2].

**Intentionally simplified:**

- Evidence level, population MAF, and cancer-association flag are reduced to a single caller-supplied record; **consequence:** the user must perform the underlying database/literature lookups (the guideline lists 12 evidence pieces [1]) and pass the resulting summary, rather than the library curating them.

**Not implemented:**

- Curated population/somatic/germline databases (gnomAD, COSMIC, ClinVar, HGMD) and in-silico predictors (SIFT, PolyPhen2, CADD); **users should rely on:** external annotation pipelines to populate the evidence inputs and the caller-supplied COSMIC catalog.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Caller-supplied evidence inputs | Assumption | User must supply database-derived facts | accepted | The guideline's tiering rule is applied verbatim to supplied evidence. |
| 2 | Tier III/IV discriminator (MAF ≥ 0.01 OR no cancer assoc.) | Assumption | Determines III vs IV when no evidence level | accepted | Direct reading of Figure 2 boxes / Tables 6–7 [1]. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Level A but high MAF | Tier I | Categorized by evidence level [1, Fig 2]. |
| No evidence, MAF = 0.01 exactly | Tier IV | Cutoff is "≥ 1%" (inclusive) [1, Table 7]. |
| No evidence, MAF = 0.0099, cancer assoc. | Tier III | Below cutoff, has association [1, Table 6]. |
| No evidence, low MAF, no cancer assoc. | Tier IV | "No published evidence of cancer association" [1, Fig 2]. |
| MAF NaN / < 0 / > 1 | `ArgumentOutOfRangeException` | input validation. |
| COSMIC catalog miss | `null` | external content; do not fabricate [2]. |

### 6.2 Limitations

The classification is only as good as the caller-supplied evidence; it does not itself determine evidence levels, query databases, run in-silico predictors, or distinguish somatic vs germline origin. It does not implement the ClinGen/CGC/VICC oncogenicity standard (a separate guideline). Tumor-type context (a variant can be Tier I in one tumor and Tier II off-label in another) must be reflected by the caller in the supplied evidence level.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var input = new OncologyAnalyzer.CancerVariantAnnotationInput(
    Gene: "BRAF", ProteinChange: "p.V600E",
    EvidenceLevel: OncologyAnalyzer.ClinicalEvidenceLevel.A,
    PopulationMaf: 0.0, HasCancerAssociation: true);

OncologyAnalyzer.VariantTier tier = OncologyAnalyzer.ClassifyVariantTier(input);
// tier == VariantTier.TierI_StrongClinicalSignificance (Level A ⇒ Tier I)

var cosmic = new Dictionary<(string, string), string> { [("BRAF", "p.V600E")] = "COSV56056643" };
string? id = OncologyAnalyzer.GetCOSMICAnnotation(input, cosmic); // "COSV56056643"
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_AnnotateCancerVariants_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnnotateCancerVariants_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`
- Evidence: [ONCO-ANNOT-001-Evidence.md](../../../docs/Evidence/ONCO-ANNOT-001-Evidence.md)
- Related algorithms: [Somatic_Mutation_Calling](../Oncology/Somatic_Mutation_Calling.md), [Driver_Mutation_Detection](../Oncology/Driver_Mutation_Detection.md)

## 8. References

1. Li MM, Datto M, Duncavage EJ, Kulkarni S, Lindeman NI, Roy S, Tsimberidou AM, Vnencak-Jones CL, Wolff DJ, Younes A, Nikiforova MN. 2017. Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists. Journal of Molecular Diagnostics 19(1):4–23. https://doi.org/10.1016/j.jmoldx.2016.10.002
2. Tate JG, Bamford S, Jubb HC, et al. 2019. COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Research 47(D1):D941–D947. https://doi.org/10.1093/nar/gky1015
