# Clinical Actionability Assessment (OncoKB Therapeutic Levels of Evidence)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-ACTION-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-15 |

## 1. Overview

This algorithm classifies the *therapeutic actionability* of a somatic variant by reducing the variant's
curated biomarker–drug associations to their highest OncoKB therapeutic level of evidence. OncoKB stratifies
treatment implications into an ordered set of levels — sensitivity Levels 1, 2, 3A, 3B, 4 and resistance
Levels R1, R2 — by how strongly an alteration predicts response (or resistance) to a therapy [1][2]. The
operation is a specification-driven ranking: it deterministically selects the maximum level under a fixed
precedence order. It is distinct from AMP/ASCO/CAP variant tiering (ONCO-ANNOT-001), which classifies
clinical significance rather than therapeutic actionability. The OncoKB knowledgebase content itself is
caller-supplied; this library ranks the supplied levels and does not embed the curated database.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Precision oncology matches tumor alterations to targeted therapies. OncoKB (the Precision Oncology Knowledge
Base, Memorial Sloan Kettering) annotates somatic alterations and stratifies their treatment implications
"by the level of evidence that a specific molecular alteration is predictive of drug response on the basis of
US Food and Drug Administration labeling, National Comprehensive Cancer Network guidelines, disease-focused
expert group recommendations, and scientific literature" [1].

### 2.2 Core Model

OncoKB defines seven therapeutic levels (verbatim definitions [2]):

| Level | Category | Definition |
|-------|----------|------------|
| 1 | Standard Care | FDA-recognized biomarker predictive of response to an FDA-approved drug in this indication |
| 2 | Standard Care | Standard care biomarker recommended by the NCCN or other professional guidelines predictive of response to an FDA-approved drug in this indication |
| 3A | Investigational | Compelling clinical evidence supports the biomarker as being predictive of response to a drug in this indication |
| 3B | Investigational | Standard care or investigational biomarker predictive of response to an FDA-approved or investigational drug in another indication |
| 4 | Hypothetical | Compelling biological evidence supports the biomarker as being predictive of response to a drug |
| R1 | Standard Care Resistance | Standard care biomarker predictive of resistance to an FDA-approved drug in this indication |
| R2 | Investigational Resistance | Compelling clinical evidence supports the biomarker as being predictive of resistance to a drug |

The *actionability* of a variant is the maximum level over its associations under the combined precedence
order **R1 > 1 > 2 > 3A > 3B > 4 > R2** [4]. Two single-axis maxima are also defined: highest sensitivity
level uses **1 > 2 > 3A > 3B > 4**, and highest resistance level uses **R1 > R2** [4]. Levels 1, 2, R1 are
"standard care"; Levels 3A, 3B, 4, R2 are investigational/hypothetical [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Highest combined level = max of associations under R1 > 1 > 2 > 3A > 3B > 4 > R2 | OncoKB HIGHEST_LEVEL order [4] |
| INV-02 | Highest sensitive level uses 1 > 2 > 3A > 3B > 4, ignoring R1/R2 | OncoKB HIGHEST_SENSITIVE_LEVEL [4] |
| INV-03 | Highest resistance level uses R1 > R2, ignoring sensitivity levels | OncoKB HIGHEST_RESISTANCE_LEVEL [4] |
| INV-04 | Zero leveled associations ⇒ None on all axes (not actionable) | Annotator leaves HIGHEST_LEVEL empty [4]; see 5.4 |
| INV-05 | One assessment per input variant, input order preserved | Implementation contract |

### 2.5 Comparison with Related Methods

| Aspect | OncoKB therapeutic levels (this) | AMP/ASCO/CAP tiers (ONCO-ANNOT-001) |
|--------|----------------------------------|--------------------------------------|
| Axis classified | Therapeutic actionability (drug response/resistance) | Clinical significance (Tier I–IV) |
| Output | Level 1/2/3A/3B/4/R1/R2 or None | Tier I/II/III/IV |
| Source | Chakravarty 2017 [1] | Li 2017 |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| variant | `VariantActionabilityInput` | required | Gene/protein change + caller-supplied leveled drug associations | `Associations` non-null (may be empty) |
| variants | `IEnumerable<VariantActionabilityInput>` | required | Batch of variants for `AssessActionability` | non-null |
| level | `OncoKbLevel` | required | A single level (for `IsStandardCare`) | enum value |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| HighestSensitiveLevel | `OncoKbLevel` | Max sensitivity level (1 > 2 > 3A > 3B > 4) or None |
| HighestResistanceLevel | `OncoKbLevel` | Max resistance level (R1 > R2) or None |
| HighestCombinedLevel | `OncoKbLevel` | Max over both axes (R1 > 1 > 2 > 3A > 3B > 4 > R2) or None |
| IsActionable | `bool` | True iff HighestCombinedLevel ≠ None |

### 3.3 Preconditions and Validation

`AssessActionability(null)` and `ClassifyActionabilityLevel`/`GetTherapyRecommendations` on a variant whose
`Associations` is null throw `ArgumentNullException`. `VariantActionabilityInput` rejects a null
`associations` list at construction. An empty association list is valid and yields `None`/`NotActionable`.
The classifier is case- and content-agnostic about drug names; it ranks only the supplied `OncoKbLevel`
values.

## 4. Algorithm

### 4.1 High-Level Steps

1. For each variant, scan its caller-supplied associations.
2. Track the maximum level under the combined order (for `HighestCombinedLevel`).
3. Track the maximum among sensitivity levels and among resistance levels separately.
4. A variant with no leveled association yields `None` on every axis (not actionable).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The combined precedence order **R1 > 1 > 2 > 3A > 3B > 4 > R2** [4] is encoded directly in the integer order
of the `OncoKbLevel` enum (None lowest … R1 highest), so `CompareLevels` is an integer comparison of enum
values. Sensitivity-set {1,2,3A,3B,4} and resistance-set {R1,R2} are the membership filters for the
single-axis maxima [4]; standard-care set {1,2,R1} is from the SOP grouping [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ClassifyActionabilityLevel` | O(k) | O(1) | k = associations on the variant |
| `AssessActionability` | O(n·k) | O(n) | n variants; matches Registry "O(n × k)" |
| `GetTherapyRecommendations` | O(k log k) | O(k) | sort by level |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.AssessActionability(...)`: batch per-variant assessment (combined/sensitive/resistance).
- `OncologyAnalyzer.ClassifyActionabilityLevel(...)`: single-variant highest combined level.
- `OncologyAnalyzer.GetTherapyRecommendations(...)`: associations ordered most-actionable first.
- `OncologyAnalyzer.CompareLevels(...)`: combined-order comparator.
- `OncologyAnalyzer.IsStandardCare(...)`: standard-care grouping predicate.

### 5.2 Current Behavior

The ranking is a single linear scan per variant. No substring/pattern search is involved, so the repository
suffix tree is **not** applicable (it is for sequence occurrence search, not enum ranking). The knowledgebase
is caller-supplied; the library performs no database lookup.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All seven levels (1, 2, 3A, 3B, 4, R1, R2) with verbatim definitions from OncoKB Levels V2 [2].
- Combined order R1 > 1 > 2 > 3A > 3B > 4 > R2; sensitive order 1 > 2 > 3A > 3B > 4; resistance order R1 > R2 [4].
- Standard-care grouping {1, 2, R1} [3].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Diagnostic (Dx) and prognostic (Px) level systems [4]; **users should rely on:** no current alternative — scope is therapeutic actionability only.
- The curated OncoKB content itself (biomarker→level lookup); **users should rely on:** caller-supplied associations (framework boundary).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | No-association ⇒ `None` (NotActionable) | Assumption | Behavior for a variant with no leveled drug match | accepted | Matches annotator empty HIGHEST_LEVEL [4]; name is ours (A1). |
| 2 | Caller-supplied knowledgebase | Assumption | Library ranks levels, does not curate | accepted | Framework boundary (A2). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No leveled associations | `None` on all axes; `IsActionable=false` | INV-04 [4] |
| Both sensitive and resistance associations | each axis reported; combined = max | separate axes [3][4] |
| {1, R1} | combined = R1 (R1 > 1) | combined order [4] |
| Null `variants` / null `Associations` | `ArgumentNullException` | validation contract |

### 6.2 Limitations

Therapeutic axis only (no Dx/Px). Does not curate biomarkers, evaluate tumor-type matching, or resolve
indication context — it ranks the levels the caller supplies. Conflicting curated data is out of scope; per
the SOP, Levels 1/2/R1 are fixed by guideline inclusion [3].

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var variant = new OncologyAnalyzer.VariantActionabilityInput(
    "BRAF", "p.V600E",
    new[]
    {
        new OncologyAnalyzer.TherapyAssociation("Dabrafenib", OncologyAnalyzer.OncoKbLevel.Level1),
        new OncologyAnalyzer.TherapyAssociation("OtherDrug",  OncologyAnalyzer.OncoKbLevel.Level3A),
    });

var level = OncologyAnalyzer.ClassifyActionabilityLevel(variant); // Level1 (1 > 3A)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_AssessActionability_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AssessActionability_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [ONCO-ACTION-001-Evidence.md](../../../docs/Evidence/ONCO-ACTION-001-Evidence.md)
- Related algorithms: [Cancer_Variant_Annotation](Cancer_Variant_Annotation.md)

## 8. References

1. Chakravarty D, Gao J, Phillips SM, et al. 2017. OncoKB: A Precision Oncology Knowledge Base. JCO Precision Oncology 2017:1-16. https://doi.org/10.1200/PO.17.00011
2. OncoKB. Therapeutic Levels of Evidence (V2). Memorial Sloan Kettering Cancer Center. https://www.oncokb.org/content/files/levelOfEvidence/V2/LevelsOfEvidence.pdf
3. OncoKB. Curation Standard Operating Procedure v3. https://sop.oncokb.org/static/sop/OncoKB_Curation_Standard_Operating_Procedure_v3.pdf
4. oncokb-annotator. README — annotation columns HIGHEST_LEVEL / HIGHEST_SENSITIVE_LEVEL / HIGHEST_RESISTANCE_LEVEL. GitHub. https://github.com/oncokb/oncokb-annotator
