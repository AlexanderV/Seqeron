# Homozygous (Deep) Deletion Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Copy-Number Alteration |
| Test Unit ID | ONCO-CNA-003 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Homozygous (deep) deletion detection identifies copy-number segments that have lost **both** alleles — total copy number 0 — across a set of arm-anchored copy-number segments, and maps those segments to the recurrently deleted tumour-suppressor genes resident on their chromosome arms. It is a specification-driven, deterministic filter built on the discrete copy-number classification of ONCO-CNA-001: a segment is a homozygous deletion exactly when its hard-threshold integer copy number is 0 (the cBioPortal "−2" Deep Deletion / DeepDeletion state) [1][2][3]. It is used to flag complete loss of tumour-suppressor loci (e.g. CDKN2A, PTEN, TP53), which is a hallmark of tumour-suppressor inactivation [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Somatic copy-number alterations span a spectrum from single-copy (heterozygous) loss to complete (homozygous) loss. A **homozygous deletion** is the loss of all copies at a locus — "zero copies of both alleles in the tumour cells" — and requires two independent deletion events; it is rare and recurrently targets tumour-suppressor genes [1]. Discrete copy-number platforms encode this as the deepest loss state: cBioPortal's "−2" / "Deep Deletion" is "a deep loss, possibly a homozygous deletion", whereas "−1" / "Shallow Deletion" is a single-copy (heterozygous) loss that is **not** homozygous [2][3].

### 2.2 Core Model

Given a segment's mean log2 copy ratio, an integer copy number is called by CNVkit's hard-threshold method (`absolute_threshold`): the copy number is the index of the first ascending cutoff the log2 value is less than or equal to, counting up from 0 [4]. With the default cutoffs (−1.1, −0.25, 0.2, 0.7) [4]:

- log2 ≤ −1.1 ⇒ integer CN **0** ⇒ DeepDeletion ⇒ **homozygous deletion**;
- −1.1 < log2 ≤ −0.25 ⇒ integer CN 1 ⇒ Loss (heterozygous, single-copy) ⇒ not homozygous;
- otherwise CN ≥ 2 ⇒ Neutral / Gain / Amplification ⇒ not homozygous.

A segment is therefore a homozygous deletion **iff** its integer copy number equals 0 [1][2][4]. Tumour-suppressor mapping is by chromosome arm: a homozygous deletion on an arm reports the panel gene(s) on that arm, by NCBI Gene cytogenetic location [5].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A segment is reported iff its integer copy number is exactly 0 | Definition: homozygous = total CN 0 [1]; cBioPortal −2 [2]; CNVkit integer CN [4] |
| INV-02 | A single-copy (heterozygous) loss (integer CN 1) is never reported | cBioPortal −1 ≠ −2 [2][3]; one allele remains [1] |
| INV-03 | The result is a subset of the input in input order (order-preserving filter) | Filter semantics; mirrors `DetectFocalAmplifications` |
| INV-04 | Tumour-suppressor mapping is by arm, in fixed panel order, each gene reported once | NCBI Gene arms [5]; `HashSet` of arms + ordered panel |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| segments / deletions | `IEnumerable<CopyNumberArmSegment>` | required | Arm-anchored copy-number segments (Arm, Start, End, ArmLength, Log2Ratio) | non-null; each segment ArmLength > 0 and End > Start |
| thresholds | `IReadOnlyList<double>?` | null → (−1.1, −0.25, 0.2, 0.7) | Four strictly ascending log2 cutoffs | exactly 4, strictly ascending |
| ploidy | `double` | 2 | Reference (germline) ploidy | > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (DetectHomozygousDeletions) | `IReadOnlyList<CopyNumberArmSegment>` | The CN-0 segments, in input order |
| (IsHomozygousDeletion) | `bool` | true iff the segment's integer CN is 0 |
| (IdentifyDeletedTumorSuppressors) | `IReadOnlyList<string>` | Distinct tumour-suppressor symbols on deleted arms, in panel order |

### 3.3 Preconditions and Validation

Null `segments`/`deletions` → `ArgumentNullException`. A segment with non-positive `ArmLength` or `End ≤ Start` → `ArgumentException`. Non-positive `ploidy` → `ArgumentOutOfRangeException`; `thresholds` not four strictly ascending values → `ArgumentException` (both via `CallCopyNumber`). A NaN log2 ratio is a CNVkit no-call returning the neutral reference copy number (rounded ploidy), so it is **not** reported as a homozygous deletion [4]. Coordinates are bp; arm labels are chromosome number + p/q (e.g. "17p").

## 4. Algorithm

### 4.1 High-Level Steps

1. For each segment, validate it and compute its integer copy number via CNVkit `absolute_threshold` (`CallCopyNumber`).
2. Report the segment iff that integer copy number is 0 (homozygous / DeepDeletion); preserve input order.
3. For gene mapping, collect the distinct arms of the deletion segments and emit each panel tumour suppressor whose arm is present, in panel order.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Tumour-suppressor panel and arms (NCBI Gene cytogenetic locations) [5]:

| Gene | Cytoband | Arm |
|------|----------|-----|
| TP53 | 17p13.1 | 17p |
| RB1 | 13q14.2 | 13q |
| CDKN2A | 9p21.3 | 9p |
| PTEN | 10q23.31 | 10q |
| BRCA1 | 17q21.31 | 17q |
| BRCA2 | 13q13.1 | 13q |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| DetectHomozygousDeletions | O(n) | O(k) | n segments; k = number of CN-0 segments; O(1) per-segment classification |
| IdentifyDeletedTumorSuppressors | O(n + g) | O(a) | a distinct arms; g = panel size (6) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.DetectHomozygousDeletions(segments, thresholds?, ploidy?)`: order-preserving filter of CN-0 segments.
- `OncologyAnalyzer.IsHomozygousDeletion(segment, thresholds?, ploidy?)`: predicate, integer CN == 0.
- `OncologyAnalyzer.IdentifyDeletedTumorSuppressors(deletions)`: arm → tumour-suppressor panel mapping.

### 5.2 Current Behavior

Reuses the ONCO-CNA-002 `CopyNumberArmSegment` record and `ValidateArmSegment`, and the ONCO-CNA-001 `CallCopyNumber` integer-CN calling, so the same segment input drives focal-amplification and homozygous-deletion detection. No new copy-number threshold is introduced; the homozygous-deletion call is defined solely as integer CN 0. **Search reuse:** the unit is a numeric filter / dictionary lookup, not substring/pattern matching, so the repository suffix tree is not applicable (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Homozygous deletion = total copy number 0 (both alleles lost) [1]; called as integer CN 0 via CNVkit `absolute_threshold` [4] and cBioPortal Deep Deletion ("−2") [2][3].
- Tumour-suppressor arms from NCBI Gene cytogenetic locations [5].

**Intentionally simplified:**

- Homozygous status is inferred from total integer copy number (CN 0), not from allele-specific copy number; **consequence:** a copy-neutral LOH or an allele-specific zero with a retained other allele is not distinguished here (total-CN model, consistent with cBioPortal discrete calls).

**Not implemented:**

- Purity/ploidy correction of the discrete calls beyond the `ploidy` parameter; **users should rely on:** purity/ploidy estimation units (ONCO-PURITY/PLOIDY) upstream before classification.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Single-copy loss (CN 1, e.g. log2 = −0.5) | Not reported | cBioPortal −1 (heterozygous) ≠ homozygous [2][3] |
| log2 exactly −1.1 | Reported (CN 0) | CNVkit "≤ each threshold in sequence" [4] |
| log2 just above −1.1 | Not reported (CN 1) | CNVkit threshold boundary [4] |
| NaN log2 | Not reported (neutral no-call) | CNVkit no-call → reference CN [4] |
| Empty input | Empty result | Filter of empty set |
| Null input | `ArgumentNullException` | Validation |
| Deletion on non-panel arm | No gene reported | Closed panel |

### 6.2 Limitations

Uses total copy number, not allele-specific copy number; cannot separate homozygous deletion from copy-neutral LOH. Discrete calls are putative and sensitive to tumour purity/ploidy [2]; the gene panel is the fixed six-gene tumour-suppressor list (TP53, RB1, CDKN2A, PTEN, BRCA1, BRCA2) and does not annotate other deleted loci.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
var segs = new[]
{
    new OncologyAnalyzer.CopyNumberArmSegment("9p", 0, 1_000, 40_000_000, -2.0), // CN 0
    new OncologyAnalyzer.CopyNumberArmSegment("3p", 0, 1_000, 90_000_000, -0.5), // CN 1
};
var hom = OncologyAnalyzer.DetectHomozygousDeletions(segs);          // -> [9p segment]
var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(hom);   // -> ["CDKN2A"]
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_DetectHomozygousDeletions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_DetectHomozygousDeletions_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-CNA-003-Evidence.md](../../../docs/Evidence/ONCO-CNA-003-Evidence.md)
- Related algorithms: [Cancer Copy Number Alteration Classification](./Copy_Number_Alteration_Classification.md), [Focal Amplification Detection](./Focal_Amplification_Detection.md)

## 8. References

1. Cheng J, Demeulemeester J, Wedge DC, et al. 2017. Pan-cancer analysis of homozygous deletions in primary tumours uncovers rare tumour suppressors. Nature Communications 8:1221. https://pmc.ncbi.nlm.nih.gov/articles/PMC5663922/
2. cBioPortal. Discrete Copy Number data file format. https://docs.cbioportal.org/file-formats/ (accessed 2026-06-14)
3. cBioPortal. FAQ — meaning of Amplification / Gain / Deep Deletion / Shallow Deletion / −2..2. https://docs.cbioportal.org/user-guide/faq/ (accessed 2026-06-14)
4. Talevich E, Shain AH, Botton T, Bastian BC. 2016. CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol 12(4):e1004873. `cnvlib/call.py` `absolute_threshold`. https://cnvkit.readthedocs.io/
5. NCBI Gene cytogenetic locations (accessed 2026-06-14): TP53 https://www.ncbi.nlm.nih.gov/gene/7157 ; RB1 https://www.ncbi.nlm.nih.gov/gene/5925 ; CDKN2A https://www.ncbi.nlm.nih.gov/gene/1029 ; PTEN https://www.ncbi.nlm.nih.gov/gene/5728 ; BRCA1 https://www.ncbi.nlm.nih.gov/gene/672 ; BRCA2 https://www.ncbi.nlm.nih.gov/gene/675
