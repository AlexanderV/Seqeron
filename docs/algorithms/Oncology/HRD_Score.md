# Homologous Recombination Deficiency (HRD) Score

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-HRD-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

The HRD score is a composite genomic-instability biomarker used to identify tumours with homologous-recombination repair deficiency (e.g. *BRCA1/2*, *PALB2*), which predicts response to platinum chemotherapy and PARP inhibitors. It is the unweighted sum of three independent genomic-scar metrics: loss of heterozygosity (LOH), telomeric allelic imbalance (TAI), and large-scale state transitions (LST) [1]. A tumour is classified HRD-high when the sum is at or above a fixed clinical cutoff of 42 [1]. This implementation is specification-driven and exact: it computes the composite sum and the threshold classification from the three already-computed component counts.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Homologous-recombination-deficient tumours accumulate characteristic genome-wide "scars". Three independently developed metrics quantify these scars: HRD-LOH (Abkevich 2012) [2], NtAI/TAI (Birkbak 2012) [3], and LST (Popova 2012) [4]. Telli et al. (2016) combined them into a single score and established the clinical HRD-high cutoff used by the Myriad myChoice CDx assay [1][5].

### 2.2 Core Model

The combined HRD score is the **unweighted sum** of the three component counts [1]:

```
HRD = LOH + TAI + LST
```

where each component is a non-negative integer count of genomic events:

- **LOH** — number of LOH regions longer than 15 Mb but shorter than a whole chromosome [2][6].
- **TAI** — number of allelic-imbalance regions that extend to a sub-telomere but do not cross the centromere; allelic imbalance is when "the copy number of the two alleles were not equal, and at least one allele was present" [3].
- **LST** — number of chromosomal breaks between adjacent regions each ≥ 10 Mb, after filtering regions < 3 Mb [4].

Classification [1]:

```
status = HRD-high   if HRD ≥ 42
         HRD-negative otherwise
```

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | HRD = LOH + TAI + LST | unweighted sum [1] |
| INV-02 | Sum is order-independent (commutative) | integer addition; "unweighted" [1] |
| INV-03 | status = HRD-high iff score ≥ 42 (inclusive) | Telli 2016 cutoff "≥42" [1] |
| INV-04 | score ≥ 0 and each component ≥ 0 | components are event counts [2][3][4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| loh | int | required | HRD-LOH component count | ≥ 0 |
| tai | int | required | TAI/NtAI component count | ≥ 0 |
| lst | int | required | LST component count | ≥ 0 |
| score | int | required | combined HRD score for classification | ≥ 0 |
| components | HrdComponents | required | bundle of (Loh, Tai, Lst) for `DetectHRD` | each ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| CalculateHRDScore | int | combined score = loh + tai + lst |
| ClassifyHRDStatus | HrdStatus | HrdHigh (≥ 42) or HrdNegative |
| DetectHRD | HrdResult | components, summed score, and status |

### 3.3 Preconditions and Validation

Any negative component count throws `ArgumentOutOfRangeException` (`CalculateHRDScore`, `DetectHRD`). A negative score throws `ArgumentOutOfRangeException` (`ClassifyHRDStatus`). The cutoff comparison is inclusive at 42.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate that each component count (and the score, for classification) is non-negative.
2. Sum the three components: `score = loh + tai + lst`.
3. Classify: `HrdHigh` if `score ≥ 42`, else `HrdNegative`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Constant | Value | Source |
|----------|-------|--------|
| HRD-high cutoff | 42 (inclusive) | Telli 2016 [1]; Stewart 2022 [5] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateHRDScore / ClassifyHRDStatus / DetectHRD | O(1) | O(1) | three integer additions + one comparison |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CalculateHRDScore(int, int, int)`: returns the unweighted sum.
- `OncologyAnalyzer.ClassifyHRDStatus(int)`: returns `HrdStatus` against the 42 cutoff.
- `OncologyAnalyzer.DetectHRD(HrdComponents)`: end-to-end sum + classification.

### 5.2 Current Behavior

Takes the three component counts as already-computed inputs. The cutoff is exposed as the public constant `HrdHighScoreThreshold = 42`. No search/matching is involved, so the repository suffix tree is not applicable (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- HRD = LOH + TAI + LST as an unweighted sum [1].
- HRD-high classification at score ≥ 42, boundary inclusive [1][5].

**Intentionally simplified:**

- Component computation from raw segmented copy-number/allelic data (15 Mb LOH segmentation [2][6], sub-telomere/centromere TAI geometry [3], Popova LST smoothing/filtering [4]) is **not** performed here; the three counts are taken as inputs. **Consequence:** callers must supply valid component counts (produced by ONCO-LOH-001 / ONCO-CNA-001 once implemented).

**Not implemented:**

- Per-segment LOH/TAI/LST derivation; **users should rely on:** the dependent units ONCO-LOH-001 and ONCO-CNA-001.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Score exactly 42 | HRD-high | inclusive cutoff [1] |
| Score 41 | HRD-negative | below cutoff [1] |
| All components 0 (near-diploid) | score 0, HRD-negative | low-signal tumour; sum well-defined |
| Negative component or score | ArgumentOutOfRangeException | counts are non-negative |

### 6.2 Limitations

Quality of the result depends entirely on the upstream component counts; this unit does not validate their biological derivation. The 42 cutoff is the myChoice/Telli-2016 value validated primarily in breast and ovarian carcinoma [1][5]; tissue-specific cutoffs are out of scope.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var result = OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(Loh: 20, Tai: 15, Lst: 12));
// result.Score  == 47
// result.Status == OncologyAnalyzer.HrdStatus.HrdHigh   (47 ≥ 42)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_CalculateHRDScore_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateHRDScore_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-HRD-001-Evidence.md](../../../docs/Evidence/ONCO-HRD-001-Evidence.md)

## 8. References

1. Telli ML, Timms KM, Reid J, et al. 2016. Homologous Recombination Deficiency (HRD) Score Predicts Response to Platinum-Containing Neoadjuvant Chemotherapy in Patients with Triple-Negative Breast Cancer. Clin Cancer Res 22(15):3764–3773. https://pubmed.ncbi.nlm.nih.gov/26957554/
2. Abkevich V, Timms KM, Hennessy BT, et al. 2012. Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. Br J Cancer 107(10):1776–1782. https://www.nature.com/articles/bjc2012451
3. Birkbak NJ, Wang ZC, Kim JY, et al. 2012. Telomeric allelic imbalance indicates defective DNA repair and sensitivity to DNA-damaging agents. Cancer Discov 2(4):366–375. https://pmc.ncbi.nlm.nih.gov/articles/PMC3806629/
4. Popova T, Manié E, Rieunier G, et al. 2012. Ploidy and large-scale genomic instability consistently identify basal-like breast carcinomas with BRCA1/2 inactivation. Cancer Res 72(21):5454–5462. https://aacrjournals.org/cancerres/article/72/21/5454/576090/
5. Stewart MD, Merino Vega D, Arend RC, et al. 2022. Homologous Recombination Deficiency: Concepts, Definitions, and Assays. Oncologist 27(3):167–174. https://pmc.ncbi.nlm.nih.gov/articles/PMC8914493/
6. Christinat Y. oncoscanR `score_loh` documentation. https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html
