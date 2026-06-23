# Homologous Recombination Deficiency (HRD) Score

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-HRD-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-23 |

## 1. Overview

The HRD score is a composite genomic-instability biomarker used to identify tumours with homologous-recombination repair deficiency (e.g. *BRCA1/2*, *PALB2*), which predicts response to platinum chemotherapy and PARP inhibitors. It is the unweighted sum of three independent genomic-scar metrics: loss of heterozygosity (LOH), telomeric allelic imbalance (TAI), and large-scale state transitions (LST) [1]. A tumour is classified HRD-high when the sum is at or above a fixed clinical cutoff of 42 [1]. This implementation computes the composite sum and the threshold classification, and additionally **derives the HRD-LOH component end-to-end from allele-specific copy-number segments** (Abkevich 2012 / scarHRD `calc.hrd` [2][7]); the TAI and LST components remain caller-supplied because their faithful derivation requires scarHRD's exact binary centromere coordinate table (see §5.3 / §6.2).

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
| INV-05 | `DetectHRD(segments, tai, lst).Components.Loh == DetectLOH(segments).Score` | LOH derived from segments, not supplied [2][7] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| loh | int | required | HRD-LOH component count | ≥ 0 |
| tai | int | required | TAI/NtAI component count | ≥ 0 |
| lst | int | required | LST component count | ≥ 0 |
| score | int | required | combined HRD score for classification | ≥ 0 |
| components | HrdComponents | required | bundle of (Loh, Tai, Lst) for `DetectHRD` | each ≥ 0 |
| segments | IEnumerable&lt;AlleleSpecificSegment&gt; | required | allele-specific CN segments LOH is derived from (`DetectHRD(segments,tai,lst)`) | non-null; End > Start; CN ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| CalculateHRDScore | int | combined score = loh + tai + lst |
| ClassifyHRDStatus | HrdStatus | HrdHigh (≥ 42) or HrdNegative |
| DetectHRD | HrdResult | components, summed score, and status |
| DetectHRD(segments,tai,lst) | HrdResult | LOH derived from segments + caller-supplied TAI/LST, summed score, and status |

### 3.3 Preconditions and Validation

Any negative component count throws `ArgumentOutOfRangeException` (`CalculateHRDScore`, `DetectHRD`). A negative score throws `ArgumentOutOfRangeException` (`ClassifyHRDStatus`). The cutoff comparison is inclusive at 42. The segment-driven `DetectHRD(segments, tai, lst)` throws `ArgumentNullException` for null segments, `ArgumentException` for a segment with non-positive length or negative copy number (via `DetectLOH`), and `ArgumentOutOfRangeException` for a negative `tai`/`lst`.

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
- `OncologyAnalyzer.DetectHRD(HrdComponents)`: end-to-end sum + classification from supplied counts.
- `OncologyAnalyzer.DetectHRD(IEnumerable<AlleleSpecificSegment>, int tai, int lst)`: derives the HRD-LOH count from segments via `DetectLOH` (Abkevich 2012 / scarHRD `calc.hrd`), then sums with caller-supplied TAI/LST and classifies.

### 5.2 Current Behavior

`CalculateHRDScore` / `ClassifyHRDStatus` / `DetectHRD(HrdComponents)` take the three counts as already-computed inputs. The segment overload `DetectHRD(segments, tai, lst)` derives LOH end-to-end from allele-specific segments (delegating to the scarHRD-faithful `DetectLOH`) while TAI/LST stay caller-supplied (see §5.3). The cutoff is exposed as the public constant `HrdHighScoreThreshold = 42`. LOH derivation reuses the per-chromosome grouping + same-state merge already in the LOH path; no substring/pattern search is involved, so the repository suffix tree is not applicable (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- HRD = LOH + TAI + LST as an unweighted sum [1] (scarHRD `scar_score.R`: `sum_HRD0 <- res_lst + res_hrd + res_ai[1]` [7]).
- HRD-high classification at score ≥ 42, boundary inclusive [1][5].
- **HRD-LOH derived from allele-specific segments**: `DetectHRD(segments, tai, lst)` derives LOH via `DetectLOH` — LOH regions with minor==0 & major≠0, length > 15 Mb, whole-chromosome-LOH excluded, after same-state merge (Abkevich 2012 [2]; scarHRD `calc.hrd.R` [7]; oncoscanR `score_loh` [6]).

**Intentionally simplified:**

- TAI and LST are **not** derived from segments; they are taken as caller-supplied counts. **Consequence:** for those two components callers must supply valid counts from an external tool (e.g. scarHRD). See "Not implemented" for why.

**Not implemented:**

- **Per-segment TAI derivation** (scarHRD `calc.ai_new` [7]) and **per-segment LST derivation** (scarHRD `calc.lst` [7]); **users should rely on:** an external scarHRD run for these two counts. Reason: both require scarHRD's exact per-build centromere/telomere `chrominfo` coordinate table — TAI's telomeric-vs-interstitial classification (`chrominfo[i,2]`/`[i,3]`) and LST's p/q-arm split (`q.arm[1,3] <- chrominfo[i,3]`, `p.arm[nrow,4] <- chrominfo[i,2]`) are sensitive to those coordinates. That table ships only as binary `R/sysdata.rda` and could not be retrieved as a verifiable numeric reference in this session; deriving TAI/LST from an unverified centromere table would not reproduce scarHRD, so it is deliberately not attempted (Evidence §scarHRD point 4).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Score exactly 42 | HRD-high | inclusive cutoff [1] |
| Score 41 | HRD-negative | below cutoff [1] |
| All components 0 (near-diploid) | score 0, HRD-negative | low-signal tumour; sum well-defined |
| Negative component or score | ArgumentOutOfRangeException | counts are non-negative |

### 6.2 Limitations

The HRD-LOH component is derived from segments, but TAI and LST are taken as caller-supplied counts; the quality of the final score therefore still depends on the upstream TAI/LST values. Their in-library derivation is blocked on scarHRD's exact centromere `chrominfo` table (binary, not retrievable as a verifiable reference; see §5.3). The 42 cutoff is the myChoice/Telli-2016 value validated primarily in breast and ovarian carcinoma [1][5]; tissue-specific cutoffs are out of scope.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Supplied-components path:
var result = OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(Loh: 20, Tai: 15, Lst: 12));
// result.Score  == 47
// result.Status == OncologyAnalyzer.HrdStatus.HrdHigh   (47 ≥ 42)

// Segment-driven path: LOH derived from allele-specific segments, TAI/LST supplied:
var fromSegments = OncologyAnalyzer.DetectHRD(segments, tai: 25, lst: 16);
// fromSegments.Components.Loh is derived via DetectLOH (scarHRD calc.hrd)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_CalculateHRDScore_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateHRDScore_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ONCO-HRD-001-Evidence.md](../../../docs/Evidence/ONCO-HRD-001-Evidence.md)
- Related algorithm: [Loss of Heterozygosity (HRD-LOH)](Loss_Of_Heterozygosity.md) (ONCO-LOH-001) — the derived LOH component

## 8. References

1. Telli ML, Timms KM, Reid J, et al. 2016. Homologous Recombination Deficiency (HRD) Score Predicts Response to Platinum-Containing Neoadjuvant Chemotherapy in Patients with Triple-Negative Breast Cancer. Clin Cancer Res 22(15):3764–3773. https://pubmed.ncbi.nlm.nih.gov/26957554/
2. Abkevich V, Timms KM, Hennessy BT, et al. 2012. Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. Br J Cancer 107(10):1776–1782. https://www.nature.com/articles/bjc2012451
3. Birkbak NJ, Wang ZC, Kim JY, et al. 2012. Telomeric allelic imbalance indicates defective DNA repair and sensitivity to DNA-damaging agents. Cancer Discov 2(4):366–375. https://pmc.ncbi.nlm.nih.gov/articles/PMC3806629/
4. Popova T, Manié E, Rieunier G, et al. 2012. Ploidy and large-scale genomic instability consistently identify basal-like breast carcinomas with BRCA1/2 inactivation. Cancer Res 72(21):5454–5462. https://aacrjournals.org/cancerres/article/72/21/5454/576090/
5. Stewart MD, Merino Vega D, Arend RC, et al. 2022. Homologous Recombination Deficiency: Concepts, Definitions, and Assays. Oncologist 27(3):167–174. https://pmc.ncbi.nlm.nih.gov/articles/PMC8914493/
6. Christinat Y. oncoscanR `score_loh` documentation. https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html
7. Sztupinszki Z, Diossy M, Krzystanek M, et al. 2018. scarHRD — reference implementation (`R/calc.hrd.R`, `R/calc.ai_new.R`, `R/calc.lst.R`, `R/scar_score.R`). GitHub. https://github.com/sztup/scarHRD
