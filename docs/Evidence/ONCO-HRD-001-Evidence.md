# Evidence Artifact: ONCO-HRD-001

**Test Unit ID:** ONCO-HRD-001
**Algorithm:** Homologous Recombination Deficiency (HRD) composite genomic-scar score
**Date Collected:** 2026-06-14

---

## Online Sources

### Telli et al. (2016), Clin Cancer Res — combined HRD score definition and cutoff

**URL:** https://pubmed.ncbi.nlm.nih.gov/26957554/
**Accessed:** 2026-06-14 (retrieved via WebFetch on the PubMed abstract page)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Composite definition (verbatim):** "combined homologous recombination deficiency (HRD) score, an unweighted sum of LOH, TAI, and LST scores". → The composite HRD score = LOH + TAI + LST, with no weighting.
2. **HRD-high cutoff (verbatim):** "HR deficiency, defined as HRD score ≥42 or BRCA1/2 mutation". → The genomic-scar cutoff for HRD positivity is ≥ 42 (boundary inclusive).

### Stewart et al. (2022), "Homologous Recombination Deficiency: Concepts, Definitions, and Assays" (review) — component definitions

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8914493/
**Accessed:** 2026-06-14 (retrieved via WebFetch)
**Authority rank:** 1 (peer-reviewed review summarising primary sources)

**Key Extracted Points:**

1. **gLOH/HRD-LOH (verbatim):** "regions of intermediate size (>15 MB and < whole chromosome)" — attributed to Abkevich. → LOH score counts LOH regions longer than 15 Mb but shorter than the whole chromosome.
2. **TAI (verbatim):** "the number of regions with allelic imbalance which extend to the sub-telomere but not cross the centromere" — attributed to Birkbak.
3. **LST (verbatim):** "chromosome breaks (translocations, inversions, or deletions)" — attributed to Popova.
4. **Combination (verbatim):** myChoice CDx "determines HR status by … genomic instability—measured by the evaluation of a combination of molecular measures to derive a genomic instability (gLOH + TAI + LST)". → confirms the three-way sum.
5. **Cutoff (verbatim):** trials used "a cutoff of 42" to define HR status. → confirms 42 independently of Telli 2016.

### Birkbak et al. (2012), Cancer Discovery — TAI / NtAI primary definition

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3806629/
**Accessed:** 2026-06-14 (retrieved via WebFetch on PMC full text)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Allelic imbalance (verbatim):** "Allelic imbalance was defined as any time the copy number of the two alleles were not equal, and at least one allele was present."
2. **Telomeric requirement (verbatim):** "Telomeric AI and telomeric CNA are defined as regions that extend to one of the sub-telomeres but do not cross the centromere."

### oncoscanR (Christinat) `score_loh` reference implementation — HRD-LOH operational definition

**URL:** https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html
**Accessed:** 2026-06-14 (retrieved via WebFetch)
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **LOH counting rule (verbatim):** "All LOH segments larger than 15Mb but excluding chromosome with a global LOH alteration" — based on Abkevich et al., Br J Cancer 2012. → confirms 15 Mb minimum and whole-chromosome exclusion.

### Birkbak/Popova secondary summary — npj Precision Oncology (2022)

**URL:** https://www.nature.com/articles/s41698-022-00339-8 (search snippet via WebSearch)
**Accessed:** 2026-06-14 (retrieved via WebSearch; page itself behind an IDP redirect)
**Authority rank:** 1 (peer-reviewed)

**Key Extracted Points:**

1. **Sum + cutoff (verbatim snippet):** "The HRD score is an unweighted sum of LOH, TAI, and LST scores" and the Telli2016 selection criterion is "HRD-AIs ≥ 42 considered HRD". → corroborates Telli 2016.

### Popova et al. (2012), Cancer Research — LST primary definition

**URL:** https://aacrjournals.org/cancerres/article/72/21/5454/576090/ (snippet via WebSearch; rdrr.io calc.lst doc)
**Accessed:** 2026-06-14 (retrieved via WebSearch)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **LST (verbatim snippet):** "LST was defined as a chromosomal breakpoint (change in copy number or allelic content) between adjacent regions each of at least 10 megabases (Mb) obtained after smoothing and filtering <3 Mb small-scale copy number variation."

---

## Documented Corner Cases and Failure Modes

### From Telli et al. (2016)

1. **Boundary at 42:** the cutoff is inclusive (score ≥ 42 is HRD-high); a score of exactly 42 is HRD-high and 41 is HRD-negative.

### From the per-unit NOTE / Stewart et al. (2022)

1. **Near-diploid / low-signal tumours:** small component counts produce a low sum that is HRD-negative (e.g. 0 + 0 + 0 = 0). The arithmetic is well-defined for any non-negative counts.
2. **Non-negative counts:** the three components are event counts (LOH regions / telomeric AIs / LSTs), so each is a non-negative integer; a negative input is invalid.

---

## Test Datasets

### Dataset: Telli 2016 cutoff boundary

**Source:** Telli ML et al. (2016), Clin Cancer Res 22(15):3764–3773 — "HRD score ≥42".

| Parameter | Value |
|-----------|-------|
| HRD-high cutoff | 42 (inclusive) |
| Composite rule | HRD = LOH + TAI + LST (unweighted) |

### Dataset: worked boundary examples (derived arithmetically from the cited sum + cutoff)

**Source:** Sum definition (Telli 2016) applied to component triples; not from any published table.

| LOH | TAI | LST | Score | Status |
|-----|-----|-----|-------|--------|
| 14 | 14 | 14 | 42 | HRD-high (boundary) |
| 14 | 13 | 14 | 41 | HRD-negative |
| 20 | 15 | 12 | 47 | HRD-high |
| 5 | 4 | 3 | 12 | HRD-negative |
| 0 | 0 | 0 | 0 | HRD-negative (near-diploid) |

---

## Assumptions

1. **ASSUMPTION: component-count input shape** — This unit implements only the retrievable composite-sum + threshold classification, taking the three component counts (LOH, TAI, LST) as already-computed integer inputs (per the unit NOTE). Computing the components from raw segmented copy-number/allelic data (the 15 Mb LOH segmentation, NtAI sub-telomere/centromere geometry, Popova LST smoothing) is out of scope here and is left to ONCO-LOH-001 / ONCO-CNA-001. This is an API-shape decision, not a correctness-affecting parameter: the sum and the 42 cutoff are fully source-backed.

---

## Recommendations for Test Coverage

1. **MUST Test:** composite score = LOH + TAI + LST for representative triples — Evidence: Telli 2016 ("unweighted sum").
2. **MUST Test:** classification at the boundary — score 42 → HRD-high, 41 → HRD-negative — Evidence: Telli 2016 ("≥42").
3. **MUST Test:** end-to-end `DetectHRD` returns the summed score and matching status — Evidence: Telli 2016.
4. **SHOULD Test:** near-diploid all-zero input → score 0, HRD-negative — Rationale: documented low-signal edge case.
5. **SHOULD Test:** negative component / negative score rejected with `ArgumentOutOfRangeException` — Rationale: counts are non-negative.
6. **COULD Test:** ordering invariance of the sum (LOH+TAI+LST is symmetric) — Rationale: confirms unweighted/commutative sum.

---

## References

1. Telli ML, Timms KM, Reid J, et al. (2016). Homologous Recombination Deficiency (HRD) Score Predicts Response to Platinum-Containing Neoadjuvant Chemotherapy in Patients with Triple-Negative Breast Cancer. Clin Cancer Res 22(15):3764–3773. https://pubmed.ncbi.nlm.nih.gov/26957554/
2. Abkevich V, Timms KM, Hennessy BT, et al. (2012). Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. Br J Cancer 107(10):1776–1782. https://www.nature.com/articles/bjc2012451
3. Birkbak NJ, Wang ZC, Kim JY, et al. (2012). Telomeric allelic imbalance indicates defective DNA repair and sensitivity to DNA-damaging agents. Cancer Discov 2(4):366–375. https://pmc.ncbi.nlm.nih.gov/articles/PMC3806629/
4. Popova T, Manié E, Rieunier G, et al. (2012). Ploidy and large-scale genomic instability consistently identify basal-like breast carcinomas with BRCA1/2 inactivation. Cancer Res 72(21):5454–5462. https://aacrjournals.org/cancerres/article/72/21/5454/576090/
5. Stewart MD, Merino Vega D, Arend RC, et al. (2022). Homologous Recombination Deficiency: Concepts, Definitions, and Assays. Oncologist 27(3):167–174. https://pmc.ncbi.nlm.nih.gov/articles/PMC8914493/
6. Christinat Y. oncoscanR `score_loh` documentation. https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html

---

## Change History

- **2026-06-14**: Initial documentation.
