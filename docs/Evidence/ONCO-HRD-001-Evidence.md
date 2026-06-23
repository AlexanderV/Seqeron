# Evidence Artifact: ONCO-HRD-001

**Test Unit ID:** ONCO-HRD-001
**Algorithm:** Homologous Recombination Deficiency (HRD) composite genomic-scar score
**Date Collected:** 2026-06-14 (updated 2026-06-23 for segment-driven LOH derivation)

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

### scarHRD reference implementation (Sztupinszki et al. 2018) — derivation feasibility of LOH/TAI/LST

**URLs / retrieval (2026-06-23):**

- `R/calc.hrd.R` — https://raw.githubusercontent.com/sztup/scarHRD/master/R/calc.hrd.R (WebFetch)
- `R/calc.ai_new.R` (TAI) — https://raw.githubusercontent.com/sztup/scarHRD/master/R/calc.ai_new.R (WebFetch)
- `R/calc.lst.R` (LST) — https://raw.githubusercontent.com/sztup/scarHRD/master/R/calc.lst.R (WebFetch)
- `R/scar_score.R` (build switch + sum) — https://raw.githubusercontent.com/sztup/scarHRD/master/R/scar_score.R (WebFetch)
- `R/shrink.seg.ai.R` (smoothing) — https://raw.githubusercontent.com/sztup/scarHRD/master/R/shrink.seg.ai.R (WebFetch)
- repo tree (file list incl. binary `R/sysdata.rda`, 706 bytes) — https://api.github.com/repos/sztup/scarHRD/git/trees/master?recursive=1 (WebFetch)

**Authority rank:** 3 (reference implementation cross-checking the primary papers).

**Key Extracted Points:**

1. **HRD-LOH (`calc.hrd.R`, verbatim):** whole-chromosome exclusion `if(all(segSamp[segSamp[,2]==j,nB]==0))` → `chrDel`; LOH `segLOH <- segSamp[segSamp[,nB]==0 & segSamp[,nA]!=0,]`; size filter `segLOH[,4]-segLOH[,3] > sizelimit1`; whole-chr removal `!segLOH[,2] %in% chrDel`; preceded by `shrink.seg.ai.wrapper`. → reproducible from the repo's `AlleleSpecificSegment` (Major/Minor CN) with NO centromere table needed. This is exactly what `DetectLOH` implements.
2. **TAI (`calc.ai_new.R`) needs the centromere table:** signature `calc.ai_new(seg, chrominfo, min.size=1e6, cont=0, ploidyByChromosome=TRUE, shrink=TRUE)`. Telomeric-vs-interstitial classification is decided by the centromere coordinates `chrominfo[i,2]`/`chrominfo[i,3]` (e.g. `if(sample.chrom.seg[1,'AI']==2 & nrow!=1 & sample.chrom.seg[1,4] < chrominfo[i,2])` downgrades the first segment 2→1). It also requires per-chromosome ploidy (`ploidyByChromosome`) and an aberrant-cell-fraction column — modelling state absent from `AlleleSpecificSegment`.
3. **LST (`calc.lst.R`) needs the centromere table:** signature `calc.lst(seg, chrominfo, nA=7, chr.arm='no')`. With `chr.arm='no'` the centromere positions split each chromosome into p/q arms and set arm boundaries `q.arm[1,3] <- chrominfo[i,3]`, `p.arm[nrow(p.arm),4] <- chrominfo[i,2]`. The 3 Mb iterative shrink (`while(length(n.3mb)>0){...shrink.seg.ai...}`), the 10 Mb flag `(p.arm[,4]-p.arm[,3]) >= 10e6`, and the breakpoint rule `if(q.arm[k,9]==1 & q.arm[(k-1),9]==1 & (q.arm[k,3]-q.arm[(k-1),4]) < 3e6)` are algorithmically reproducible, but the arm split itself depends on the exact centromere coordinates.
4. **The centromere `chrominfo` table is binary, unretrievable as citable text:** `scar_score.R` selects `chrominfo = chrominfo_grch38 / chrominfo_grch37 / chrominfo_mouse`; those objects live only in `R/sysdata.rda` (706-byte binary). They could NOT be retrieved as a verifiable numeric table in this session. Public centromere tables (CNAqc `chr_coordinates_GRCh38`, rCGH `hg38`) render truncated on their doc pages and come from a different lineage, so they cannot be confirmed to match scarHRD's embedded coordinates. Since TAI's telomeric classification and LST's arm split are sensitive to these exact coordinates, deriving TAI/LST from an unverified centromere table would not reproduce scarHRD within tolerance.
5. **Combination (`scar_score.R`, verbatim):** `sum_HRD0 <- res_lst + res_hrd + res_ai[1]` — confirms the unweighted three-way sum (`res_hrd`=LOH, `res_ai[1]`=TAI, `res_lst`=LST).

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

1. **ASSUMPTION: TAI / LST remain caller-supplied** — The HRD-LOH component is now DERIVED end-to-end from allele-specific segments (`DetectHRD(segments, tai, lst)` → `DetectLOH`), faithful to Abkevich 2012 / scarHRD `calc.hrd` (no centromere table needed). TAI (`calc.ai_new`) and LST (`calc.lst`) are NOT derived: both require scarHRD's exact per-build centromere/telomere `chrominfo` coordinate table, which ships only as binary `R/sysdata.rda` and could not be retrieved as a verifiable numeric table in this session (point 4 under the scarHRD source above). Public centromere tables come from a different lineage and could not be cross-confirmed to match scarHRD, and TAI's telomeric classification + LST's p/q-arm split are sensitive to those exact coordinates. Per the conditional guard, TAI/LST are therefore left caller-supplied rather than approximated from an unverified table — this is a documented blocker on the centromere coordinates, not an invented constant. The sum, the 42 cutoff, and the LOH derivation are fully source-backed.

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
