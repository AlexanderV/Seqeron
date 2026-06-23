# Evidence Artifact: ONCO-HRD-001

**Test Unit ID:** ONCO-HRD-001
**Algorithm:** Homologous Recombination Deficiency (HRD) composite genomic-scar score
**Date Collected:** 2026-06-14 (updated 2026-06-23 for segment-driven LOH derivation; updated 2026-06-23 again for segment-driven TAI + LST derivation, ONCO-HRD-001 fix)

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
4. **The centromere coordinates are PUBLIC and citable (RESOLVED 2026-06-23):** `scar_score.R` selects `chrominfo = chrominfo_grch38 / chrominfo_grch37`; `rownames(chrominfo) = chrominfo$chr` keys it by chromosome, and the calc functions read `chrominfo[i,2]` (centromere start = p-arm boundary) and `chrominfo[i,3]` (centromere end = q-arm boundary). These per-chromosome centromere boundaries are the cytoBand `acen` regions, retrievable as citable text from the UCSC Genome Browser API (see the "UCSC … cytoBand" sources below) and cross-verified against the NCBI GRC modeled-centromere table. The exact-coordinate sensitivity remains, but the coordinates are no longer unretrievable — they are embedded here as a published per-chromosome table for GRCh38 and GRCh37, exactly like the existing WGD/ploidy chromosome-size table. TAI and LST are therefore now derivable.
5. **Combination (`scar_score.R`, verbatim):** `sum_HRD0 <- res_lst + res_hrd + res_ai[1]` — confirms the unweighted three-way sum (`res_hrd`=LOH, `res_ai[1]`=TAI, `res_lst`=LST).

### scarHRD `calc.ai_new.R` (TAI) — verbatim algorithm retrieved 2026-06-23

**URL:** https://raw.githubusercontent.com/sztup/scarHRD/master/R/calc.ai_new.R
**Accessed:** 2026-06-23 (WebFetch, full source)
**Authority rank:** 3

**Key Extracted Points:**

1. **Pre-filter:** `seg <- seg[seg[,4]-seg[,3] >= min.size,]` — drop segments shorter than `min.size` (default **1e6 = 1 Mb**) before AI assignment.
2. **AI state (even/diploid path):** `sample.chrom.seg[,'AI'] <- c(0,2)[match(seg[,7]==seg[,8], c('TRUE','FALSE'))]` — AI = 0 when A-allele CN (col 7 = major) equals B-allele CN (col 8 = minor); AI = 2 (allelic imbalance) otherwise.
3. **Telomeric downgrade (TAI definition):** first segment with AI==2, `nrow!=1`, and `sample.chrom.seg[1,4] < chrominfo[i,2]` (its END before the centromere start) → AI set to 1 (p-telomeric). Last segment with AI==2, `nrow!=1`, and `sample.chrom.seg[nrow,3] > chrominfo[i,3]` (its START after the centromere end) → AI set to 1 (q-telomeric).
4. **Whole-chr:** `if(nrow==1 & AI[1]!=0) AI[1] <- 3` — single-segment imbalanced chromosome is whole-chromosome AI, not telomeric.
5. **TAI count:** `no.events[j,1] <- nrow(seg[seg[,'AI']==1,])` — TAI = number of AI==1 segments; `scar_score` reports this as "Telomeric AI" (`res_ai[1]`).

### scarHRD `calc.lst.R` (LST) — verbatim algorithm retrieved 2026-06-23

**URL:** https://raw.githubusercontent.com/sztup/scarHRD/master/R/calc.lst.R
**Accessed:** 2026-06-23 (WebFetch, full source)
**Authority rank:** 3

**Key Extracted Points:**

1. **Autosomes only:** `chroms <- chroms[!chroms %in% c(23,24,'chr23','chr24','chrX','chrx','chrY','chry')]`.
2. **Arm split (chr.arm='no'):** `p.arm <- seg[seg[,3] <= chrominfo[i,2],]` (start ≤ centromere start); `q.arm <- seg[seg[,4] >= chrominfo[i,3],]` (end ≥ centromere end). Each arm is shrink-merged then clamped to the centromere (`p.arm[last,4] <- chrominfo[i,2]`, `q.arm[1,3] <- chrominfo[i,3]`). Chromosomes with `< 2` segments are skipped.
3. **3 Mb smoothing:** iteratively remove arm segments with length `< 3e6` (re-merging after each removal) until none remain.
4. **Break rule:** flag each arm segment `1` if `length >= 10e6` else `0`; for adjacent pair (k-1,k) count one LST when both flags == 1 AND `(seg[k,3]-seg[k-1,4]) < 3e6`. Sample LST = sum of counted breaks.

### UCSC Genome Browser — hg38 cytoBand acen (GRCh38 centromere coordinates)

**URL:** https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cytoBand;chrom=chrN (per chromosome chr1..chr22)
**Accessed:** 2026-06-23 (WebFetch, per-chromosome to avoid summarizer truncation)
**Authority rank:** 5 (authoritative coordinate database)

**Key Extracted Points:**

1. Centromere [start,end] = [p11 acen `chromStart`, q11 acen `chromEnd`]. GRCh38 per-chromosome (bp), chrominfo col2=start, col3=end:
   chr1 121700000–125100000; chr2 91800000–96000000; chr3 87800000–94000000; chr4 48200000–51800000; chr5 46100000–51400000; chr6 58500000–62600000; chr7 58100000–62100000; chr8 43200000–47200000; chr9 42200000–45500000; chr10 38000000–41600000; chr11 51000000–55800000; chr12 33200000–37800000; chr13 16500000–18900000; chr14 16100000–18200000; chr15 17500000–20500000; chr16 35300000–38400000; chr17 22700000–27400000; chr18 15400000–21500000; chr19 24200000–28100000; chr20 25700000–30400000; chr21 10900000–13000000; chr22 13700000–17400000.

### UCSC Genome Browser — hg19 cytoBand acen (GRCh37 centromere coordinates)

**URL:** https://api.genome.ucsc.edu/getData/track?genome=hg19;track=cytoBand;chrom=chrN (per chromosome chr1..chr22)
**Accessed:** 2026-06-23
**Authority rank:** 5

**Key Extracted Points:**

1. GRCh37 per-chromosome centromere [start,end] (bp):
   chr1 121500000–128900000; chr2 90500000–96800000; chr3 87900000–93900000; chr4 48200000–52700000; chr5 46100000–50700000; chr6 58700000–63300000; chr7 58000000–61700000; chr8 43100000–48100000; chr9 47300000–50700000; chr10 38000000–42300000; chr11 51600000–55700000; chr12 33300000–38200000; chr13 16300000–19500000; chr14 16100000–19100000; chr15 15800000–20700000; chr16 34600000–38600000; chr17 22200000–25800000; chr18 15400000–19000000; chr19 24400000–28600000; chr20 25600000–29400000; chr21 10900000–14300000; chr22 12200000–17900000.

### NCBI GRC modeled centromeres (cross-verification of GRCh38)

**URL:** https://www.ncbi.nlm.nih.gov/grc/human
**Accessed:** 2026-06-23 (WebFetch)
**Authority rank:** 2 (assembly standard owner)

**Key Extracted Points:**

1. GRCh38 modeled CEN1 = 122,026,460–125,184,587; CEN21 = 10,864,561–12,915,808 — agree to cytoband resolution with the UCSC acen bounds embedded here (chr1 121,700,000–125,100,000; chr21 10,900,000–13,000,000), confirming the cytoBand acen coordinates are the correct centromere region.

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

1. **ASSUMPTION (resolved): TAI / LST are now DERIVED from segments.** The prior blocker (centromere `chrominfo` table unretrievable) is resolved: the per-chromosome centromere boundaries are the UCSC cytoBand `acen` regions, retrieved as citable text and cross-verified against the NCBI GRC modeled-centromere table, and embedded here as a published GRCh38/GRCh37 table (like the existing WGD chromosome-size table). LOH, TAI and LST are all derived from `AlleleSpecificSegment` in `DetectHRD(segments)`.
2. **ASSUMPTION: even-ploidy / standard allelic-imbalance path for TAI** — `AlleleSpecificSegment` carries major/minor allele CN but not scarHRD's ASCAT per-sample ploidy / aberrant-cell-fraction columns (`seg[,9]`,`seg[,10]`). The implemented AI rule is scarHRD's default even/diploid path: AI present ⟺ major ≠ minor (the literal `seg[,7]==seg[,8]` test). The odd-ploidy ploidy-normalized branch is not reproduced (it needs the absent ploidy column). This is the dominant path and matches Birkbak's "regions of allelic imbalance"; documented as an intentional simplification in §5.3 of the algorithm doc.

---

## Recommendations for Test Coverage

1. **MUST Test:** composite score = LOH + TAI + LST for representative triples — Evidence: Telli 2016 ("unweighted sum").
2. **MUST Test:** classification at the boundary — score 42 → HRD-high, 41 → HRD-negative — Evidence: Telli 2016 ("≥42").
3. **MUST Test:** end-to-end `DetectHRD` returns the summed score and matching status — Evidence: Telli 2016.
4. **SHOULD Test:** near-diploid all-zero input → score 0, HRD-negative — Rationale: documented low-signal edge case.
5. **SHOULD Test:** negative component / negative score rejected with `ArgumentOutOfRangeException` — Rationale: counts are non-negative.
6. **COULD Test:** ordering invariance of the sum (LOH+TAI+LST is symmetric) — Rationale: confirms unweighted/commutative sum.
7. **MUST Test (TAI):** first-segment imbalance ending before centromere start (p-telomeric) and last-segment imbalance starting after centromere end (q-telomeric) → counted; interstitial imbalance and single-segment whole-chr imbalance → NOT counted; sub-1 Mb segment dropped. — Evidence: `calc.ai_new`.
8. **MUST Test (TAI):** first imbalanced segment whose end is ≥ centromere start (crosses the centromere) → NOT telomeric. — Evidence: `calc.ai_new`.
9. **MUST Test (LST):** two adjacent ≥10 Mb segments separated by <3 Mb on the same arm → 1 break; one neighbour <10 Mb → 0; gap ≥3 Mb → 0; iterative 3 Mb smoothing exposes a break → counted; <2 segments skipped; sex chromosomes ignored. — Evidence: `calc.lst`.
10. **MUST Test (coords):** embedded GRCh38/GRCh37 centromere table matches UCSC cytoBand acen values. — Evidence: UCSC cytoBand.
11. **SHOULD Test (E2E):** `DetectHRD(segments)` derives LOH+TAI+LST and classifies at the 42 cutoff. — Evidence: scar_score.R sum + Telli 2016.

---

## References

1. Telli ML, Timms KM, Reid J, et al. (2016). Homologous Recombination Deficiency (HRD) Score Predicts Response to Platinum-Containing Neoadjuvant Chemotherapy in Patients with Triple-Negative Breast Cancer. Clin Cancer Res 22(15):3764–3773. https://pubmed.ncbi.nlm.nih.gov/26957554/
2. Abkevich V, Timms KM, Hennessy BT, et al. (2012). Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. Br J Cancer 107(10):1776–1782. https://www.nature.com/articles/bjc2012451
3. Birkbak NJ, Wang ZC, Kim JY, et al. (2012). Telomeric allelic imbalance indicates defective DNA repair and sensitivity to DNA-damaging agents. Cancer Discov 2(4):366–375. https://pmc.ncbi.nlm.nih.gov/articles/PMC3806629/
4. Popova T, Manié E, Rieunier G, et al. (2012). Ploidy and large-scale genomic instability consistently identify basal-like breast carcinomas with BRCA1/2 inactivation. Cancer Res 72(21):5454–5462. https://aacrjournals.org/cancerres/article/72/21/5454/576090/
5. Stewart MD, Merino Vega D, Arend RC, et al. (2022). Homologous Recombination Deficiency: Concepts, Definitions, and Assays. Oncologist 27(3):167–174. https://pmc.ncbi.nlm.nih.gov/articles/PMC8914493/
6. Christinat Y. oncoscanR `score_loh` documentation. https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html
7. Sztupinszki Z, et al. scarHRD reference implementation (R): R/calc.ai_new.R, R/calc.lst.R, R/scar_score.R. https://github.com/sztup/scarHRD (accessed 2026-06-23)
8. UCSC Genome Browser, cytoBand track (acen = centromere) for hg38 and hg19. https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cytoBand and ...genome=hg19... (accessed 2026-06-23)
9. NCBI Genome Reference Consortium — GRCh38 modeled centromeres. https://www.ncbi.nlm.nih.gov/grc/human (accessed 2026-06-23)

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-06-23**: Segment-driven HRD-LOH derivation.
- **2026-06-23**: Segment-driven HRD-TAI and HRD-LST derivation; embedded UCSC cytoBand centromere coordinate table (GRCh38 + GRCh37), cross-verified vs NCBI GRC modeled centromeres; resolved the centromere-table blocker.
