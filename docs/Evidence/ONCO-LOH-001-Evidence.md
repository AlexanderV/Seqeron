# Evidence Artifact: ONCO-LOH-001

**Test Unit ID:** ONCO-LOH-001
**Algorithm:** Loss of Heterozygosity (LOH) detection and HRD-LOH genomic-scar score
**Date Collected:** 2026-06-14

---

## Online Sources

### Abkevich et al. (2012), Br J Cancer — HRD-LOH primary definition

**URL:** https://pubmed.ncbi.nlm.nih.gov/23047548/
**Accessed:** 2026-06-14 (retrieved via WebFetch on the PubMed abstract page)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **HRD score definition (verbatim):** "The homologous recombination deficiency (HRD) score was defined as the number of these regions observed in a tumour sample." — where "these regions" are the intermediate-size LOH regions established earlier in the abstract.
2. **Intermediate-size LOH correlates with HR deficiency (verbatim):** "Loss of heterozygosity regions of intermediate size were observed more frequently in tumours with defective BRCA1 or BRCA2" (P=10⁻¹¹). → only intermediate-size LOH (not whole-chromosome LOH) is counted.
3. **Full citation:** Abkevich V, Timms KM, Hennessy BT, et al. Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. Br J Cancer. 2012;107(10):1776–1782. doi:10.1038/bjc.2012.451. PMID:23047548.

### scarHRD (sztup) — reference implementation: HRD-LOH definition and size cutoff

**URL:** https://github.com/sztup/scarHRD/blob/master/scarHRD.md (raw: https://raw.githubusercontent.com/sztup/scarHRD/master/scarHRD.md)
**Accessed:** 2026-06-14 (retrieved via WebFetch)
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **HRD-LOH definition (verbatim):** "The definition of a sample's HRD-LOH score is the number of 15 Mb exceeding LOH regions which do not cover the whole chromosome." → count LOH regions with length > 15 Mb that are not whole-chromosome.
2. **15 Mb cutoff rationale (verbatim):** "the cut-off of 15 Mb approximately in the middle of the interval was arbitrarily selected for further analysis."
3. **Abkevich citation reproduced:** "Abkevich, V., K. M. Timms, B. T. Hennessy, et al. 2012. ... Br. J. Cancer 107 (10): 1776–82."

### scarHRD (sztup) — reference implementation: `calc.hrd` source (exact counting algorithm)

**URL:** https://github.com/sztup/scarHRD/blob/master/R/calc.hrd.R (raw: https://raw.githubusercontent.com/sztup/scarHRD/master/R/calc.hrd.R)
**Accessed:** 2026-06-14 (retrieved via WebFetch)
**Authority rank:** 3 (reference implementation)

**Key Extracted Points (verbatim R code):**

1. **Whole-chromosome exclusion:** for each chromosome `j`, `if(all(segSamp[segSamp[,2] == j,nB] == 0)) chrDel <- c(chrDel, j)` — a chromosome where ALL segments have minor copy number `nB == 0` is a "global LOH" chromosome and is excluded.
2. **LOH segment condition (verbatim):** `segLOH <- segSamp[segSamp[,nB] == 0 & segSamp[,nA] != 0,,drop=F]` — a segment is LOH when the **minor** allele copy number (`nB`) is 0 AND the **major** allele copy number (`nA`) is non-zero (so it is allelic loss, not a homozygous deletion).
3. **Size filter (verbatim):** `segLOH <- segLOH[segLOH[,4]-segLOH[,3] > sizelimit1,,drop=F]` — keep only segments whose length (end − start, columns 4 − 3) is **strictly greater than** `sizelimit1`.
4. **Apply chromosome exclusion (verbatim):** `segLOH <- segLOH[!segLOH[,2] %in% chrDel,,drop=F]` then `output[i] <- nrow(segLOH)` — drop LOH segments on global-LOH chromosomes; the score is the count of remaining LOH segments.
5. **Major-allele cap before merge (verbatim):** `segSamp[segSamp[,nA] > 1,nA] <- 1` then `segSamp <- shrink.seg.ai.wrapper(segSamp)` — major copy number > 1 is capped to 1 before adjacent equal-state segments are merged (so the LOH/non-LOH state, not absolute copy number, drives merging).

### scarHRD (sztup) — reference implementation: input column order and `sizelimitLOH`

**URL:** https://github.com/sztup/scarHRD/blob/master/R/scar_score.R (raw: https://raw.githubusercontent.com/sztup/scarHRD/master/R/scar_score.R)
**Accessed:** 2026-06-14 (retrieved via WebFetch)
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Input columns (verbatim doc):** "1st column: sample name, 2nd column: chromosome, 3rd column: segmentation start, 4th column: segmentation end, 5th column: total copynumber, 6th column: copy number of A allele, 7th column: copy number of B allele."
2. **Column shift before scoring (verbatim):** `seg[,9]<-seg[,8]; seg[,8]<-seg[,7]; seg[,7]<-seg[,6]` — so within `calc.hrd` column 7 (`nA`) = A allele (major) and column 8 (`nB`) = B allele (minor) copy number.
3. **Size limit default:** `sizelimitLOH = 15e6` (15,000,000 bp) passed as `calc.hrd(seg, sizelimit1 = sizelimitLOH)`.

### oncoscanR (Christinat) `score_loh` — independent corroboration of the counting rule

**URL:** https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html
**Accessed:** 2026-06-14 (retrieved via WebFetch)
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Counting rule (verbatim):** "All LOH segments larger than 15Mb but excluding chromosome with a global LOH alteration."
2. **Merging (verbatim):** "merges overlapping or adjacent LOH segments (separated by 1bp)."
3. **Scientific basis (verbatim):** "Procedure based on the paper from Abkevich et al., Br J Cancer 2012 (PMID: 23047548)."

---

## Documented Corner Cases and Failure Modes

### From scarHRD `calc.hrd`

1. **Homozygous deletion is NOT LOH:** a segment with minor == 0 AND major == 0 fails the `nA != 0` condition; it is not counted (both alleles lost, not loss of *heterozygosity*).
2. **Whole-chromosome LOH excluded:** if every segment on a chromosome has minor == 0, the whole chromosome is in `chrDel` and none of its LOH segments are counted (Abkevich: whole-chromosome LOH does not correlate with HRD).
3. **Strict size boundary:** the comparison is `> sizelimit1`, so a segment of length exactly 15,000,000 bp is NOT counted; only length strictly greater than 15 Mb is counted.
4. **Length = end − start:** length is computed as end − start (columns 4 − 3), not end − start + 1.

### From oncoscanR `score_loh`

1. **Adjacent same-state merge:** overlapping or adjacent (≤ 1 bp gap) LOH segments are merged before the size filter, so two adjacent < 15 Mb LOH pieces forming a > 15 Mb LOH region count as one region.

---

## Test Datasets

### Dataset: scarHRD `calc.hrd` rule applied to synthetic allele-specific segments

**Source:** scarHRD `calc.hrd` (sztup), R/calc.hrd.R — the LOH condition `minor==0 & major!=0`, size filter `end-start > 15e6`, whole-chromosome exclusion. Lengths chosen to straddle the 15 Mb boundary.

| Chr | Start | End | Major (A) | Minor (B) | Length (bp) | LOH? | Counted? | Reason |
|-----|-------|-----|-----------|-----------|-------------|------|----------|--------|
| 1 | 0 | 20,000,000 | 1 | 0 | 20,000,000 | yes | **yes** | minor=0, major≠0, >15 Mb, chr1 has a non-LOH segment |
| 1 | 20,000,000 | 60,000,000 | 1 | 1 | 40,000,000 | no | no | minor≠0 (het retained) |
| 2 | 0 | 10,000,000 | 2 | 0 | 10,000,000 | yes | no | LOH but length ≤ 15 Mb |
| 3 | 0 | 16,000,000 | 1 | 0 | 16,000,000 | yes | no | whole chr3 (all segments minor=0) → excluded |
| 4 | 0 | 30,000,000 | 0 | 0 | 30,000,000 | no | no | homozygous deletion (major=0) |
| 5 | 0 | 15,000,000 | 1 | 0 | 15,000,000 | yes | no | length exactly 15 Mb → not > 15 Mb |
| 5 | 15,000,000 | 50,000,000 | 1 | 1 | 35,000,000 | no | no | minor≠0 |

**Expected HRD-LOH score for this dataset = 1** (only the chr1 20 Mb LOH segment qualifies).

### Dataset: per-chromosome LOH fraction (invariant 0 ≤ fraction ≤ 1)

**Source:** Derived from the LOH-segment definition (length-weighted fraction of a chromosome under LOH). Bounds follow because LOH-segment lengths are a subset of total covered length.

| Chr | Segments (start–end, minor) | LOH bp | Total bp | Fraction |
|-----|------------------------------|--------|----------|----------|
| 1 | (0–20M, m=0), (20M–60M, m=1) | 20M | 60M | 0.333… |
| 2 | (0–50M, m=1) | 0 | 50M | 0.0 |
| 3 | (0–40M, m=0) | 40M | 40M | 1.0 |

---

## Assumptions

1. **ASSUMPTION: LOH-fraction definition** — Abkevich 2012 / scarHRD define only the *count* of qualifying LOH regions (the HRD-LOH score), not a per-chromosome "LOH fraction". The Registry lists `CalculateLOHFraction(chromosome)` with invariant `0 ≤ LOH_fraction ≤ 1`. We define LOH fraction as the length-weighted fraction of a chromosome's covered length that lies under LOH segments (minor==0 & major!=0), which is the natural quantity satisfying the stated invariant and is consistent with the segment-based LOH definition. The *correctness-affecting* LOH-segment criterion itself is fully source-backed (minor==0 & major!=0); only the aggregation into a 0–1 fraction is a definitional/API choice and is not contradicted by any source.
2. **ASSUMPTION: input shape** — `DetectLOH` takes allele-specific copy-number **segments** (chromosome, start, end, major CN, minor CN), i.e. the scarHRD `seg`-table shape, rather than raw `tumorVcf, normalVcf` text. The retrievable algorithm (Abkevich/scarHRD/oncoscanR) operates on segmented allele-specific copy number; raw-VCF parsing and the upstream segmentation/BAF model are out of scope for this unit (handled upstream, cf. ONCO-CNA-001). This is an API-shape decision; the counting logic is unchanged.

---

## Recommendations for Test Coverage

1. **MUST Test:** HRD-LOH count on the synthetic dataset = 1 — Evidence: scarHRD `calc.hrd` (LOH condition + 15 Mb filter + whole-chromosome exclusion).
2. **MUST Test:** segment of length exactly 15 Mb is NOT counted (strict `>` boundary) — Evidence: `segLOH[,4]-segLOH[,3] > sizelimit1`.
3. **MUST Test:** homozygous deletion (minor=0, major=0) is NOT counted — Evidence: `nA != 0` clause.
4. **MUST Test:** whole-chromosome LOH (all segments minor=0) excluded — Evidence: `chrDel` logic / Abkevich "< whole chromosome".
5. **MUST Test:** a heterozygous-retained segment (minor≠0) is NOT LOH — Evidence: `nB == 0` clause.
6. **MUST Test:** `CalculateLOHFraction` returns length-weighted LOH fraction in [0,1] for representative chromosomes — Evidence: invariant + LOH definition.
7. **SHOULD Test:** adjacent LOH segments forming > 15 Mb after merge count as one region — Evidence: oncoscanR "merges adjacent LOH segments".
8. **SHOULD Test:** null / empty segment input handled (empty → score 0, fraction 0).
9. **COULD Test:** ordering invariance — shuffling input segment order does not change the count.

---

## References

1. Abkevich V, Timms KM, Hennessy BT, et al. (2012). Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. Br J Cancer 107(10):1776–1782. https://doi.org/10.1038/bjc.2012.451 (PMID:23047548, https://pubmed.ncbi.nlm.nih.gov/23047548/)
2. Sztupinszki Z, et al. scarHRD R package — `calc.hrd` and `scar_score`. https://github.com/sztup/scarHRD/blob/master/R/calc.hrd.R , https://github.com/sztup/scarHRD/blob/master/scarHRD.md
3. Christinat Y. oncoscanR `score_loh` documentation. https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html

---

## Change History

- **2026-06-14**: Initial documentation.
