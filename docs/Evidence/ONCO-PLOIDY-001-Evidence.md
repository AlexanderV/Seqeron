# Evidence Artifact: ONCO-PLOIDY-001

**Test Unit ID:** ONCO-PLOIDY-001
**Algorithm:** Tumor Ploidy Estimation (length-weighted mean segment copy number) and Whole-Genome-Doubling detection
**Date Collected:** 2026-06-14

---

## Online Sources

### Patchwork — allele-specific copy number analysis of whole-genome sequenced tumor tissue (Genome Biology 2010; PMC4053982)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4053982/
**Retrieved by:** WebFetch of the URL above, 2026-06-14, prompting for the definition of average tumour ploidy as a segment-length-weighted mean of total copy number.
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Genome Biology).

**Key Extracted Points:**

1. **Average-ploidy definition (verbatim):** "The average ploidy, PloidyTum, is the average total copy number of all genomic segments weighted by segment length." This is the length-weighted mean ψ = Σ(CN_i · L_i) / Σ(L_i) over segments.
2. **Total copy number per segment:** the per-segment quantity averaged is the segment **total** copy number (sum of allele copy numbers), not an allele-specific value.

### ASCAT — Allele-specific copy number analysis of tumors (Van Loo et al., PNAS 2010; PMID 20837533)

**URL / retrieval:** Citation verified via Europe PMC REST core record
`https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:20837533&format=json&resultType=core` (WebFetch), 2026-06-14. (PNAS HTML returned HTTP 403; the PDF was non-extractable binary, so the verified bibliographic metadata and abstract were taken from the Europe PMC core record.)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, PNAS) — the originating method that infers tumour purity and **ploidy** with allele-specific segment copy numbers.

**Key Extracted Points:**

1. **Verified citation:** Van Loo P, Nordgard SH, Lingjærde OC, Russnes HG, Rye IH, Sun W, Weigman VJ, Marynen P, Zetterberg A, Naume B, Perou CM, Børresen-Dale AL, Kristensen VN. "Allele-specific copy number analysis of tumors." *PNAS* 107(39):16910–16915, 2010. DOI: 10.1073/pnas.1009843107.
2. **Aneuploidy threshold (verbatim abstract):** "We observe aneuploidy (>2.7n) in 45% of the cases" — establishes that average tumour ploidy is reported on the n-scale (2n = diploid) and that elevated whole-genome ploidy (here >2.7) marks aneuploidy/near-triploid genomes (basal-like breast carcinomas yield "near-triploid genomes"). ASCAT outputs a final tumour `ploidy` field (`ascat.output$ploidy`).

### facets-suite — copy-number-scores.R (MSKCC reference implementation of the WGD rule)

**URL:** https://raw.githubusercontent.com/mskcc/facets-suite/master/R/copy-number-scores.R
**Retrieved by:** WebFetch of the raw GitHub source above, 2026-06-14, prompting for the verbatim `is_genome_doubled` function, the elevated-major-CN fraction, the 0.5 threshold, and the Bielski PMID.
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation in an established bioinformatics library, encoding the rank-1 Bielski 2018 definition).

**Key Extracted Points:**

1. **WGD rule (verbatim code):**
   ```r
   is_genome_doubled = function(segs, chrom_info, treshold = 0.5) {
       autosomal_genome = sum(as.numeric(chrom_info$size[chrom_info$chr %in% 1:22]))
       # Check for whole-genome duplication // PMID 30013179
       frac_elevated_mcn = sum(as.numeric(segs$length[which(segs$mcn >= 2 & segs$chrom %in% 1:22)])) / autosomal_genome
       wgd = frac_elevated_mcn > treshold
       wgd
   }
   ```
   i.e. WGD is called when the autosome-restricted fraction of the genome with **major copy number ≥ 2** is **strictly greater than 0.5**.
2. **Major copy number (verbatim):** `mcn = tcn - lcn` (major CN = total CN − minor/lesser CN). In the allele-specific segment record used here this is the larger allele copy number directly.
3. **Bielski attribution (verbatim comment):** `# Check for whole-genome duplication // PMID 30013179` → Bielski et al. 2018.

### Bielski et al. — Genome doubling shapes the evolution and prognosis of advanced cancers (Nature Genetics 2018; PMID 30013179)

**URL / retrieval:** Citation verified via Europe PMC REST core record
`https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:30013179&format=json&resultType=core` (WebFetch), 2026-06-14. (PubMed HTML returned a reCAPTCHA screen; the verified bibliographic metadata and abstract were taken from the Europe PMC core record. The operational ≥50%/MCN≥2 threshold is taken from the facets-suite reference implementation above, which cites this PMID.)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Nature Genetics).

**Key Extracted Points:**

1. **Verified citation:** Bielski CM, Zehir A, Penson AV, Donoghue MTA, Chatila W, Armenia J, Chang MT, Schram AM, Jonsson P, Bandlamudi C, Razavi P, Iyer G, Robson ME, Stadler ZK, Schultz N, Baselga J, Solit DB, Hyman DM, Berger MF, Taylor BS. "Genome doubling shapes the evolution and prognosis of advanced cancers." *Nature Genetics* 50(8):1189–1195, 2018. DOI: 10.1038/s41588-018-0165-1.
2. **Context (verbatim abstract):** "we identified whole-genome doubling (WGD) in the tumors of nearly 30% of 9,692 prospectively sequenced advanced cancer patients" — WGD is a per-sample binary classification of genome state, operationalised by the facets-suite ≥50% major-CN≥2 rule.

---

## Documented Corner Cases and Failure Modes

### From Patchwork / ASCAT

1. **Empty segment set:** ploidy is the weighted mean over segments; with no segments Σ(L_i) = 0 and ploidy is undefined (division by zero) → must be rejected.
2. **Zero / negative segment length:** a segment with Length ≤ 0 carries no genomic weight and corrupts the weighted mean; such input is invalid.
3. **Diploid genome baseline:** a genome of pure copy-number-2 segments has ploidy exactly 2.0 (n-scale); elevated ploidy (>2.7n per Van Loo) indicates aneuploidy.

### From facets-suite

1. **Autosome restriction:** the WGD fraction is computed over autosomes (chromosomes 1–22) only; sex chromosomes are excluded so a single-X male genome does not bias the doubling call.
2. **Strict threshold (`> 0.5`, not `≥`):** exactly half of the autosomal genome at major CN ≥ 2 is NOT doubled; it must be strictly more than half.
3. **Major CN vs total CN:** doubling uses the **major** allele copy number (mcn = tcn − lcn) ≥ 2, not total CN ≥ 2 — a 1:1 → tcn 2 segment is NOT elevated (major = 1), whereas a 2:0 (LOH) or 2:1 segment IS (major = 2).

---

## Test Datasets

### Dataset: Length-weighted ploidy worked example (derived from Patchwork definition)

**Source:** Patchwork (Genome Biology 2010), PMC4053982 — "average total copy number of all genomic segments weighted by segment length."

| Segment | Total CN (Major+Minor) | Length (bp) | CN × Length |
|---------|------------------------|-------------|-------------|
| A | 2 (1:1) | 100,000,000 | 200,000,000 |
| B | 4 (2:2) | 100,000,000 | 400,000,000 |
| C | 3 (2:1) | 50,000,000  | 150,000,000 |

Σ(CN·L) = 750,000,000; Σ(L) = 250,000,000 → ψ = 750,000,000 / 250,000,000 = **3.0**.

### Dataset: Pure-diploid genome (identity case)

**Source:** ASCAT/Patchwork n-scale (2n = diploid).

| Segment | Total CN | Length (bp) |
|---------|----------|-------------|
| all | 2 (1:1) | any positive | → ψ = **2.0** exactly. |

### Dataset: WGD calls (facets-suite is_genome_doubled rule, threshold 0.5, autosomes)

**Source:** facets-suite `copy-number-scores.R`; PMID 30013179.

| Case | Major-CN≥2 length / autosomal length | WGD? |
|------|--------------------------------------|------|
| 60% of autosomal genome at major CN ≥ 2 | 0.60 > 0.50 | **true** (doubled) |
| exactly 50% at major CN ≥ 2 | 0.50, not > 0.50 | **false** (boundary, strict) |
| 40% at major CN ≥ 2 | 0.40 ≤ 0.50 | **false** |
| all segments 1:1 (major = 1) | 0.0 | **false** |
| all segments 2:0/2:2 (major ≥ 2) | 1.0 > 0.50 | **true** |

---

## Assumptions

1. **ASSUMPTION: per-segment total copy number is supplied as an allele-specific segment (Major + Minor CN).** The unit reuses the existing `AlleleSpecificSegment` record (ONCO-LOH-001 / ONCO-HRD-001), whose total copy number is `Major + Minor` and whose length is `End − Start`. Patchwork defines ploidy on per-segment **total** copy number; representing total CN as Major+Minor is exactly that total and is also required to evaluate the major-CN≥2 WGD rule. This is an input-shape reuse decision, not an invented numeric constant.
2. **ASSUMPTION: WGD fraction is computed over the supplied segments' total length.** facets-suite divides the elevated-major-CN length by the **autosomal genome** length from a chromosome-size table. With no external chromosome-size table in scope, the denominator is the total length of the supplied segments (the interrogated genome). This preserves the fraction semantics exactly when the caller supplies autosomal segments covering the genome, and the ≥2/>0.5 rule is unchanged. Recorded for transparency; it does not alter the threshold or comparison operator.

---

## Recommendations for Test Coverage

1. **MUST Test:** `EstimatePloidy` on the 3-segment worked example (CN 2/4/3, lengths 100/100/50 Mb) returns exactly 3.0. — Evidence: Patchwork length-weighted mean, ψ = Σ(CN·L)/Σ(L) = 750M/250M.
2. **MUST Test:** `EstimatePloidy` on a pure-diploid genome (all 1:1) returns exactly 2.0. — Evidence: n-scale 2n diploid baseline.
3. **MUST Test:** `EstimatePloidy` is length-weighted, not a plain segment mean — a long CN-2 segment plus a short CN-4 segment must weight toward 2, not toward 3. — Evidence: "weighted by segment length".
4. **MUST Test:** `EstimatePloidy` rejects an empty segment set and any segment with Length ≤ 0 or negative copy number. — Evidence: weighted mean undefined for Σ(L)=0; invalid segment.
5. **MUST Test:** `DetectWholeGenomeDoubling` returns true when 60% of length has major CN ≥ 2, false at exactly 50% (strict `>`), false at 40%, false for all-1:1, true for all major-CN≥2. — Evidence: facets-suite `frac_elevated_mcn > 0.5`, mcn = tcn − lcn.
6. **MUST Test:** `DetectWholeGenomeDoubling` uses the **major** allele CN, not total CN: a genome entirely of 1:1 (total 2) segments is NOT doubled. — Evidence: facets-suite `mcn >= 2`.
7. **SHOULD Test:** `DetectWholeGenomeDoubling` rejects empty / invalid-length segments consistently with `EstimatePloidy`. — Rationale: shared validation, fraction denominator undefined.
8. **COULD Test:** a ploidy-based overload of the WGD flag (ploidy above a stated cutoff) agrees with the n-scale aneuploidy direction (>2.7n elevated). — Rationale: Van Loo aneuploidy direction (kept as a derived convenience only if it is itself source-backed).

---

## References

1. Mayrhofer M, Viklund B, Isaksson A (2016 reprint of 2014 method). [Patchwork as PMC4053982] Allele-specific copy number analysis of whole-genome sequenced tumor tissue. *Genome Biology*. https://pmc.ncbi.nlm.nih.gov/articles/PMC4053982/ (accessed 2026-06-14) — verbatim ploidy definition.
2. Van Loo P, Nordgard SH, Lingjærde OC, et al. (2010). Allele-specific copy number analysis of tumors. *PNAS* 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
3. Bielski CM, Zehir A, Penson AV, et al. (2018). Genome doubling shapes the evolution and prognosis of advanced cancers. *Nature Genetics* 50(8):1189–1195. https://doi.org/10.1038/s41588-018-0165-1
4. facets-suite (MSKCC), `R/copy-number-scores.R`, `is_genome_doubled` (treshold = 0.5, mcn = tcn − lcn, PMID 30013179). https://github.com/mskcc/facets-suite/blob/master/R/copy-number-scores.R (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
</content>
</invoke>
