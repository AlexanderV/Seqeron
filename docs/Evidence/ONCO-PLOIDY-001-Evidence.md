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

### UCSC Genome Browser — hg38.chrom.sizes / hg19.chrom.sizes (reference chromosome-size tables)

**URL (GRCh38):** https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes
**URL (GRCh37):** https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
**Retrieved by:** WebFetch of the two raw UCSC files above, 2026-06-22, prompting for the exact integer base-pair size of chr1–chr22, chrX, chrY, chrM.
**Cross-verification (GRCh38):** Ensembl REST `https://rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json` (WebFetch, 2026-06-22) — assembly "GRCh38.p14"; chr1 = 248,956,422; chr21 = 46,709,983; chr22 = 50,818,468; chrX = 156,040,895 — identical to the UCSC values.
**Accessed:** 2026-06-22
**Authority rank:** 5 (well-maintained genome database, UCSC/Ensembl assembly metadata) — the canonical published chromosome lengths of the named human reference assemblies.

**Key Extracted Points:**

1. **GRCh38 / hg38 autosome lengths (bp), chr1…chr22:** 248,956,422; 242,193,529; 198,295,559; 190,214,555; 181,538,259; 170,805,979; 159,345,973; 145,138,636; 138,394,717; 133,797,422; 135,086,622; 133,275,309; 114,364,328; 107,043,718; 101,991,189; 90,338,345; 83,257,441; 80,373,285; 58,617,616; 64,444,167; 46,709,983; 50,818,468. **Σ(chr1–22) = 2,875,001,522 bp** (the WGD denominator).
2. **GRCh37 / hg19 autosome lengths (bp), chr1…chr22:** 249,250,621; 243,199,373; 198,022,430; 191,154,276; 180,915,260; 171,115,067; 159,138,663; 146,364,022; 141,213,431; 135,534,747; 135,006,516; 133,851,895; 115,169,878; 107,349,540; 102,531,392; 90,354,753; 81,195,210; 78,077,248; 59,128,983; 63,025,520; 48,129,895; 51,304,566. **Σ(chr1–22) = 2,881,033,286 bp**.
3. **facets-suite denominator:** `autosomal_genome = sum(chrom_info$size[chrom_info$chr %in% 1:22])` — exactly the autosomal sum above, taken from a reference chromosome-size table (parameterised by genome build hg19/hg18/hg38), NOT from the interrogated segments.

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

### Dataset: WGD calls against the reference chromosome-size table (facets-suite is_genome_doubled, autosomes)

**Source:** facets-suite `copy-number-scores.R` (`autosomal_genome = sum(chrom_info$size[chr %in% 1:22])`,
`frac_elevated_mcn = sum(length where mcn≥2 & chrom %in% 1:22) / autosomal_genome`, `wgd = frac_elevated_mcn > 0.5`);
PMID 30013179; UCSC hg38.chrom.sizes. **GRCh38 autosomal genome = 2,875,001,522 bp; half = 1,437,500,761 bp.**

| Case (denominator = GRCh38 Σchr1–22) | Major-CN≥2 autosomal length / 2,875,001,522 | WGD? |
|--------------------------------------|---------------------------------------------|------|
| 1,437,500,762 bp at major CN ≥ 2 (half + 1) | > 0.5 (strict) | **true** (doubled) |
| 1,437,500,761 bp at major CN ≥ 2 (exactly half) | = 0.5, not > 0.5 | **false** (boundary, strict) |
| 1,437,500,760 bp at major CN ≥ 2 (half − 1) | < 0.5 | **false** |
| 100 Mb fully amplified (major ≥ 2), genome not tiled | 100M / 2.875G ≈ 0.035 | **false** (no supplied-segment bias) |
| all 1:1 autosomal segments (major = 1) | 0.0 | **false** |
| chrX/chrY amplified, no autosomal elevation | 0.0 (sex chromosomes excluded) | **false** |

**Legacy supplied-segment-length variant** (`DetectWholeGenomeDoublingFromSuppliedLength`, denominator = Σ supplied
segment length): retains the pre-fix semantics — 60% of supplied length at major CN ≥ 2 → true; exactly 50% → false.

---

## Assumptions

1. **ASSUMPTION: per-segment total copy number is supplied as an allele-specific segment (Major + Minor CN).** The unit reuses the existing `AlleleSpecificSegment` record (ONCO-LOH-001 / ONCO-HRD-001), whose total copy number is `Major + Minor` and whose length is `End − Start`. Patchwork defines ploidy on per-segment **total** copy number; representing total CN as Major+Minor is exactly that total and is also required to evaluate the major-CN≥2 WGD rule. This is an input-shape reuse decision, not an invented numeric constant.
2. **RESOLVED (2026-06-22): WGD fraction denominator is now the reference autosomal genome length, not supplied-segment length.** facets-suite divides the elevated-major-CN length by `autosomal_genome = sum(chrom_info$size[chr %in% 1:22])` — a reference chromosome-size table. The previous assumption (supplied-segment-length denominator) is replaced by the embedded UCSC `hg38.chrom.sizes` / `hg19.chrom.sizes` tables (cross-verified against Ensembl GRCh38.p14), selected via a `ReferenceGenome { GRCh38, GRCh37 }` parameter (default GRCh38). Only autosomal (chr1–22) segments contribute to the numerator (facets-suite `chrom %in% 1:22`); sex chromosomes and contigs are excluded. The legacy supplied-segment-length behaviour remains available via `DetectWholeGenomeDoublingFromSuppliedLength`. No invented constants: every chromosome length is the published value from the retrieved UCSC table.

---

## Recommendations for Test Coverage

1. **MUST Test:** `EstimatePloidy` on the 3-segment worked example (CN 2/4/3, lengths 100/100/50 Mb) returns exactly 3.0. — Evidence: Patchwork length-weighted mean, ψ = Σ(CN·L)/Σ(L) = 750M/250M.
2. **MUST Test:** `EstimatePloidy` on a pure-diploid genome (all 1:1) returns exactly 2.0. — Evidence: n-scale 2n diploid baseline.
3. **MUST Test:** `EstimatePloidy` is length-weighted, not a plain segment mean — a long CN-2 segment plus a short CN-4 segment must weight toward 2, not toward 3. — Evidence: "weighted by segment length".
4. **MUST Test:** `EstimatePloidy` rejects an empty segment set and any segment with Length ≤ 0 or negative copy number. — Evidence: weighted mean undefined for Σ(L)=0; invalid segment.
5. **MUST Test:** `DetectWholeGenomeDoubling` flips at the 0.5 boundary computed against the **GRCh38 autosomal genome** (2,875,001,522 bp): true at half+1 bp at major CN ≥ 2, false at exactly half (strict `>`), false at half−1 bp. — Evidence: facets-suite `frac_elevated_mcn > 0.5`, denominator `sum(chrom_info$size[chr %in% 1:22])`.
6. **MUST Test:** the embedded GRCh38 and GRCh37 autosome length tables equal the authoritative UCSC `hg38.chrom.sizes` / `hg19.chrom.sizes` values exactly; the autosomal genome sums equal 2,875,001,522 / 2,881,033,286 bp. — Evidence: UCSC chrom.sizes (Ensembl-cross-verified).
7. **MUST Test:** `DetectWholeGenomeDoubling` uses the **major** allele CN, not total CN: a genome entirely of 1:1 (total 2) segments is NOT doubled. — Evidence: facets-suite `mcn >= 2`.
8. **MUST Test:** a small fully-amplified region (e.g. 100 Mb all major ≥ 2) is NOT WGD against the reference genome (no supplied-segment bias); non-autosomal (chrX/chrY) segments are excluded from the numerator; "chr"-prefixed autosomes are recognised. — Evidence: facets-suite `chrom %in% 1:22`, reference denominator.
9. **SHOULD Test:** `DetectWholeGenomeDoubling` rejects invalid-length / negative-CN / null segments (shared validation); an empty set returns false (numerator 0 over a fixed denominator). The legacy `DetectWholeGenomeDoublingFromSuppliedLength` retains the supplied-length denominator (60% → true, exactly 50% → false, empty → throws). — Rationale: shared validation; both overloads source-backed.
10. **COULD Test:** the GRCh37 selector uses the hg19 denominator and can disagree with GRCh38 near the boundary. — Rationale: build-dependent denominator (facets-suite `genome` parameter).

---

## References

1. Mayrhofer M, Viklund B, Isaksson A (2016 reprint of 2014 method). [Patchwork as PMC4053982] Allele-specific copy number analysis of whole-genome sequenced tumor tissue. *Genome Biology*. https://pmc.ncbi.nlm.nih.gov/articles/PMC4053982/ (accessed 2026-06-14) — verbatim ploidy definition.
2. Van Loo P, Nordgard SH, Lingjærde OC, et al. (2010). Allele-specific copy number analysis of tumors. *PNAS* 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
3. Bielski CM, Zehir A, Penson AV, et al. (2018). Genome doubling shapes the evolution and prognosis of advanced cancers. *Nature Genetics* 50(8):1189–1195. https://doi.org/10.1038/s41588-018-0165-1
4. facets-suite (MSKCC), `R/copy-number-scores.R`, `is_genome_doubled` (treshold = 0.5, mcn = tcn − lcn, `autosomal_genome = sum(chrom_info$size[chr %in% 1:22])`, PMID 30013179). https://github.com/mskcc/facets-suite/blob/master/R/copy-number-scores.R (accessed 2026-06-14, re-fetched 2026-06-22 for the denominator)
5. UCSC Genome Browser, `hg38.chrom.sizes` (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes) and `hg19.chrom.sizes` (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes) — reference chromosome lengths (accessed 2026-06-22).
6. Ensembl REST, assembly metadata for *Homo sapiens* GRCh38.p14 (https://rest.ensembl.org/info/assembly/homo_sapiens) — cross-verification of GRCh38 chromosome lengths (accessed 2026-06-22).

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-06-22**: Limitation fix (ONCO-PLOIDY-001) — WGD genome fraction now uses a reference chromosome-size table (UCSC hg38/hg19, Ensembl-cross-verified) as the denominator per facets-suite `autosomal_genome`, replacing the supplied-segment-length denominator. Added `ReferenceGenome` selector; resolved Assumption #2; legacy supplied-length behaviour kept as `DetectWholeGenomeDoublingFromSuppliedLength`.
</content>
</invoke>
