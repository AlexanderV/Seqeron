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
3. **facets-suite denominator (CORRECTED 2026-09-28):** `autosomal_genome = sum(chrom_info$size[chrom_info$chr %in% 1:22])` where `chrom_info` is `sample_chrom_info = get_sample_genome(segs, genome)` (`calculate_fraction_cna`), i.e. per chromosome `size = max(end) − min(start)` of the **sample's own segments** — the interrogated genome. The `genome` build object only supplies `centromere`. The earlier reading ("reference chromosome-size table, NOT the interrogated segments") was wrong. The UCSC autosomal sums above are used only by the explicit reference-assembly option `DetectWholeGenomeDoubling(segments, ReferenceGenome)`; since F32 (2026-10-09) the default `DetectWholeGenomeDoubling(segments)` uses the `get_sample_genome` span. Tables re-verified 2026-09-28 against bedtools2 `genomes/human.hg38.genome` / `human.hg19.genome` (raw.githubusercontent.com; UCSC/Ensembl hosts blocked): all 22 values identical, Σ = 2,875,001,522 / 2,881,033,286.

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

### Dataset: WGD calls against the reference chromosome-size table (explicit `ReferenceGenome` option; facets-suite numerator/threshold, autosomes)

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

**facets-suite call** (default `DetectWholeGenomeDoubling(segments)` since F32, = `DetectWholeGenomeDoublingFromSuppliedLength`; denominator = Σ over autosomes of
max(end) − min(start), `get_sample_genome`). Expected values from a line-by-line Python port of `parse_segs` /
`get_sample_genome` / `is_genome_doubled` (review 2026-09-28):

| Segments (chrom, start, end, A:B) | frac | WGD |
|---|---|---|
| (1, 0, 60M, 2:0), (2, 0, 40M, 1:1) | 0.6 | true |
| (1, 0, 50M, 2:1), (2, 0, 50M, 1:1) | 0.5 | false (strict) |
| (1, 0, 60M, 2:2), (1, 100M, 140M, 1:1) | 0.428571 | false (gap counted in span; old supplied-length rule 0.6 → true) |
| (1, 0, 40M, 2:2), (1, 40M, 100M, 1:1), (X, 0, 100M, 2:2) | 0.4 | false (old rule 0.7 → true) |
| (1, 10M, 70M, 2:2), (1, 70M, 110M, 1:1) | 0.6 | true |
| (1, 0, 60M, 1:2), (2, 0, 40M, 1:1) | 0.6 | true (mcn = tcn − lcn = 3 − 1) |
| (X, 0, 1000, 2:2) only | NaN (R `NA`) | ArgumentException |

### Dataset: default WGD = facets-suite R (F32, review 2026-10-09)

**Source:** mskcc/facets-suite master `R/copy-number-scores.R` (raw.githubusercontent.com, fetched 2026-10-09), sourced
verbatim in R; `calculate_fraction_cna` → `is_genome_doubled(segs, get_sample_genome(segs, genome), treshold = 0.5)`.
`purrr::map_dfr` shimmed as `do.call(rbind, lapply(...))` (purrr unavailable for the installed R); `parse_segs`'
`length = end − start`, `lcn = ifelse(tcn <= 1, 0, lcn)`, `mcn = tcn − lcn` reproduced; chrX given as 23. The
reference-assembly fraction is elevated / 2,875,001,522 (GRCh38). Script: `map_dfr` shim + `source("copy-number-scores.R")`.

| Genome (chrom, start, end, A:B) | elevated | facets span | facets frac (R, %.17g) | facets WGD | GRCh38 frac | GRCh38 WGD |
|---|---|---|---|---|---|---|
| W1 (1, 0, 150M, 2:2), (1, 150M, 248,956,422, 1:1) | 150,000,000 | 248,956,422 | 0.60251508595347658 | TRUE | 0.052173885423077007 | FALSE |
| W2 gapped (1, 0, 80M, 2:2), (1, 120M, 140M, 1:1), (2, 0, 50M, 2:1), (X, 0, 150M, 2:2) | 130,000,000 | 190,000,000 | 0.68421052631578949 | TRUE | 0.04521736736666674 | FALSE |
| W3 every autosome end-to-end, chr1–11 2:2, chr12–22 1:1 | 1,943,767,673 | 2,875,001,522 | 0.67609274573455336 | TRUE | 0.67609274573455336 | TRUE |
| W4 every autosome [10M, size − 10M), chr1–8 2:2, rest 1:1 | 1,376,488,912 | 2,435,001,522 | 0.56529283434263089 | TRUE | 0.47877849853882615 | FALSE |
| W5 (1, 0, 50M, 2:2), (1, 50M, 100M, 1:1) | 50,000,000 | 100,000,000 | 0.5 | FALSE | 0.017391295141025668 | FALSE |
| W6 (1, 0, 60M, 2:2), (1, 100M, 140M, 1:1) | 60,000,000 | 140,000,000 | 0.42857142857142855 | FALSE | 0.020869554169230801 | FALSE |

W1, W2 (gapped) and W4 (telomere-trimmed, whole genome) give different calls under the two denominators; the C#
default and `DetectWholeGenomeDoublingFromSuppliedLength` reproduce the facets column, the `ReferenceGenome.GRCh38`
overload the GRCh38 column (test `DetectWholeGenomeDoubling_Default_MatchesFacetsSuiteR`).

### Dataset: ASCAT `ascat.metrics` WGD / GI / LOH (F30, review 2026-10-09)

Source: ASCAT `R/ascat.metrics.R` (https://raw.githubusercontent.com/VanLoo-lab/ascat/master/ASCAT/R/ascat.metrics.R),
**executed** in R 4.3.3 (file sourced; `ascat.metrics(inputObj, outputObj)` called on minimal ASCAT objects,
`sexchromosomes` = X/Y). Rules read from the code: autosomes only; `modeAllele` weight `(endpos−startpos)/1e6`
(no +1), `round`, cap `y[y>5]=5`, `tapply` groups ascending, stable `order(decreasing=TRUE)` + `which.max` ⇒ ties →
smaller value; WGD mode_majA 0 → NA, 1 → 0, 2 → 1, 3–5 → "1+"; `computeGIscore` size `endpos−startpos+1`,
baseline 1:1 (WGD 0) or 2:2 (WGD 1 and "1+"), `round(…, 4)`; LOH = `round(Σsize[nMinor=0]/Σsize, 4)`.

| Genome (chr, start, end, nMajor:nMinor) | mode_minA | mode_majA | WGD | GI | LOH |
|---|---|---|---|---|---|
| G1 (1,1,50M,1:1) (2,1,40M,1:1) (3,1,30M,2:1) | 1 | 1 | 0 | 0.25 | 0 |
| G2 (1,1,60M,2:2) (2,1,50M,2:0) (3,1,20M,1:1) (4,1,10M,3:2) | 2 | 2 | 1 | 0.5714 | 0.3571 |
| G3 (1,1,70M,3:1) (2,1,30M,2:2) (3,1,20M,2:2) | 1 | 3 | 1+ | 0.5833 | 0 |
| G4 tie (1,0,10M,2:2) (1,10M,20M,2:1) (2,0,20M,1:1) | 1 | 1 | 0 | 0.5 | 0 |
| G5 (1,1,60M,0:0) (2,1,40M,1:1) | 0 | 0 | NA | NA | 0.6 |
| G6 (1,1,30M,1:1) (2,1,20M,2:2) (X,1,150M,2:2) | 1 | 1 | 0 | 0.4 | 0 |
| G7 cap (1,1,8M,6:1) (2,1,8M,7:0) (3,1,10M,2:2) (4,1,5M,1:1) | 1 | 5 | 1+ | 0.6774 | 0.2581 |
| G8 tie 3/4 (1,100,5000100,4:2) (2,100,5000100,3:1) (3,100,3000100,2:2) | 2 | 3 | 1+ | 0.7692 | 0 |
| G9 (chr1,12345,23456789,2:2) (chr2,777,9876543,1:1) (chr5,1,13579246,2:0) (chrY,1,50M,1:0) | 2 | 2 | 1 | 0.5001 | 0.2895 |

Differential check: 3000 random genomes (seed 20261009; 1–12 segments, chr 1–22/X/Y, nMajor 0–8; 108 with an exact
major-mode tie) — `ComputeAscatGenomeMetrics` vs R `ascat.metrics`: **0 mismatches** (all five columns bit-identical).

### Dataset: ASCAT probe-count-weighted ploidy (F31, review 2026-10-09)

Source: ASCAT `R/ascat.runAscat.R` — `ploidy = sum((nA+nB) * s[, "length"]) / sum(s[, "length"])` (runASCAT, l. 283;
`s[, "length"]` = probes per `make_segments` segment) and the reported `ploidy = mean(nA+nB, na.rm=TRUE)` over
probes (ascat.runAscat, l. 98). R 4.3.3:

| Total CN per segment | Probes | R `mean(rep(cn, probes))` | bp-weighted (100/100/50 Mb) |
|---|---|---|---|
| 2, 4, 3 | 1000, 200, 800 | 2.6000000000000001 | 3.0 |
| 2, 4, 3 | 1, 1, 1 | 3 | 3.0 |
| 2, 3, 5 | 1, 1, 1 | 3.3333333333333335 | — |

The same 3000 random genomes (probes 1–500 per segment): `EstimatePloidy(segments, probeCounts)` = R `mean(rep(…))`
= R `sum(cn*len)/sum(len)` bit-for-bit (0 mismatches).

---

## Assumptions

1. **ASSUMPTION: per-segment total copy number is supplied as an allele-specific segment (Major + Minor CN).** The unit reuses the existing `AlleleSpecificSegment` record (ONCO-LOH-001 / ONCO-HRD-001), whose total copy number is `Major + Minor` and whose length is `End − Start`. Patchwork defines ploidy on per-segment **total** copy number; representing total CN as Major+Minor is exactly that total and is also required to evaluate the major-CN≥2 WGD rule. This is an input-shape reuse decision, not an invented numeric constant.
2. **SUPERSEDED (2026-10-09, F32): the default WGD denominator is the facets-suite `get_sample_genome` span; the reference autosomal length below is the explicit `ReferenceGenome` option.** *Historical (2026-06-22):* WGD fraction denominator is now the reference autosomal genome length, not supplied-segment length. facets-suite divides the elevated-major-CN length by `autosomal_genome = sum(chrom_info$size[chr %in% 1:22])` — a reference chromosome-size table. The previous assumption (supplied-segment-length denominator) is replaced by the embedded UCSC `hg38.chrom.sizes` / `hg19.chrom.sizes` tables (cross-verified against Ensembl GRCh38.p14), selected via a `ReferenceGenome { GRCh38, GRCh37 }` parameter (default GRCh38). Only autosomal (chr1–22) segments contribute to the numerator (facets-suite `chrom %in% 1:22`); sex chromosomes and contigs are excluded. The legacy supplied-segment-length behaviour remains available via `DetectWholeGenomeDoublingFromSuppliedLength`. No invented constants: every chromosome length is the published value from the retrieved UCSC table.

---

## Recommendations for Test Coverage

1. **MUST Test:** `EstimatePloidy` on the 3-segment worked example (CN 2/4/3, lengths 100/100/50 Mb) returns exactly 3.0. — Evidence: Patchwork length-weighted mean, ψ = Σ(CN·L)/Σ(L) = 750M/250M.
2. **MUST Test:** `EstimatePloidy` on a pure-diploid genome (all 1:1) returns exactly 2.0. — Evidence: n-scale 2n diploid baseline.
3. **MUST Test:** `EstimatePloidy` is length-weighted, not a plain segment mean — a long CN-2 segment plus a short CN-4 segment must weight toward 2, not toward 3. — Evidence: "weighted by segment length".
4. **MUST Test:** `EstimatePloidy` rejects an empty segment set and any segment with Length ≤ 0 or negative copy number. — Evidence: weighted mean undefined for Σ(L)=0; invalid segment.
5. **MUST Test:** `DetectWholeGenomeDoubling(segments, ReferenceGenome.GRCh38)` flips at the 0.5 boundary computed against the **GRCh38 autosomal genome** (2,875,001,522 bp): true at half+1 bp at major CN ≥ 2, false at exactly half (strict `>`), false at half−1 bp. — Evidence: facets-suite `frac_elevated_mcn > 0.5`, denominator `sum(chrom_info$size[chr %in% 1:22])`.
6. **MUST Test:** the embedded GRCh38 and GRCh37 autosome length tables equal the authoritative UCSC `hg38.chrom.sizes` / `hg19.chrom.sizes` values exactly; the autosomal genome sums equal 2,875,001,522 / 2,881,033,286 bp. — Evidence: UCSC chrom.sizes (Ensembl-cross-verified).
7. **MUST Test:** `DetectWholeGenomeDoubling` uses the **major** allele CN, not total CN: a genome entirely of 1:1 (total 2) segments is NOT doubled. — Evidence: facets-suite `mcn >= 2`.
8. **MUST Test:** a small fully-amplified region (e.g. 100 Mb all major ≥ 2) is NOT WGD under the explicit `ReferenceGenome` option (the default facets call judges it on its own span) (no supplied-segment bias); non-autosomal (chrX/chrY) segments are excluded from the numerator; "chr"-prefixed autosomes are recognised. — Evidence: facets-suite `chrom %in% 1:22`, reference denominator.
9. **SHOULD Test:** `DetectWholeGenomeDoubling` rejects invalid-length / negative-CN / null segments (shared validation); an empty set returns false under the `ReferenceGenome` option (numerator 0 over a fixed denominator) and throws under the default facets call (0/0 → NA, F32). The legacy `DetectWholeGenomeDoublingFromSuppliedLength` retains the supplied-length denominator (60% → true, exactly 50% → false, empty → throws). — Rationale: shared validation; both overloads source-backed.
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
- **2026-10-09**: FIN-B24 F32 — default `DetectWholeGenomeDoubling(segments)` now the facets-suite call (`get_sample_genome` span), R-confirmed on 6 genomes (dataset "default WGD = facets-suite R"); reference-assembly denominator kept as the explicit `DetectWholeGenomeDoubling(segments, ReferenceGenome)` overload.
- **2026-09-28**: Review 2026-09 (B24): facets-suite denominator re-read — `chrom_info = get_sample_genome(segs)` (interrogated span), not a reference table. `DetectWholeGenomeDoublingFromSuppliedLength` made facets-exact (autosome-only numerator and span denominator); mcn = max(Major, Minor); length sums in double. `DetectWholeGenomeDoubling` (reference-table denominator) kept and documented as a deviation.
- **2026-06-22**: Limitation fix (ONCO-PLOIDY-001) — WGD genome fraction now uses a reference chromosome-size table (UCSC hg38/hg19, Ensembl-cross-verified) as the denominator per facets-suite `autosomal_genome`, replacing the supplied-segment-length denominator. Added `ReferenceGenome` selector; resolved Assumption #2; legacy supplied-length behaviour kept as `DetectWholeGenomeDoublingFromSuppliedLength`.
</content>
</invoke>
