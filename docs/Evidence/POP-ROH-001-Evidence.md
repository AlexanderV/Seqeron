# Evidence Artifact: POP-ROH-001

**Test Unit ID:** POP-ROH-001
**Algorithm:** Runs of Homozygosity (ROH) detection and genomic inbreeding coefficient F_ROH
**Date Collected:** 2026-06-13

---

## Online Sources

### McQuillan et al. (2008) — Runs of Homozygosity in European Populations

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2556426/
**Accessed:** 2026-06-13 (fetched full text via WebFetch; PubMed abstract at https://pubmed.ncbi.nlm.nih.gov/18760389/ fetched separately)
**Authority rank:** 1 (peer-reviewed primary paper, Am J Hum Genet)

**Key Extracted Points:**

1. **F_ROH formula (verbatim):** "Froh = ∑Lroh/Lauto", where ∑Lroh is the total length of all an individual's runs of homozygosity above a specified minimum length, and Lauto is the length of the autosomal genome covered by SNPs, excluding centromeres.
2. **L_AUTO numeric value used:** 2,673,768 kb (≈ 2,674 Mb) of SNP-covered autosomal genome.
3. **F_ROH interpretation:** F_ROH is "the proportion of the autosomal genome in runs of homozygosity above a specified length"; it correlated with pedigree-derived inbreeding in Orkney (r = 0.86, from the PubMed abstract fetch).
4. **ROH length thresholds explored:** ≥ 0.5 Mb, ≥ 1.5 Mb, ≥ 5 Mb; the abstract notes ROHs up to ~4 Mb occur in outbred individuals, and 1.5 Mb was particularly informative for distinguishing endogamy levels.

### Chang et al. (2015) — PLINK 1.9 "Runs of homozygosity" documentation

**URL:** https://www.cog-genomics.org/plink/1.9/ibd
**Accessed:** 2026-06-13 (fetched via WebFetch; search query "PLINK --homozyg runs of homozygosity sliding window parameters homozyg-snp homozyg-kb documentation")
**Authority rank:** 3 (reference implementation; PLINK 1.9, Chang CC et al. 2015 GigaScience 4:7)

**Key Extracted Points:**

1. **--homozyg-snp default = 100:** only runs containing at least 100 SNPs are noted.
2. **--homozyg-kb default = 1000:** only runs of total length ≥ 1000 kb (= 1,000,000 bp) are noted.
3. **--homozyg-window-snp default = 50:** the scanning window contains 50 SNPs.
4. **--homozyg-window-het default = 1:** a scanning-window hit can contain at most 1 heterozygous call.
5. **--homozyg-window-missing default = 5:** at most 5 missing calls per window.
6. **--homozyg-window-threshold default = 0.05:** a SNP's scanning-window hit rate must be ≥ 0.05 to be ROH-eligible.
7. **--homozyg-gap default = 1000 kb:** if two consecutive SNPs are > 1000 kb apart they cannot be in the same ROH.
8. **--homozyg-density default = 50 kb/SNP:** a ROH must average at least one SNP per 50 kb.

### Marras et al. (2015) consecutive-runs method — detectRUNS package docs/vignette

**URL:** https://www.rdocumentation.org/packages/detectRUNS/versions/0.9.6/topics/consecutiveRUNS.run and https://cran.r-project.org/web/packages/detectRUNS/vignettes/detectRUNS.vignette.html
**Accessed:** 2026-06-13 (both fetched via WebFetch)
**Authority rank:** 3 (reference implementation of Marras G et al. 2015, Anim Genet 46(2):110-121)

**Key Extracted Points:**

1. **Window-free scan:** the consecutive-runs method "directly scans the genome SNP by SNP"; it examines each SNP to decide whether it belongs to a run (no sliding window).
2. **maxOppRun:** maximum number of opposite (heterozygous) genotypes permitted within a run; exceeding it ends the run.
3. **maxMissRun:** maximum number of missing genotypes allowed; exceeding it ends the run.
4. **minSNP:** minimum number of homozygous SNPs required for a run to be retained (runs below this are discarded).
5. **minLengthBps:** minimum physical length in bp for a valid run.
6. **maxGap:** maximum physical distance between consecutive SNPs; exceeding it breaks the run.
7. **Run termination:** a run continues while SNPs meet the criteria and terminates when any threshold (opposite/missing/gap) is violated; below-minimum stretches are discarded.

---

## Documented Corner Cases and Failure Modes

### From Marras et al. (2015) / detectRUNS

1. **Tolerated opposite genotypes:** a small number (maxOppRun) of heterozygous calls is allowed inside a run to absorb genotyping error; only crossing the threshold breaks the run.
2. **Gap break:** an inter-SNP gap larger than maxGap breaks a run even when all genotypes are homozygous.
3. **Length / count filtering:** stretches not meeting minSNP and minLengthBps are discarded.

### From PLINK 1.9 (Chang et al. 2015)

1. **Two independent thresholds:** both the SNP count (--homozyg-snp) and the physical length (--homozyg-kb) must be satisfied; passing one alone is insufficient.

---

## Test Datasets

### Dataset: McQuillan F_ROH constants

**Source:** McQuillan et al. (2008), PMC2556426

| Parameter | Value |
|-----------|-------|
| L_AUTO (SNP-covered autosome) | 2,673,768 kb |
| F_ROH definition | ΣL_roh / L_auto |
| Example: ΣL_roh = 20 Mb, L_auto = 100 Mb | F_ROH = 0.20 |

### Dataset: PLINK --homozyg default thresholds

**Source:** Chang et al. (2015), PLINK 1.9 documentation

| Parameter | Default |
|-----------|---------|
| min SNPs (--homozyg-snp) | 100 |
| min length (--homozyg-kb) | 1000 kb (1,000,000 bp) |
| max het per window (--homozyg-window-het) | 1 |
| max gap (--homozyg-gap) | 1000 kb (1,000,000 bp) |

---

## Assumptions

1. **ASSUMPTION: Genotype encoding 0/1/2** — The library encodes genotypes as 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate. This matches the additive (allele-dosage) convention used by PLINK `.raw`/`--recodeA` output and elsewhere in `PopulationGeneticsAnalyzer` (e.g., `CalculateMAF`, `CalculateLD`); a "1" is therefore the opposite (heterozygous) genotype that counts against `maxOppRun`. This is an API encoding convention, not a correctness-affecting algorithm parameter — the consecutive-runs rule (Marras 2015) operates on homozygous vs. opposite, independent of how each is encoded.
2. **ASSUMPTION: missing-genotype handling out of scope** — PLINK and detectRUNS also bound missing calls (`maxMissRun` / `--homozyg-window-missing`). The current `FindROH` input is `(Position, Genotype)` with no missing-data sentinel, so missing handling is not modeled; any non-`1` genotype is treated as homozygous. Documented as a limitation, not invented behavior.

---

## Recommendations for Test Coverage

1. **MUST Test:** Single uninterrupted homozygous run reported once with exact Start/End/SnpCount — Evidence: Marras 2015 consecutive scan.
2. **MUST Test:** One tolerated heterozygote (≤ maxOppRun) keeps a single run — Evidence: Marras 2015 maxOppRun.
3. **MUST Test:** A heterozygote beyond tolerance splits into two runs at the right boundaries — Evidence: Marras 2015 run termination.
4. **MUST Test:** Gap > maxGap breaks an all-homozygous run — Evidence: PLINK --homozyg-gap; Marras maxGap.
5. **MUST Test:** Run with < minSnps discarded; run shorter than minLength discarded — Evidence: PLINK --homozyg-snp / --homozyg-kb.
6. **MUST Test:** F_ROH = ΣL_roh / L_auto exact value (0.20 worked example; 1.0 whole-genome) — Evidence: McQuillan 2008.
7. **SHOULD Test:** Unsorted input ordered internally; homozygous-alt (genotype 2) counts as homozygous; leading heterozygotes skipped — Rationale: contract robustness.
8. **COULD Test:** Invalid arguments (null, negative thresholds, non-positive genome length) — Rationale: documented failure modes.

---

## References

1. McQuillan R, Leutenegger A-L, Abdel-Rahman R, et al. (2008). Runs of Homozygosity in European Populations. American Journal of Human Genetics 83(3):359-372. https://pmc.ncbi.nlm.nih.gov/articles/PMC2556426/ (DOI: 10.1016/j.ajhg.2008.08.007)
2. Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience 4:7. PLINK 1.9 "Runs of homozygosity" documentation: https://www.cog-genomics.org/plink/1.9/ibd (DOI: 10.1186/s13742-015-0047-8)
3. Marras G, Gaspa G, Sorbolini S, et al. (2015). Analysis of runs of homozygosity and their relationship with inbreeding in five cattle breeds. Animal Genetics 46(2):110-121. detectRUNS docs: https://cran.r-project.org/web/packages/detectRUNS/vignettes/detectRUNS.vignette.html ; https://www.rdocumentation.org/packages/detectRUNS/versions/0.9.6/topics/consecutiveRUNS.run (DOI: 10.1111/age.12259)

---

## Change History

- **2026-06-13**: Initial documentation.
