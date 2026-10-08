# Evidence Artifact: ONCO-SOMATIC-001

**Test Unit ID:** ONCO-SOMATIC-001
**Algorithm:** Somatic Mutation Calling (tumor vs matched normal classification)
**Date Collected:** 2026-06-14

---

## Online Sources

### Strelka (Saunders et al. 2012) — original somatic small-variant caller

**URL:** https://academic.oup.com/bioinformatics/article/28/14/1811/218573 (citation confirmed via https://pubmed.ncbi.nlm.nih.gov/22581179/)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Bioinformatics)

**Retrieved how:** WebSearch "Strelka2 somatic variant calling tumor normal allele frequency definition paper" → WebFetch of the Oxford Academic article page; full citation cross-checked by WebFetch of the PubMed record (PMID 22581179).

**Key Extracted Points:**

1. **Formal somatic state:** The paper defines the somatic state as **S = {(f_t, f_n): f_t ≠ f_n}** — the tumor allele frequency f_t and normal allele frequency f_n must differ. (Quoted from fetched article text.)
2. **Normal genotype restriction:** Calls are restricted to a **homozygous-reference normal genotype**: quality scores reflect the joint **P(S, G_n = 'ref/ref' | D)**. The raw somatic probability alone is not used because it detects many variants in LOH / copy-number regions.
3. **Generative model:** The normal sample is modeled as a mixture of diploid germline variation with noise; the tumor sample is modeled as a mixture of the normal sample plus somatic variation.
4. **Quality / thresholds:** Quality scores and filtration thresholds are empirically chosen, not absolute probability calibrations; no single universal minimum quality threshold is mandated.

### Strelka2 (Kim et al. 2018) — successor caller

**URL:** https://www.nature.com/articles/s41592-018-0051-x
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Nature Methods)

**Retrieved how:** WebSearch as above → WebFetch of the Nature article (followed two IDP redirect hops to the final `?error=cookies_not_supported` URL).

**Key Extracted Points:**

1. **Continuous allele-frequency model:** Strelka2 represents continuous allele frequencies for both tumor and normal samples, leveraging the expected genotype structure of the normal.
2. **Somatic classification:** A variant is classified as somatic by comparing VAF in tumor vs normal at each site and requiring the somatic criteria (somatic LOD and somatic VAF) to pass.
3. **Citation:** Kim S, Scheffler K, Halpern AL, et al. (2018). Nature Methods 15:591–594. DOI 10.1038/s41592-018-0051-x.

### GATK Mutect2 (Benjamin et al. 2019) — reference implementation documentation

**URL:** https://raw.githubusercontent.com/broadinstitute/gatk/master/docs/mutect/mutect.tex
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation documentation, Broad Institute GATK)

**Retrieved how:** WebSearch "Mutect2 somatic variant calling tumor normal VAF threshold GATK documentation" and "FilterMutectCalls germline filter ..."; the HTML help page returned HTTP 403, so the canonical LaTeX model source was fetched via raw.githubusercontent.com.

**Key Extracted Points:**

1. **Germline filter model:** Mutect2 compares unnormalized probabilities of (germline het in normal), (germline hom-alt), and (somatic, tumor only). The germline error probability is the normalized sum of the germline possibilities.
2. **Role of matched normal:** **"If we have no matched normal, ℓ_n = 1."** A somatic variant shows strong evidence in the tumor (ℓ_t) while the normal shows essentially no evidence beyond baseline (ℓ_n ≈ 1). Given a matched normal, Mutect2 skips emitting variants clearly present in the germline (i.e. present in the matched normal).
3. **Tumor-normal prior:** From the GATK Mutect2 article (via search snippet), `--af-of-alleles-not-in-resource` defaults to 1e-6 for tumor-normal mode, 5e-8 for tumor-only, 4e-3 for mitochondrial.
4. **Diploid het reference fraction:** FilterMutectCalls distinguishes variants whose allele fractions differ significantly from the diploid heterozygous fraction of **1/2** (used for germline vs somatic in impure samples).

### Yan et al. (2021) — VAF limit of detection

**URL:** https://www.nature.com/articles/s41598-021-91142-1
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Scientific Reports)

**Retrieved how:** WebSearch "somatic variant tumor VAF 5% normal VAF below 1% threshold ..." → WebFetch of the Nature article (followed IDP redirect to the final URL).

**Key Extracted Points:**

1. **Tumor LoD (verbatim):** "WES has a mutation limit of detection (LoD) at variant allele frequencies (VAF) of 5%. Putative mutations called at ≤ 5% VAF are frequently due to sequencing errors, therefore reporting these subclonal mutations incurs risk of significant false positives."
2. **Normal baseline band:** Search-level evidence (npj Precision Oncology, DRAGEN) places normal/baseline variant thresholds in the sub-percent band (≈0.01%–0.5%); a 1% ceiling for "absent in normal" is a conservative noise band consistent with the ref/ref restriction.

---

## Documented Corner Cases and Failure Modes

### From Benjamin et al. 2019 (Mutect2)

1. **Tumor-only mode:** With no matched normal, ℓ_n = 1 — classification relies on the tumor only; the germline filter cannot use a matched normal.
2. **Low tumor purity / impure samples:** Allele fractions drift away from the diploid 1/2; FilterMutectCalls separates somatic from germline by deviation from the 1/2 fraction.

### From Saunders et al. 2012 (Strelka)

1. **LOH / copy-number regions:** Raw somatic probability over-calls in LOH regions; restricting to ref/ref normal genotype mitigates this.

### From Yan et al. 2021

1. **Sub-5% VAF subclonal calls:** Frequently sequencing errors; high false-positive risk → treated as not detected at the standard LoD.

---

## Test Datasets

### Dataset: Synthetic tumor/normal variant panel (derived from the somatic-state definition)

**Source:** Derived directly from Saunders et al. (2012) S = {(f_t, f_n): f_t ≠ f_n} with ref/ref normal, and Yan et al. (2021) 5% tumor LoD. VAF = altReads / totalReads.

| Variant | Tumor alt/total | Tumor VAF f_t | Normal alt/total | Normal VAF f_n | Expected status | Somatic score = max(0, f_t − f_n) |
|---------|-----------------|---------------|------------------|----------------|-----------------|-----------------------------------|
| A (clear somatic) | 25/100 | 0.25 | 0/100 | 0.00 | Somatic | 0.25 |
| B (germline het)  | 48/100 | 0.48 | 50/100 | 0.50 | Germline | 0.00 |
| C (sub-LoD)       | 2/100  | 0.02 | 0/100 | 0.00 | NotDetected | 0.00 |
| D (CHIP-like in normal) | 30/100 | 0.30 | 3/100 | 0.03 | Germline | 0.00 |
| E (tumor-only, no normal coverage) | 20/100 | 0.20 | 0/0 | 0.00 | Somatic | 0.20 |
| F (at tumor threshold) | 5/100 | 0.05 | 0/100 | 0.00 | Somatic | 0.05 |
| G (at normal threshold) | 30/100 | 0.30 | 1/100 | 0.01 | Somatic | 0.29 |

---

## Assumptions

1. **ASSUMPTION: Normal "absent" ceiling = 1% VAF** — Strelka restricts to a ref/ref normal genotype but does not publish a single VAF cutoff (its decision is Bayesian). A concrete rule-based cutoff is required; 1% sits in the sub-percent normal baseline band reported by Yan et al. (2021) and is exposed as a configurable parameter (`normalVafThreshold`), so the correctness-critical decision is parameterized rather than hard-wired to an invented value.
2. **ASSUMPTION: Somatic score = max(0, f_t − f_n)** — No authoritative source publishes this exact closed-form scalar (Mutect2/Strelka use Bayesian LOD scores not retrievable in full). The monotone separation score is a transparent, bounded [0,1] surrogate for "present in tumor, absent in normal"; it is documented as a simplification, not presented as a caller LOD.

---

## Recommendations for Test Coverage

1. **MUST Test:** Clear somatic variant (high tumor VAF, zero normal VAF) classified Somatic with score = f_t − f_n. — Evidence: Saunders et al. (2012) S = {f_t ≠ f_n}, ref/ref normal.
2. **MUST Test:** Germline variant (high tumor VAF, comparable normal VAF) classified Germline. — Evidence: Mutect2 germline filter (Benjamin et al. 2019).
3. **MUST Test:** Sub-5% tumor VAF classified NotDetected. — Evidence: Yan et al. (2021) 5% LoD.
4. **MUST Test:** Tumor-only mode (no normal coverage) classified Somatic. — Evidence: Mutect2 ℓ_n = 1 when no matched normal.
5. **MUST Test:** Boundary at tumor threshold (f_t = 0.05) is present; boundary at normal threshold (f_n = 0.01) is absent. — Evidence: thresholds above.
6. **MUST Test:** FilterGermlineVariants returns only somatic calls. — Evidence: Mutect2 skips germline.
7. **SHOULD Test:** Low-VAF normal contamination (CHIP-like, f_n just above 1%) → Germline. — Rationale: documented contamination case.
8. **COULD Test:** Custom thresholds change classification. — Rationale: parameter validation.

---

## References

1. Saunders CT, Wong WSW, Swamy S, Becq J, Murray LJ, Cheetham RK (2012). Strelka: accurate somatic small-variant calling from sequenced tumor-normal sample pairs. Bioinformatics 28(14):1811–1817. https://doi.org/10.1093/bioinformatics/bts271 (PMID 22581179)
2. Kim S, Scheffler K, Halpern AL, et al. (2018). Strelka2: fast and accurate calling of germline and somatic variants. Nature Methods 15:591–594. https://doi.org/10.1038/s41592-018-0051-x
3. Benjamin D, et al. / Broad Institute (2019). GATK Mutect2 model documentation (mutect.tex). https://raw.githubusercontent.com/broadinstitute/gatk/master/docs/mutect/mutect.tex
4. Yan YH, Chen SX, Cheng LY, et al. (2021). Confirming putative variants at ≤ 5% allele frequency using allele enrichment and Sanger sequencing. Scientific Reports 11:11640. https://doi.org/10.1038/s41598-021-91142-1

---

## 2026-09 review — Mutect2 somatic likelihoods model (sources opened 2026-09-28)

- `https://raw.githubusercontent.com/broadinstitute/gatk/master/docs/mutect/mutect.tex` — §Somatic Likelihoods Model: evidence
  ln P(R|A) ≈ g(α) − g(β) + Σ_ra z̄_ra(ln ℓ_ra − ln z̄_ra); z̄_ra ∝ f̃_a ℓ_ra, ln f̃_a = ψ(β_a) − ψ(Σβ); TLOD = log evidence ratio of
  allele set with vs without the alt. §Germline filter: "If we have no matched normal, ℓ_n = 1".
- `.../tools/walkers/mutect/SomaticLikelihoodsEngine.java` — CONVERGENCE_THRESHOLD 0.001, NEGLIGIBLE_RESPONSIBILITY 1e-10, xLogx cut-off 1e-8, flat initial posterior.
- `.../tools/walkers/mutect/SomaticGenotypingEngine.java` — somaticLogOdds (flat pseudocounts when minAF = 0), diploidAltLogOdds
  (hom-ref − Σ(logsumexp(ref,alt)+ln½)), emission TLOD > emit-lod, genotype only if !hasNormal ∨ NLOD > normal-lod, Q45 mismapping normalisation.
- `.../tools/walkers/mutect/M2ArgumentCollection.java` — DEFAULT_EMISSION_LOG_10_ODDS 3.0, DEFAULT_NORMAL_LOG_10_ODDS 2.2, minAF 0.
- `.../utils/pairhmm/LoglessPairHMM.java` / `PairHMM.java` — mismatch prior qualToErrorProb/TRISTATE_CORRECTION (3.0).
- `.../utils/genotyper/AlleleLikelihoods.java` — normalizeLikelihoods floors each allele at best + cap;
  `.../haplotypecaller/LikelihoodEngineArgumentCollection.java` — phredScaledGlobalReadMismappingRate = 45.

Reference dataset (per-read Python port of the above, scipy digamma/gammaln; exact integral as sanity check):

| ref/alt, Q | TLOD (port) | TLOD (exact ∫, flat prior) | NLOD |
|---|---|---|---|
| 75/25, 30 | 61.54246801643632 | 61.5429 | −56.83 |
| 98/2, 30 | 1.2661263748027127 | 1.2692 | 23.14 |
| 97/3, 30 | 3.2295883512710644 | 3.2318 | 19.658450780899575 |
| 97/3, 20 | 0.3263769386302512 | — | 22.54 |
| 52/48, 30 | 135.9216837925928 | 135.9220 | −136.79 |
| 80/20, 50 (cap active; uncapped 86.81) | 67.26788516897696 | — | −59.89837377162566 |
| 100/0, 30 | −2.0042635241874462 | — | 30.088511009736504 |
| 7/0, 8/0 (normal), 30 | — | — | 2.106195770681557 / 2.407080880778922 |

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-28**: Added Mutect2 somatic likelihoods model evidence (campaign 2026-09 review).
