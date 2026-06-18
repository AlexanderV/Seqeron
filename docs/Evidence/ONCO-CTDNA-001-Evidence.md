# Evidence Artifact: ONCO-CTDNA-001

**Test Unit ID:** ONCO-CTDNA-001
**Algorithm:** ctDNA Analysis (Poisson limit-of-detection, tumor-fraction estimation, mean variant-allele-fraction summarization)
**Date Collected:** 2026-06-15

---

## Online Sources

### Newman et al. 2014 — CAPP-Seq (Nature Medicine), PMC full text

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4016134/ (fetched 2026-06-15; original https://www.nature.com/articles/nm.3519 redirects to an authenticated host)
**Retrieved by:** WebSearch "Newman CAPP-Seq 2014 Nature Medicine ctDNA detection circulating tumor DNA limit of detection Poisson" → WebFetch of the PMC URL.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Limit of detection (mutant allele fraction):** CAPP-Seq "accurately detected defined inputs of NSCLC DNA at fractional abundances between 0.025% and 10%", with "96% specificity for mutant allele fractions down to ~0.02%". (Establishes that ctDNA detection is parameterized by a *mutant allele fraction* d.)
2. **Background error rate:** "mean and median background rates of 0.006% and 0.0003%, respectively" — the analytic noise floor below which an allele fraction cannot be distinguished from error.
3. **Observed ctDNA fractions in patients:** "Fractions of ctDNA detected in plasma by SNV and/or indel reporters ranged from ~0.02% to 3.2%, with a median of ~0.1% in pre-treatment samples" — ctDNA level is summarized as a fraction across SNV/indel reporters.

### USPTO Patent US 11,085,084 B2 — "Identification and use of circulating nucleic acids" (Poisson detection model)

**URL:** https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11085084
**Retrieved by:** WebSearch "ctDNA limit of detection Poisson model minimum mutant molecules required input genome equivalents…" and WebSearch '"1 - exp" … ctDNA detection probability Poisson "genome equivalents" mutant allele fraction'. The detection equations were returned verbatim by the search extraction of this document.
**Accessed:** 2026-06-15
**Authority rank:** 2 (formal specification of the detection model; restates the Avanzini et al. *Science Advances* 2020 shedding-model derivation, which is paywalled at science.org/HTTP 403 and could not be opened this session).

**Key Extracted Points (verbatim):**

1. **Poisson mean:** "the probability of observing a single tumor reporter in cfDNA follows a Poisson distribution with mean λ = n × d", where **n** = number of sequenced genome equivalents and **d** = detection limit (fraction of ctDNA molecules / mutant allele fraction).
2. **Single-reporter detection probability:** "Given one reporter, the probability x of detecting ≥1 ctDNA molecule is equal to 1 − Poisson(λ), which simplifies to: x = 1 − e^(−nd)".
3. **k-reporter detection probability:** "For k independent tumor reporters, the probability p of detecting ≥1 ctDNA molecule is given by: p = 1 − e^(−ndk)".
4. **Low-burden regime:** "In the low-burden regime (λ < 3), detection is governed by Poisson sampling, such that small changes in plasma input or recovery can shift results across the limit of detection."

### Multiplex ddPCR cfDNA quantification — Alcaide et al. 2020 (Scientific Reports), PMC full text

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7387491/ (open access)
**Retrieved by:** WebSearch 'cfDNA "6.6 pg" diploid genome "3.3 pg" haploid genome equivalents per nanogram 303 copies' → WebFetch.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points (verbatim):**

1. **Mass→genome-equivalents conversion:** "Assuming that 1 ng of cfDNA roughly contains 303 haploid genome equivalents, this corresponds to 440/303 = 1.45 ng of cfDNA." → 1 ng ≈ 303 haploid genome equivalents.

### Standardisation of cfDNA measurement — Devonshire et al. 2014 (Anal Bioanal Chem), PMC full text

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4182654/ (open access)
**Retrieved by:** WebSearch (same query as above) → WebFetch.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points (verbatim):**

1. **Mass of one haploid genome:** "one copy is … a single human haploid genome that is calculated as 3.3 pg." → one haploid genome equivalent = 3.3 pg; 1000 pg ÷ 3.3 pg ≈ 303 copies/ng (consistent with Alcaide et al.).

### MRD review — Pessoa et al. 2023 (Mol Cancer / Clin Transl Med), PMC full text

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC10314661/ (open access)
**Retrieved by:** WebSearch "ctDNA limit of detection Poisson model minimum mutant molecules required input genome equivalents…" → WebFetch.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed review citing primaries)

**Key Extracted Points (verbatim):**

1. **Worked molecule count:** "At a VAF of 0.1%, consistent with localized malignancy or post-treatment MRD, this equates to only 15 molecules of tumor DNA" from "~15,000 haploid genome equivalents" (i.e. expected mutant count = n·d = 15,000 × 0.001 = 15; corroborates λ = n·d).

### Tumor fraction ≈ 2 × VAF for a clonal heterozygous SNV (copy-neutral diploid)

**Retrieved by:** WebSearch 'tumor fraction equals two times variant allele frequency clonal heterozygous mutation ctDNA … VAF ichorCNA'. The diploid-heterozygous identity is the same one already cited in this repository for `OncologyAnalyzer.EstimatePurityFromVaf` from CNAqc (Antonello et al. 2024, *Genome Biology* 25:38, https://doi.org/10.1186/s13059-024-03170-5): for m = 1, n_tot = 2 the expected-VAF relation v = m·π / [2(1−π) + π·n_tot] reduces to v = π/2, hence fraction = 2·v.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed; CNAqc relation).

**Key Extracted Points:**

1. **TF = 2 × VAF:** "for a clonal heterozygous mutation in ctDNA, the tumor fraction equals two times the variant allele frequency … because a heterozygous mutation is present on only one of two DNA copies, the observed VAF is half the actual cellular prevalence (tumor fraction)."

---

## Documented Corner Cases and Failure Modes

### From Newman et al. 2014 (CAPP-Seq)

1. **Below background floor:** allele fractions at or under the analytic background (median 0.0003%, mean 0.006%) cannot be reliably distinguished from sequencing error; the validated detection range starts at 0.025%.
2. **Sub-detection ctDNA fraction:** ctDNA fractions reach as low as ~0.02%; a fraction below the assay limit of detection is reported as not detected, not as zero ctDNA.

### From Patent US 11,085,084 / Avanzini shedding model

1. **Poisson-limited low-input regime (λ < 3):** when n·d < 3 the expected mutant-molecule count is small, so detection is stochastic and a true-positive may be missed (false negative) purely from sampling — detection is a probability, not a guarantee.
2. **Zero genome equivalents / zero allele fraction:** λ = 0 ⇒ P(detect) = 1 − e⁰ = 0 (cannot detect with no molecules or no tumor signal).

---

## Test Datasets

### Dataset: CAPP-Seq detection range (Newman et al. 2014)

**Source:** Newman et al. 2014, Nat Med 20(5):548–554, PMC4016134.

| Parameter | Value |
|-----------|-------|
| Validated detection range (mutant allele fraction) | 0.025% – 10% |
| Specificity floor | ~0.02% at 96% specificity |
| Median pre-treatment ctDNA fraction | ~0.1% |
| Mean / median background error | 0.006% / 0.0003% |

### Dataset: Poisson detection worked example (Pessoa et al. 2023)

**Source:** PMC10314661 (corroborates λ = n·d).

| Parameter | Value |
|-----------|-------|
| Genome equivalents n | 15,000 |
| Mutant allele fraction d (VAF) | 0.001 (0.1%) |
| Expected mutant molecules λ = n·d | 15 |
| P(detect ≥1) = 1 − e^(−15) | 0.99999969… (≈ 1) |

### Dataset: Mass→molecule conversion (Devonshire 2014 / Alcaide 2020)

| Parameter | Value |
|-----------|-------|
| Mass of one haploid genome equivalent | 3.3 pg |
| Haploid genome equivalents per ng cfDNA | ≈ 303 (1000 / 3.3) |

---

## Assumptions

1. **ASSUMPTION: Detection-decision rule from the Poisson model.** The patent/Avanzini model gives the *probability* p = 1 − e^(−ndk). The literature does not fix a single universal probability threshold at which a variant is *declared* detected. The implementation therefore (a) returns the exact probability p, and (b) offers a deterministic detectability test against a caller-supplied probability threshold (default 0.95, the assay-standard 95% sensitivity convention) plus the requirement that the expected mutant count λ = n·d·k ≥ 1 (at least one mutant molecule must be expected). Only the probability value itself is non-assumption; the default 0.95 confidence is flagged as an assumption (same 0.95 convention already used by `CalculateVAFConfidenceInterval` in this class). Changing the threshold changes only the boolean detect flag, not the returned probability.

---

## Recommendations for Test Coverage

1. **MUST Test:** `DetectionProbability` returns 1 − e^(−n·d·k); e.g. n=15000, d=0.001, k=1 ⇒ λ=15, p=1−e⁻¹⁵ — Evidence: Patent US11085084 (p = 1 − e^(−ndk)); Pessoa et al. 2023 (n=15000,d=0.001⇒15 molecules).
2. **MUST Test:** `CalculateTumorFraction` = 2 × mean clonal heterozygous VAF (copy-neutral diploid) — Evidence: CNAqc Antonello 2024 (v = π/2 ⇒ TF = 2·v).
3. **MUST Test:** `CalculateMeanVaf` = mean of altReads/totalReads across reporters — Evidence: Newman 2014 (ctDNA fraction summarized across SNV/indel reporters).
4. **MUST Test:** genome-equivalents helper: 1 ng ⇒ ≈303 haploid GE; 3.3 pg ⇒ 1 GE — Evidence: Devonshire 2014 (3.3 pg/haploid); Alcaide 2020 (303/ng).
5. **MUST Test:** λ = 0 (n=0 or d=0) ⇒ p = 0; not detected — Evidence: Poisson model, P=1−e⁰=0.
6. **MUST Test:** below-LoD VAF (e.g. d < 0.025%) detect=false when λ<1 — Evidence: Newman 2014 detection range; λ≥1 rule.
7. **SHOULD Test:** k>1 reporters raises detection probability (p with k=10 > p with k=1 at fixed n,d) — Rationale: multi-reporter integration (Newman 2014).
8. **SHOULD Test:** input validation — negative n, d∉[0,1], k<1, null variant set, VAF>0.5 for tumor fraction — Rationale: documented domain limits.
9. **COULD Test:** monotonicity of detection probability in n and d — Rationale: model is strictly increasing in λ.

---

## References

1. Newman A.M., Bratman S.V., To J., et al. (2014). An ultrasensitive method for quantitating circulating tumor DNA with broad patient coverage. *Nature Medicine* 20(5):548–554. https://doi.org/10.1038/nm.3519 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC4016134/)
2. Avanzini S., Kurtz D.M., Chabon J.J., et al. (2020). A mathematical model of ctDNA shedding predicts tumor detection size. *Science Advances* 6(50):eabc4308. https://doi.org/10.1126/sciadv.abc4308 (paywalled this session; detection equations corroborated verbatim via US Patent 11,085,084).
3. US Patent 11,085,084 B2. Identification and use of circulating nucleic acids. https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11085084
4. Devonshire A.S., Whale A.S., Gutteridge A., et al. (2014). Towards standardisation of cell-free DNA measurement in plasma: controls for extraction efficiency, fragment size bias and quantification. *Anal Bioanal Chem*. https://pmc.ncbi.nlm.nih.gov/articles/PMC4182654/
5. Alcaide M., Cheung M., Hillman J., et al. (2020). Evaluating the quantity, quality and size distribution of cell-free DNA by multiplex droplet digital PCR. *Scientific Reports* 10:12564. https://doi.org/10.1038/s41598-020-69432-x (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC7387491/)
6. Pessoa L.S., et al. (2023). Genomic approaches to cancer and minimal residual disease detection using circulating tumor DNA. https://pmc.ncbi.nlm.nih.gov/articles/PMC10314661/
7. Antonello A., et al. (2024). CNAqc … *Genome Biology* 25:38. https://doi.org/10.1186/s13059-024-03170-5 (clonal heterozygous diploid relation v = π/2 ⇒ fraction = 2·v).
8. Wan J.C.M., Massie C., Garcia-Corbacho J., et al. (2017). Liquid biopsies come of age: towards implementation of circulating tumour DNA. *Nature Reviews Cancer* 17:223–238. https://doi.org/10.1038/nrc.2017.7 (paywalled this session; cited for context only — no value taken from it).

---

## Change History

- **2026-06-15**: Initial documentation (ONCO-CTDNA-001).
