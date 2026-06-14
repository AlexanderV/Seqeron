# Evidence Artifact: ONCO-MRD-001

**Test Unit ID:** ONCO-MRD-001
**Algorithm:** Minimal (Molecular) Residual Disease Detection — tumor-informed panel-level ctDNA MRD call
**Date Collected:** 2026-06-15

---

## Online Sources

### Signatera analytical-validation white paper (Natera)

**URL:** https://www.natera.com/wp-content/uploads/2020/11/Oncology-Clinical-A-personalized-tumor-informed-approach-to-detect-molecular-residual-disease-SGN_SR_WP.pdf
**Accessed:** 2026-06-15 (PDF fetched via WebFetch; binary saved locally and parsed with `pdftotext`)
**Authority rank:** 3 (reference implementation / assay technical specification)

**Key Extracted Points:**

1. **Panel size — 16 SNVs:** "a bespoke assay of 16 tumor-specific, clonal, somatic variants are generated for each patient." / "a proprietary algorithm is used to select a set of 16 somatic SNVs for multiplex PCR primer design." (verbatim from the white paper, parsed text).
2. **Tumor-informed workflow:** somatic clonal variants identified by whole-exome sequencing of primary tumor + matched normal (buffy coat); top 16 selected by clonality, detectability, and frequency; tracked longitudinally in plasma cfDNA by patient-specific 16-plex PCR + ultra-deep NGS.
3. **Poisson detection limit (verbatim, Figure 2 "Theoretical Limit of ctDNA Detection"):** `p = 1 - e^(-nfm)` where `p = Probability of ≥ 1 ctDNA molecule detection`, `n = haploid genome equivalents`, `f = ctDNA % (VAF)`, `m = tumor somatic mutations`. This is the same Poisson model as ONCO-CTDNA-001 `CtDnaDetectionProbability` (p = 1 − e^(−n·d·k)) with m (tracked mutations) = k (reporters).
4. **Sensitivity scales with panel size:** "The clinical sensitivity of MRD detection (range: 0.01%-0.1% VAF) is dependent on the number of somatic SNVs tracked … MRD detection of ≤0.1% VAF is compromised when targeting ≤8 clonal mutations … more than 8 [are needed]." Limit of detection in VAF = 0.01%; analytical specificity >99.5%.

### Tumor-informed ctDNA review (PMC9265001) quoting the Reinert/Signatera positivity definition

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9265001/
**Accessed:** 2026-06-15 (WebFetch)
**Authority rank:** 4 (peer-reviewed review citing primary Reinert 2019; used for the positivity rule it quotes)

**Key Extracted Points:**

1. **Positivity rule (verbatim, the review's Table 1):** "Sixteen patient-specific, somatic SNVs are selected for each patient based on the whole-exome sequencing of the tumor for interrogation in the cfDNA. Plasma samples with at least two tumor-specific SNVs are defined as ctDNA-positive."
2. **Threshold:** a plasma sample is ctDNA-positive (MRD-positive) when **≥ 2 of the 16 tracked variants are detected**; fewer than 2 detected ⇒ MRD-negative.

### Reinert et al. 2019 — primary clinical study (PubMed)

**URL:** https://pubmed.ncbi.nlm.nih.gov/31070691/
**Accessed:** 2026-06-15 (WebFetch; abstract retrieved — methodological threshold not in abstract, taken from sources above which quote the study's methods)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Study:** Reinert T, et al. (2019) "Analysis of Plasma Cell-Free DNA by Ultradeep Sequencing in Patients With Stages I to III Colorectal Cancer," JAMA Oncology 5(8):1124–1131. Tumor-informed personalized multiplex-PCR MRD assay (the assay later commercialized as Signatera).
2. **Clinical signal:** post-surgery ctDNA-positive patients ~7× more likely to relapse (HR 7.2; 95% CI 2.7–19.0; p < 0.001) — confirms the post-treatment MRD-positive call is the clinically meaningful output.

### Wan et al. 2020 — INVAR (integration of variant reads), Science Translational Medicine

**URL:** https://www.science.org/doi/10.1126/scitranslmed.aaz8084 (paywalled, 403 on fetch); corroborating descriptions via https://www.eurekalert.org/news-releases/846258 and search-result abstracts
**Accessed:** 2026-06-15 (WebFetch of EurekAlert; WebSearch for abstract content)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Integration principle:** INVAR ("INtegration of VAriant Reads") improves ctDNA sensitivity by analyzing hundreds–thousands of tumor-informed mutations and integrating signal across all loci, reaching ~1 mutant molecule per 100,000 (and parts-per-million in favorable cases).
2. **Integrated mutant allele fraction (IMAF):** ctDNA level is summarized as a depth-weighted average of the per-locus mutant fractions across the tracked loci, corrected for per-locus background error. This motivates the depth-weighted (read-pooled) panel VAF reported alongside the MRD call.

---

## Documented Corner Cases and Failure Modes

### From the Signatera white paper / Reinert 2019

1. **Exactly 1 variant detected:** below the ≥2 threshold ⇒ MRD-negative (single-variant signal is treated as not sufficient for a positive call).
2. **< 8 variants tracked:** sensitivity at ≤0.1% VAF is compromised; the assay is designed around 16 markers. (Affects sensitivity, not the calling rule itself.)
3. **No tumor markers / empty panel:** nothing to interrogate ⇒ undefined input.

### From Wan 2020 (INVAR)

1. **Per-locus background error:** real ctDNA signal must exceed background; loci with no supporting alt reads contribute 0 to detected count.

---

## Test Datasets

### Dataset: Signatera positivity rule (canonical worked cases)

**Source:** PMC9265001 Table 1 (quoting Reinert 2019); Natera white paper.

| Parameter | Value |
|-----------|-------|
| Panel size (tracked SNVs) | up to 16 |
| MRD-positive threshold | ≥ 2 tracked variants detected in plasma |
| 2 of 16 detected | MRD-positive |
| 1 of 16 detected | MRD-negative |
| 0 of 16 detected | MRD-negative |
| 3 of 16 detected | MRD-positive |

### Dataset: Poisson panel detection probability (reuses ONCO-CTDNA-001)

**Source:** Natera white paper Figure 2 (`p = 1 - e^(-nfm)`); Avanzini et al. 2020.

| Parameter | Value | p = 1 − e^(−n·f·m) |
|-----------|-------|--------------------|
| n=1000, f=0.001, m=1 | λ = 1 | 0.6321205588285577 |
| n=1000, f=0.001, m=16 | λ = 16 | 0.9999998874648379 |

---

## Assumptions

1. **ASSUMPTION: per-variant "detected" criterion** — A tracked variant is counted as *detected* in plasma when it has at least one supporting alternate read (alt reads ≥ a minimum supporting-read count, default 1), i.e. signal above the trivial-zero background. The cited sources define positivity at the *panel* level (≥2 variants) and require per-locus signal above background, but do not publish an exact universal per-locus read-count cutoff (it is instrument/error-model specific, e.g. INVAR's trinucleotide GLRT). The default ≥1 alt read is the minimal, source-consistent presence rule and is a configurable parameter; it does not change the panel-level ≥2 calling rule. Correctness-affecting only for the per-variant flag, which is exposed as a tunable threshold.

---

## Recommendations for Test Coverage

1. **MUST Test:** DetectMRD returns MRD-positive iff ≥2 of the tracked variants are detected (cases 0,1,2,3 detected). — Evidence: PMC9265001 Table 1 / Reinert 2019.
2. **MUST Test:** DetectMRD reports DetectedVariantCount and TrackedVariantCount correctly. — Evidence: Signatera white paper (16 tracked).
3. **MUST Test:** panel-level Poisson detection probability p = 1 − e^(−n·f·m) matches the existing primitive for the panel size m. — Evidence: white paper Figure 2.
4. **MUST Test:** IMAF = depth-weighted (Σ alt / Σ total over loci) mean plasma VAF across tracked loci. — Evidence: Wan 2020 (IMAF).
5. **SHOULD Test:** custom positivity threshold (e.g. ≥1, ≥3) shifts the call. — Rationale: threshold is parameterized.
6. **SHOULD Test:** TrackVariantsOverTime yields per-timepoint MRD status and flags first positive timepoint. — Rationale: longitudinal method in scope.
7. **COULD Test:** null/empty panel and invalid threshold raise the documented exceptions. — Rationale: input validation.

---

## References

1. Reinert T, Henriksen TV, Christensen E, et al. (2019). Analysis of Plasma Cell-Free DNA by Ultradeep Sequencing in Patients With Stages I to III Colorectal Cancer. *JAMA Oncology* 5(8):1124–1131. https://pubmed.ncbi.nlm.nih.gov/31070691/ (DOI: 10.1001/jamaoncol.2019.0528)
2. Natera Inc. (2020). A personalized, tumor-informed approach to detect molecular residual disease with high sensitivity and specificity (Signatera analytical-validation white paper). https://www.natera.com/wp-content/uploads/2020/11/Oncology-Clinical-A-personalized-tumor-informed-approach-to-detect-molecular-residual-disease-SGN_SR_WP.pdf
3. Wan JCM, Heider K, Gale D, et al. (2020). ctDNA monitoring using patient-specific sequencing and integration of variant reads. *Science Translational Medicine* 12(548):eaaz8084. https://www.science.org/doi/10.1126/scitranslmed.aaz8084 (DOI: 10.1126/scitranslmed.aaz8084)
4. Tie J, et al. / review: Tumor-informed ctDNA MRD assessment in colorectal cancer (quotes the Reinert/Signatera 16-SNV, ≥2-positive rule, Table 1). PMC9265001. https://pmc.ncbi.nlm.nih.gov/articles/PMC9265001/
5. Avanzini S, et al. (2020). A mathematical model of ctDNA shedding predicts tumor detection size. *Science Advances* 6(50):eabc4308 (Poisson detection model p = 1 − e^(−λ)). DOI: 10.1126/sciadv.abc4308

---

## Change History

- **2026-06-15**: Initial documentation (ONCO-MRD-001).
