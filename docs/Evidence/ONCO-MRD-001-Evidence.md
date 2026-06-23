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

### INVAR2 reference implementation — exact GLRT / IMAFv2 formulas (nrlab-CRUK/INVAR2)

**URL:** https://github.com/nrlab-CRUK/INVAR2 — files `R/shared/detectionFunctions.R` and `R/4_detection/generalisedLikelihoodRatioTest.R`
**Accessed:** 2026-06-23 (repository cloned with `git clone --depth 1 https://github.com/nrlab-CRUK/INVAR2.git`; source files read directly)
**Authority rank:** 3 (reference implementation of the Wan et al. 2020 INVAR method, by the originating Rosenfeld lab)

**Key Extracted Points (verbatim formulas):**

1. **Per-locus mixture model:** for a locus with tumour allele fraction `AF`, caller-supplied background error rate `e`, and per-sample ctDNA fraction `p`, the probability that a read is mutant is `q = AF·(1−e)·p + (1−AF)·e·p + e·(1−p)` (`calc_log_likelihood`). Equivalently `q = p·g + e·(1−p)` with `g = AF·(1−e) + (1−AF)·e`. Background `e` is the read-error rate under the null (`p = 0`).
2. **Per-locus-mean log-likelihood:** `logL = Σ[ lchoose(R,M) + M·log(q) + (R−M)·log(1−q) ] / length(R)` where `M` = mutant reads, `R` = depth (`calc_log_likelihood`).
3. **EM estimator of ctDNA fraction `p`:** E-step `Z0 = (1−g)·p / ((1−g)·p + (1−e)·(1−p))`, `Z1 = g·p / (g·p + e·(1−p))`; M-step `p = Σ(M·Z1 + (R−M)·Z0)/ΣR`; `initial_p = 0.01`, `iterations = 200` (`estimate_p_EM`).
4. **Generalised likelihood ratio (detection statistic):** `LR = logL(p̂_MLE) − logL(p = 0)` where `p̂_MLE` is the EM estimate (`calc_likelihood_ratio`). Larger `LR` ⇒ stronger ctDNA evidence; for a pure-background sample `p̂ ≈ 0` and `LR ≈ 0`.
5. **AF / signal-to-noise weighting:** the model weights each locus by its tumour AF and background via `q`; high-`AF`, low-`e` loci contribute more to the likelihood gradient, so detection is signal-to-noise weighted (not a flat read pool).
6. **IMAFv2 (background-subtracted, depth-weighted aggregate):** per context, `MEAN_AF.BS = pmax(0, MEAN_AF − BACKGROUND_AF)` then `IMAFV2 = weighted.mean(MEAN_AF.BS, TOTAL_DP)` (`calculateIMAFv2`). This is per-locus/per-context background subtraction followed by depth-weighted aggregation.
7. **Zero-background guard:** `BACKGROUND_AF = ifelse(BACKGROUND_AF > 0, BACKGROUND_AF, 1 / BACKGROUND_DP)` — a zero background is floored to one expected error in the locus depth so `log(e)` stays finite (`doMain`).
8. **Informative-locus filter:** `filter(TUMOUR_AF > 0)` — only loci with positive tumour AF contribute to the likelihood (`doMain`).

### INVAR2 reference implementation — fragment-size weighting, outlier suppression, locus-noise & background estimation (nrlab-CRUK/INVAR2)

**URL:** https://github.com/nrlab-CRUK/INVAR2 — files `R/shared/detectionFunctions.R`, `R/3_outlier_suppression/outlierSuppression.R`, `R/3_outlier_suppression/sizeCharacterisation.R`, `R/4_detection/generalisedLikelihoodRatioTest.R`, `R/1_parse/onTargetErrorRatesAndFilter.R`
**Accessed:** 2026-06-23 (repository cloned with `git clone --depth 1 https://github.com/nrlab-CRUK/INVAR2.git`; the R source files were read directly)
**Authority rank:** 3 (reference implementation of Wan et al. 2020 INVAR, originating Rosenfeld lab)

**Key Extracted Points (verbatim formulas):**

1. **Fragment-size-weighted GLRT** (`calc_log_likelihood_with_RL`): for each molecule (one read; `DP = 1`) with tumour AF `AF`, background `e`, ctDNA fraction `p`, fragment-size probability under the normal size profile `RL_PROB_0` (P0) and the mutant/tumour size profile `RL_PROB_1` (P1), and `g = AF·(1−e)+(1−AF)·e`: the wild-type read likelihood is `L0 = (1−e)·P0·(1−p) + (1−g)·P1·p` and the mutant read likelihood is `L1 = e·P0·(1−p) + g·P1·p`; `logL = Σ[ M·log(L1) + (R−M)·log(L0) ] / length(R)` with `M ∈ {0,1}` the per-molecule mutant indicator. Tumour-derived cfDNA is shorter, so a short fragment has higher P1 and is up-weighted.
2. **Size-weighted EM** (`estimate_p_EM_with_RL`): E-step `Z0 = (1−g)·P1·p / ((1−g)·P1·p + (1−e)·P0·(1−p))`, `Z1 = g·P1·p / (g·P1·p + e·P0·(1−p))`; M-step `p = Σ(M·Z1 + (R−M)·Z0)/ΣR`; `initial_p = 0.01`, `iterations = 200`. The detection statistic is `LR = logL(p̂)−logL(0)` (`calc_likelihood_ratio_with_RL`).
3. **Size profile** (`sizeCharacterisation.R` / `estimate_real_length_probability`): per fragment-length bin, `PROPORTION = COUNT / TOTAL` separately for `MUTANT` (tumour) and non-mutant (normal) reads. The `onlyWeighMutants` fall-back sets non-discriminating reads to a uniform `1/((maxFragmentLength − minFragmentLength)+1)`; defaults `MINIMUM_FRAGMENT_LENGTH = 60`, `MAXIMUM_FRAGMENT_LENGTH = 300`, smoothing bandwidth `0.03`.
4. **Patient-specific outlier suppression** (`repolish`, `outlierSuppression.R`): from loci passing locus-noise & both-strands with `AF ≤ alleleFrequencyThreshold (0.01)`, `MUTATED_READS_PER_LOCI ≤ maximumMutantReads (10)`, `TUMOUR_AF > 0`, estimate the null sample ctDNA fraction `P_ESTIMATE = max(estimate_p_EM(MUTANT, DP, TUMOUR_AF, BACKGROUND_AF), weighted.mean(AF, TUMOUR_AF))`; Bonferroni threshold `P_THRESHOLD = outlierSuppression (0.05) / n_distinct(loci)`; per locus `BINOMIAL_PROB = binom.test(x = MUTATED_READS_PER_LOCI, n = DP, p = P_ESTIMATE, alternative = "greater")$p.value` (`= 1` when `x ≤ 0`); `OUTLIER.PASS = BINOMIAL_PROB > P_THRESHOLD` — a locus is an OUTLIER (removed) when its mutant-read count is a one-sided binomial outlier above the sample estimate.
5. **Locus-noise filtering & background-error estimation** (`createLociErrorRateTable`, `onTargetErrorRatesAndFilter.R`): over control (`!PATIENT_SPECIFIC & CASE_OR_CONTROL == 'case'`) samples per locus, `BACKGROUND_AF = sum(ALT_F + ALT_R) / sum(DP)`, `N_SAMPLES` distinct control samples, `N_SAMPLES_WITH_SIGNAL` controls with `ALT_F + ALT_R > 0`; `LOCUS_NOISE.PASS = (N_SAMPLES_WITH_SIGNAL / N_SAMPLES) < proportion_of_controls (0.1) AND BACKGROUND_AF < max_background_mean_AF (0.01)`. The per-locus background error model `e` is thus ESTIMATED from control plasma (not only caller-supplied).
6. **Both-strands filter** (`onTargetErrorRatesAndFilter.R` main): `BOTH_STRANDS.PASS = ALT_F > 0 & ALT_R > 0 | AF == 0`.

---

## Documented Corner Cases and Failure Modes

### From the Signatera white paper / Reinert 2019

1. **Exactly 1 variant detected:** below the ≥2 threshold ⇒ MRD-negative (single-variant signal is treated as not sufficient for a positive call).
2. **< 8 variants tracked:** sensitivity at ≤0.1% VAF is compromised; the assay is designed around 16 markers. (Affects sensitivity, not the calling rule itself.)
3. **No tumor markers / empty panel:** nothing to interrogate ⇒ undefined input.

### From Wan 2020 (INVAR)

1. **Per-locus background error:** real ctDNA signal must exceed background; loci with no supporting alt reads contribute 0 to detected count.

### From INVAR2 (GLRT / IMAFv2)

1. **Pure-background sample:** mutant reads occurring only at the background rate ⇒ EM `p̂ ≈ 0` and `LR ≈ 0` ⇒ not detected. Background subtraction removes pure noise.
2. **Zero background, zero signal:** `q → 0` at `p = 0`; the implementation clamps `q` to `(0,1)` (and INVAR floors `e` to `1/depth`) so logs are finite.
3. **No informative locus (all tumour AF = 0):** nothing to estimate ⇒ undefined input (INVAR `filter(TUMOUR_AF > 0)` empties the table).
4. **Flat size profile:** when the mutant and normal fragment-size profiles are identical (P1 == P0) the size factor cancels and the with-RL GLRT reduces to the no-size GLRT — size weighting adds discrimination only when the distributions differ.
5. **Outlier `x ≤ 0`:** a locus with no mutant reads has binomial tail probability 1 (`ifelse(x <= 0, 1, ...)`) and is never an outlier.
6. **Empty control panel:** background cannot be estimated with no control observation ⇒ undefined input.
7. **Both-strands:** `AF == 0` (no alt reads) passes vacuously; alt on a single strand only is strand-biased ⇒ fails.

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
| n=1000, f=0.001, m=16 | λ = 16 | 0.9999998874648253 |

### Dataset: INVAR GLRT synthetic recovery (controlled injection)

**Source:** Derived from the INVAR2 formulas above (`calc_log_likelihood`, `estimate_p_EM`, `calc_likelihood_ratio`), evaluated on synthetic loci. n = 50 loci, depth R = 1000/locus, tumour AF = 0.4, background e = 0.001; mutant reads `M = round(q·R)` with `q = inj·g + e·(1−inj)`, `g = AF·(1−e)+(1−AF)·e`. Values are computed independently of the C# implementation (Python reference of the same equations).

| Injected ctDNA fraction (inj) | Mutant reads/locus M | Estimated p̂ (EM) | Likelihood ratio LR |
|-------------------------------|----------------------|-------------------|---------------------|
| 0.000 (pure background)       | 1                    | ≈ 3.3e-5 (≈ 0)    | ≈ −0.0001 (≈ 0)     |
| 0.005                         | 5                    | ≈ 0.00501         | ≈ 1.30              |
| 0.010                         | 5                    | ≈ 0.01002         | ≈ 4.06              |
| 0.020                         | 9                    | ≈ 0.02004         | ≈ 11.81             |
| 0.050                         | 21                   | ≈ 0.0501          | ≈ 44.14             |

- **Background subtraction:** pure-background (inj=0) ⇒ p̂ ≈ 0, LR ≈ 0 ⇒ not detected.
- **Recovery:** p̂ tracks the injected fraction to ~3 significant figures.
- **Monotonicity:** LR strictly increases with injected signal.

### Dataset: AF weighting boosts sensitivity (low-signal mixture)

**Source:** Same INVAR2 GLRT formulas. N = 40 loci (20 with tumour AF = 0.5, 20 with AF = 0.05), depth 2000, e = 0.002, injected ctDNA fraction 0.008; mutant reads per locus from each locus's true AF. "Weighted" uses each locus's true AF; "unweighted" replaces every AF by the panel mean AF (flat pooling, no SNR weighting).

| Model | Likelihood ratio LR |
|-------|---------------------|
| AF-weighted (true per-locus AF) | ≈ 2.66 |
| Unweighted (flat mean AF)       | ≈ 1.91 |

AF-weighting yields a strictly larger detection statistic than flat pooling on the same data ⇒ higher sensitivity at low signal.

### Dataset: fragment-size weighting boosts sensitivity (short tumour fragments)

**Source:** INVAR2 `calc_likelihood_ratio_with_RL` / `estimate_p_EM_with_RL` / `calc_log_likelihood_with_RL` (R/shared/detectionFunctions.R), evaluated on synthetic per-molecule data. 200 molecules at one synthetic panel, tumour AF = 0.4, background e = 0.002. Two fragment-size bins: SHORT (len 120) and LONG (len 170). Mutant (tumour) molecules: 10 short + 2 long; wild-type molecules: 38 short + 150 long. Tumour size profile P1 = {short 0.8, long 0.2}; normal size profile P0 = {short 0.2, long 0.8} (tumour-derived cfDNA short-enriched). Values computed independently of the C# implementation (Python port of the exact INVAR2 with-RL equations).

| Model | Likelihood ratio LR | EM p̂ |
|-------|---------------------|-------|
| Size-weighted (with-RL) | 0.19691792427890276 | 0.12042621132507245 |
| No-size (flat fragment weight) | 0.14779084591516753 | 0.1455823293172657 |

- **Size weighting:** when tumour fragments are shorter than background, the with-RL GLRT statistic exceeds the no-size statistic on the same molecules ⇒ higher sensitivity.
- **Flat profile sanity:** when P1 == P0 (no size discrimination) the with-RL GLRT reduces exactly to the no-size GLRT.

### Dataset: patient-specific outlier suppression (repolish)

**Source:** INVAR2 `repolish` (R/3_outlier_suppression/outlierSuppression.R), evaluated on synthetic loci. 10 loci, depth 1000, tumour AF = 0.4, background e = 0.001. Nine background-only loci with 1 mutant read each; one PLANTED outlier locus with 50 mutant reads. `P_ESTIMATE = max(EM, weighted.mean(AF, TUMOUR_AF)) = 0.001`; `P_THRESHOLD = 0.05 / 10 = 0.005`. One-sided binomial tail `P(X ≥ x | n = 1000, p = 0.001)`:

| Locus | Mutant reads | Binomial tail P(X ≥ x) | OUTLIER (removed)? |
|-------|--------------|------------------------|--------------------|
| clean (×9) | 1 | 0.6323045752290356 | no (tail > 0.005) |
| planted | 50 | 3.7264670792676273e-66 | yes (tail ≤ 0.005) |

- **Outlier removal recovers the true signal:** after dropping the planted locus the residual is pure background (every clean locus VAF 0.001 == background) ⇒ IMAFv2 = 0; with the outlier present IMAFv2 is inflated.

### Dataset: control-derived background-error estimation & locus-noise filter (createLociErrorRateTable)

**Source:** INVAR2 `createLociErrorRateTable` (R/1_parse/onTargetErrorRatesAndFilter.R). 20 control samples each depth 1000; `controlProportion = 0.1`, `maxBackgroundAlleleFrequency = 0.01`.

| Locus | Controls with signal | Pooled BACKGROUND_AF | LOCUS_NOISE.PASS |
|-------|----------------------|----------------------|------------------|
| clean (signal in 1, 1 alt) | 1 / 20 = 0.05 | 1/20000 = 5e-5 | pass (0.05 < 0.1 and 5e-5 < 0.01) |
| recurrent (signal in 5) | 5 / 20 = 0.25 | 5/20000 = 2.5e-4 | fail (0.25 ≥ 0.1) |
| high background (1 sample, 250 alt) | 1 / 20 = 0.05 | 250/20000 = 0.0125 | fail (0.0125 ≥ 0.01) |
| clean, injected error 0.002 (all 20 controls, 2 alt) | 20 / 20 | 40/20000 = 0.002 | (estimate recovers injected 0.002) |

- **Background recovery:** the pooled control allele fraction equals the injected per-locus error rate (0.002).

---

## Assumptions

1. **ASSUMPTION: per-variant "detected" criterion** — A tracked variant is counted as *detected* in plasma when it has at least one supporting alternate read (alt reads ≥ a minimum supporting-read count, default 1), i.e. signal above the trivial-zero background. The cited sources define positivity at the *panel* level (≥2 variants) and require per-locus signal above background, but do not publish an exact universal per-locus read-count cutoff (it is instrument/error-model specific, e.g. INVAR's trinucleotide GLRT). The default ≥1 alt read is the minimal, source-consistent presence rule and is a configurable parameter; it does not change the panel-level ≥2 calling rule. Correctness-affecting only for the per-variant flag, which is exposed as a tunable threshold.
2. **ASSUMPTION: empirical size proportion in place of KDE.** The INVAR2 size likelihood `RL_PROB_0/RL_PROB_1` is computed by a weighted kernel density estimate (`estimate_real_length_probability`, bandwidth `adjust = 0.03`) over the per-length counts. The C# `FragmentSizeProfile` uses the discrete empirical proportion `COUNT/TOTAL` per length bin (the un-smoothed estimator the KDE converges to), with the same uniform `1/((max−min)+1)` fall-back as INVAR2's `onlyWeighMutants`. The mixture/EM/GLRT-with-RL equations are ported verbatim; only the kernel SMOOTHING of the size histogram is not reproduced. Correctness-affecting only for the exact size-weight on sparsely-populated length bins; the size-likelihood-ratio direction (short fragments up-weighted) and the with-RL GLRT are exact.

---

## Recommendations for Test Coverage

1. **MUST Test:** DetectMRD returns MRD-positive iff ≥2 of the tracked variants are detected (cases 0,1,2,3 detected). — Evidence: PMC9265001 Table 1 / Reinert 2019.
2. **MUST Test:** DetectMRD reports DetectedVariantCount and TrackedVariantCount correctly. — Evidence: Signatera white paper (16 tracked).
3. **MUST Test:** panel-level Poisson detection probability p = 1 − e^(−n·f·m) matches the existing primitive for the panel size m. — Evidence: white paper Figure 2.
4. **MUST Test:** IMAF = depth-weighted (Σ alt / Σ total over loci) mean plasma VAF across tracked loci. — Evidence: Wan 2020 (IMAF).
5. **SHOULD Test:** custom positivity threshold (e.g. ≥1, ≥3) shifts the call. — Rationale: threshold is parameterized.
6. **SHOULD Test:** TrackVariantsOverTime yields per-timepoint MRD status and flags first positive timepoint. — Rationale: longitudinal method in scope.
7. **COULD Test:** null/empty panel and invalid threshold raise the documented exceptions. — Rationale: input validation.
8. **MUST Test:** `IntegratedMutantAlleleFractionV2` = depth-weighted mean of `max(0, locusVAF − background)`; a locus whose VAF ≤ background contributes 0. — Evidence: INVAR2 `calculateIMAFv2`.
9. **MUST Test:** `EstimateInvarSignal` on pure-background loci ⇒ p̂ ≈ 0, LR ≈ 0, not detected. — Evidence: INVAR GLRT synthetic-recovery dataset (inj=0).
10. **MUST Test:** `EstimateInvarSignal` recovers an injected ctDNA fraction within tolerance and reports detected. — Evidence: synthetic-recovery dataset (inj=0.01, 0.02, 0.05).
11. **MUST Test:** AF-weighted LR ≥ unweighted (flat-AF) LR on the low-signal mixture. — Evidence: AF-weighting dataset.
12. **MUST Test:** LR is monotone non-decreasing in injected signal. — Evidence: synthetic-recovery dataset.
13. **SHOULD Test:** detectionThreshold gates the detection call (high threshold ⇒ not detected even with weak signal). — Rationale: parameterised specificity knob.
14. **SHOULD Test:** out-of-range tumour AF / background and empty informative panel raise documented exceptions. — Rationale: input validation.

---

## References

1. Reinert T, Henriksen TV, Christensen E, et al. (2019). Analysis of Plasma Cell-Free DNA by Ultradeep Sequencing in Patients With Stages I to III Colorectal Cancer. *JAMA Oncology* 5(8):1124–1131. https://pubmed.ncbi.nlm.nih.gov/31070691/ (DOI: 10.1001/jamaoncol.2019.0528)
2. Natera Inc. (2020). A personalized, tumor-informed approach to detect molecular residual disease with high sensitivity and specificity (Signatera analytical-validation white paper). https://www.natera.com/wp-content/uploads/2020/11/Oncology-Clinical-A-personalized-tumor-informed-approach-to-detect-molecular-residual-disease-SGN_SR_WP.pdf
3. Wan JCM, Heider K, Gale D, et al. (2020). ctDNA monitoring using patient-specific sequencing and integration of variant reads. *Science Translational Medicine* 12(548):eaaz8084. https://www.science.org/doi/10.1126/scitranslmed.aaz8084 (DOI: 10.1126/scitranslmed.aaz8084)
4. Tie J, et al. / review: Tumor-informed ctDNA MRD assessment in colorectal cancer (quotes the Reinert/Signatera 16-SNV, ≥2-positive rule, Table 1). PMC9265001. https://pmc.ncbi.nlm.nih.gov/articles/PMC9265001/
5. Avanzini S, et al. (2020). A mathematical model of ctDNA shedding predicts tumor detection size. *Science Advances* 6(50):eabc4308 (Poisson detection model p = 1 − e^(−λ)). DOI: 10.1126/sciadv.abc4308
6. Rosenfeld lab (nrlab-CRUK). INVAR2 — restructured INVAR pipeline (reference implementation of Wan et al. 2020). `R/shared/detectionFunctions.R`, `R/4_detection/generalisedLikelihoodRatioTest.R`. https://github.com/nrlab-CRUK/INVAR2
7. Lanczos C. (1964). A precision approximation of the gamma function. *J. SIAM Numer. Anal.* 1(1):86–96 (log-gamma used for the binomial coefficient in the GLRT). DOI: 10.1137/0701008

---

## Change History

- **2026-06-15**: Initial documentation (ONCO-MRD-001).
- **2026-06-23**: Added INVAR2 reference-implementation source and the exact GLRT/EM/IMAFv2 formulas; added synthetic GLRT-recovery and AF-weighting datasets and corner cases; extended coverage recommendations for the new background-subtracted, AF-weighted estimator (`EstimateInvarSignal`, `IntegratedMutantAlleleFractionV2`).
- **2026-06-23**: Closed the residual — ported INVAR2 fragment-size weighting (`calc_likelihood_ratio_with_RL` → `EstimateInvarSignalWithSize`/`FragmentSizeProfile`/`InvarMolecule`), patient-specific outlier suppression (`repolish` → `SuppressOutlierLoci`), and locus-noise filtering + control-derived background-error estimation (`createLociErrorRateTable` → `EstimateLocusBackground`/`PassesBothStrandsFilter`). Added the size-weighting, outlier-suppression and background-estimation synthetic datasets with INVAR2-derived expected values.
