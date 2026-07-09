---
type: concept
title: "Tumor-informed MRD detection (Signatera ≥2/16 + INVAR GLRT)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-MRD-001-Evidence.md
source_commit: 13a9268a9bcc252ce0f9ce183815411f74f1ff5b
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-mrd-001-evidence
      evidence: "Test Unit ID: ONCO-MRD-001; Algorithm: Minimal (Molecular) Residual Disease Detection — tumor-informed panel-level ctDNA MRD call"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:ctdna-detection-and-tumor-fraction
      source: onco-mrd-001-evidence
      evidence: "Panel Poisson detection probability p = 1 − e^(−nfm) is the same model as ONCO-CTDNA-001 CtDnaDetectionProbability p = 1 − e^(−n·d·k) with m tracked mutations = k reporters; the MRD unit reuses the ONCO-CTDNA-001 primitive for the panel size and the across-reporter mean VAF on the same cfDNA input."
      confidence: high
      status: current
---

# Tumor-informed MRD detection (Signatera ≥2/16 + INVAR GLRT)

The **calling / statistical-detection** layer of the Oncology family for **minimal (molecular)
residual disease (MRD)** from a liquid biopsy: given a patient-specific panel of tumor-informed
somatic variants tracked in plasma cell-free DNA, return a **positive/negative MRD verdict** and a
quantified ctDNA level. It sits on top of the Poisson quantification primitive
[[ctdna-detection-and-tumor-fraction]] (ONCO-CTDNA-001) — same cfDNA input, but a **multi-variant
verdict** rather than a single-reporter probability. Two complementary detection engines: the
**Signatera-style panel positivity rule** (≥2 of 16 tracked variants) and the **INVAR
generalised-likelihood-ratio (GLRT)** statistical detector that integrates weak signal across all
tracked loci. Validated under test unit **ONCO-MRD-001**; the literature-traced record is
[[onco-mrd-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade
([[scientific-rigor|research-grade]]), **not for clinical or diagnostic use**.

## Why a distinct unit from ONCO-CTDNA-001

[[ctdna-detection-and-tumor-fraction]] answers *"is a single tumor molecule detectable, and how much
tumor is there?"* via the Poisson probability `p = 1 − e^(−n·d·k)`, tumor fraction, and mean VAF.
This unit answers *"is the patient MRD-positive?"* — a **multi-variant calling** problem with its own
machinery: the tumor-informed 16-SNV workflow, the **≥2-detected positivity rule**, longitudinal
monitoring, and the INVAR **background-subtracted, AF-/size-weighted GLRT** with EM ctDNA-fraction
estimation, outlier suppression, and control-derived background modeling. It **reuses** the Poisson
primitive (`p = 1 − e^(−nfm)`, m tracked mutations = k reporters) but is otherwise a separate algorithm.

## 1. Tumor-informed workflow and the Signatera positivity rule (`DetectMRD`)

The Signatera design (Reinert 2019 / Natera white paper): identify clonal somatic variants by
whole-exome sequencing of the primary tumor + matched normal (buffy coat); select **up to 16** by
clonality, detectability, and frequency; track them longitudinally in plasma cfDNA by patient-specific
16-plex PCR + ultra-deep NGS.

**Positivity rule (verbatim, PMC9265001 Table 1 / Reinert 2019):** a plasma sample is **ctDNA-positive
(MRD-positive) iff ≥ 2 of the 16 tracked variants are detected**; fewer than 2 detected ⇒ MRD-negative.
`DetectMRD` implements exactly this (threshold parameterized) and reports `DetectedVariantCount` /
`TrackedVariantCount`.

| Detected of 16 | MRD call |
|----------------|----------|
| 0 | negative |
| 1 | negative (single-variant signal insufficient) |
| 2 | **positive** |
| 3 | **positive** |

- **Sensitivity scales with panel size:** MRD detection at ≤0.1% VAF is compromised when tracking ≤8
  clonal mutations (>8 needed); the assay is built around 16 markers. **LoD 0.01% VAF; analytic
  specificity >99.5%.** (Affects sensitivity, not the calling rule.)
- **Clinical signal (Reinert 2019):** post-surgery ctDNA-positive patients are **~7× more likely to
  relapse** (HR 7.2; 95% CI 2.7–19.0; p < 0.001) — the post-treatment positive call is the meaningful output.

**Longitudinal monitoring (`TrackVariantsOverTime`):** applies the positivity rule at each timepoint,
yielding per-timepoint MRD status and flagging the **first positive timepoint**.

## 2. Panel Poisson detection probability (reuses ONCO-CTDNA-001)

Natera white paper Figure 2: `p = 1 − e^(−nfm)` (n = haploid genome equivalents, f = ctDNA VAF,
m = tumor somatic mutations tracked) — the **same Poisson model** as
[[ctdna-detection-and-tumor-fraction]]'s `p = 1 − e^(−n·d·k)` with m (tracked mutations) = k
(reporters). Worked: n=1000, f=0.001, m=1 (λ=1) → 0.6321; m=16 (λ=16) → 0.99999989. Integrating more
tracked variants raises panel sensitivity.

## 3. IMAF and IMAFv2 — quantifying the ctDNA level

- **IMAF (Wan 2020):** the ctDNA level summarized as a **depth-weighted** (read-pooled) mean of the
  per-locus plasma VAFs across the tracked loci (Σ alt / Σ total).
- **`IntegratedMutantAlleleFractionV2` (INVAR2 `calculateIMAFv2`):** per-context **background-subtracted,
  depth-weighted** aggregate — `MEAN_AF.BS = max(0, MEAN_AF − BACKGROUND_AF)`, then
  `IMAFV2 = weighted.mean(MEAN_AF.BS, TOTAL_DP)`. A locus whose VAF ≤ its background contributes **0**,
  so pure noise is removed.

## 4. INVAR generalised-likelihood-ratio detection (`EstimateInvarSignal`)

INVAR (Wan 2020) integrates weak tumor signal across hundreds–thousands of tumor-informed loci,
reaching ~1 mutant molecule per 100,000. The INVAR2 reference implementation gives the exact model:

- **Per-locus mixture:** the probability a read is mutant is `q = p·g + e·(1−p)` where
  `g = AF·(1−e) + (1−AF)·e`, `AF` = tumour allele fraction, `e` = per-locus background error rate,
  `p` = per-sample ctDNA fraction. Background `e` is the read-error rate under the null (`p = 0`).
- **EM estimator of `p`** (`estimate_p_EM`): E-step `Z0 = (1−g)p / ((1−g)p + (1−e)(1−p))`,
  `Z1 = g·p / (g·p + e(1−p))`; M-step `p = Σ(M·Z1 + (R−M)·Z0) / ΣR`; `initial_p = 0.01`, 200 iterations.
- **Detection statistic:** `LR = logL(p̂_MLE) − logL(p = 0)`. Larger LR ⇒ stronger ctDNA evidence; a
  pure-background sample yields `p̂ ≈ 0` and `LR ≈ 0` (not detected). A caller-supplied `detectionThreshold`
  gates the call.
- **Signal-to-noise (AF) weighting:** high-`AF`, low-`e` loci contribute more to the likelihood gradient,
  so detection is SNR-weighted, not a flat read pool — strictly more sensitive at low signal.

Synthetic-recovery oracles (computed independently of the implementation): pure background (inj=0) →
p̂≈0, LR≈0; inj=0.01 → p̂≈0.01002, LR≈4.06; inj=0.02 → LR≈11.81; inj=0.05 → LR≈44.14 — **monotone in
signal**. AF-weighted LR ≈ 2.66 > flat-mean-AF LR ≈ 1.91 on the same low-signal mixture.

## 5. Fragment-size weighting (`EstimateInvarSignalWithSize` / `FragmentSizeProfile` / `InvarMolecule`)

Tumor-derived cfDNA fragments are **shorter** than normal, so fragment length is an extra signal. The
with-RL GLRT scores each molecule under a tumour size profile `P1` and a normal profile `P0`:
`L0 = (1−e)·P0·(1−p) + (1−g)·P1·p` (wild-type read), `L1 = e·P0·(1−p) + g·P1·p` (mutant read); a short
fragment has higher `P1` and is up-weighted. The size-weighted EM mirrors §4 with `P0`/`P1` factors.

- **Oracle:** with-RL LR ≈ 0.1969 > no-size LR ≈ 0.1478 when tumour fragments are short-enriched.
- **Flat-profile sanity:** when `P1 == P0` the size factor cancels and the with-RL GLRT reduces
  **exactly** to the no-size GLRT — size weighting adds discrimination only when the distributions differ.

The size profiles are estimated per fragment-length bin as `COUNT/TOTAL` (**default**, discrete
empirical). `FragmentSizeProfile.FromKernelDensity(...)` provides an **opt-in** Gaussian-KDE-smoothed
weight (Silverman 1986 eq. 2.2a, `wᵢ = countᵢ/Σcount`, bandwidth by `bw.nrd0`
`0.9·min(σ̂, IQR/1.34)·n^(−1/5)` scaled by an `adjust` multiplier, integrated analytically over each
integer bin `[L−0.5, L+0.5]` via the standard-normal CDF `Φ`), matching INVAR2's
`estimate_real_length_probability` (`density()` + `integrate`, `adjust = 0.03`). All mixture/EM/GLRT
equations are unchanged — the KDE only changes how `P0`/`P1` are estimated.

## 6. Outlier suppression + control-derived background

- **Patient-specific outlier suppression (`SuppressOutlierLoci` / INVAR2 `repolish`):** estimate the
  sample ctDNA fraction `P_ESTIMATE`, then per locus test `BINOMIAL_PROB = P(X ≥ x | n = DP, p = P_ESTIMATE)`
  (one-sided) against a **Bonferroni** threshold `outlierSuppression(0.05) / n_loci`; a locus is an
  **outlier (removed)** when its mutant-read count is a one-sided binomial outlier above the sample
  estimate. A locus with no mutant reads (`x ≤ 0`) has tail probability 1 and is never an outlier.
  Oracle: 9 clean loci (1 alt, tail 0.632, kept) + 1 planted (50 alt, tail 3.7e-66, removed) → after
  removal the residual is pure background, IMAFv2 = 0.
- **Locus-noise filtering + background estimation (`EstimateLocusBackground` / `PassesBothStrandsFilter`):**
  the per-locus background error `e` is **estimated from control plasma**: over control samples,
  `BACKGROUND_AF = Σ(ALT_F+ALT_R)/ΣDP`; `LOCUS_NOISE.PASS = (fraction of controls with signal) < 0.1
  AND BACKGROUND_AF < 0.01`. The **both-strands filter** passes iff `ALT_F > 0 AND ALT_R > 0` (or
  `AF == 0` vacuously) — a single-strand-only signal is strand-biased and fails. Oracle: the pooled
  control allele fraction recovers the injected per-locus error (0.002); recurrent (5/20) or
  high-background (0.0125) loci fail.

## 7. Corner cases and failure modes

- **Exactly 1 variant detected** ⇒ MRD-negative (below the ≥2 threshold).
- **< 8 markers tracked** ⇒ sensitivity at ≤0.1% VAF compromised (sensitivity, not the rule).
- **Empty panel / no informative locus (all tumour AF = 0) / empty control panel** ⇒ undefined input;
  out-of-range tumour AF / background throw.
- **Pure-background sample** ⇒ EM `p̂ ≈ 0`, `LR ≈ 0` ⇒ not detected (background subtraction removes noise).
- **Flat size profile** (P1 == P0) ⇒ with-RL GLRT reduces exactly to the no-size GLRT.
- **Both-strands** — `AF == 0` (no alt reads) passes vacuously; alt on a single strand only fails.

## 8. Assumptions and scope

- **ASSUMPTION — per-variant "detected" criterion.** A tracked variant counts as *detected* at **≥ 1
  supporting alt read (default, configurable)** — the minimal source-consistent presence rule. The
  sources define positivity at the **panel** level (≥2 variants) and require per-locus signal above
  background but publish no universal per-locus read-count cutoff (instrument/error-model specific).
  Correctness-affecting only for the per-variant flag (a tunable threshold); it does not change the
  panel-level ≥2 rule.
- **RESOLVED — discrete-vs-KDE size profile** (§5): discrete `COUNT/TOTAL` is the default; the Gaussian
  KDE is opt-in; all mixture/EM/GLRT equations unchanged.

Sources are mutually consistent — the Signatera ≥2/16 positivity rule (Reinert 2019 / PMC9265001), the
Poisson panel LoD (Natera / Avanzini, shared with [[ctdna-detection-and-tumor-fraction]]), the INVAR
integrate-across-loci principle (Wan 2020), and the exact INVAR2 GLRT/EM/IMAFv2/size-weighting/outlier/
background formulas each cover a disjoint stage. **Not for clinical or diagnostic use.**
