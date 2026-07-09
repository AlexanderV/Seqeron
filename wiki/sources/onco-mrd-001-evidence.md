---
type: source
title: "Evidence: ONCO-MRD-001 (tumor-informed MRD detection — Signatera positivity + INVAR GLRT)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-MRD-001-Evidence.md
sources:
  - docs/Evidence/ONCO-MRD-001-Evidence.md
source_commit: 13a9268a9bcc252ce0f9ce183815411f74f1ff5b
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-MRD-001

The validation-evidence artifact for test unit **ONCO-MRD-001** — **Minimal (Molecular) Residual
Disease detection**, the tumor-informed panel-level ctDNA MRD call. The **twenty-third ingested unit
of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in
[[tumor-informed-mrd-detection]]; [[test-unit-registry]] tracks the unit. It is the **calling /
statistical-detection** liquid-biopsy layer that sits on top of the Poisson quantification primitive
[[ctdna-detection-and-tumor-fraction]] (ONCO-CTDNA-001) — the same cfDNA input, but here a
multi-variant *positive/negative verdict* rather than a single-reporter probability.

## What this file records

- **Online sources (six; mutually consistent — no contradictions):**
  - **Reinert et al. (2019)** *JAMA Oncology* 5(8):1124 (rank 1, primary clinical study behind the
    Signatera assay): tumor-informed personalized multiplex-PCR MRD; post-surgery ctDNA-positive
    patients **~7× more likely to relapse (HR 7.2; 95% CI 2.7–19.0; p < 0.001)** — the post-treatment
    MRD-positive call is the clinically meaningful output.
  - **Natera Signatera analytical-validation white paper (2020)** (rank 3, reference implementation /
    assay spec): a bespoke assay of **16 tumor-specific clonal somatic SNVs** per patient; selected by
    clonality/detectability/frequency from WES of primary tumor + matched normal (buffy coat); tracked
    longitudinally by patient-specific 16-plex PCR + ultra-deep NGS; **Poisson detection limit
    `p = 1 − e^(−nfm)`** (n = haploid genome equivalents, f = ctDNA VAF, m = tumor somatic mutations —
    the same model as ONCO-CTDNA-001 with m = k reporters); sensitivity scales with panel size
    (compromised at ≤0.1% VAF when tracking ≤8 mutations; **LoD 0.01% VAF; analytic specificity
    >99.5%**).
  - **Tumor-informed ctDNA review (PMC9265001)** (rank 4, quoting Reinert 2019 Table 1): the
    **positivity rule** — plasma is ctDNA-positive (MRD-positive) when **≥ 2 of the 16 tracked variants
    are detected**; fewer than 2 ⇒ MRD-negative.
  - **Wan et al. (2020) INVAR** *Sci. Transl. Med.* 12(548):eaaz8084 (rank 1): INtegration of VAriant
    Reads — analyze hundreds–thousands of tumor-informed mutations and **integrate signal across all
    loci** (down to ~1 mutant molecule per 100,000); **IMAF** = depth-weighted average of per-locus
    mutant fractions, background-corrected.
  - **INVAR2 reference implementation (nrlab-CRUK/INVAR2)** (rank 3, Rosenfeld lab, verbatim R
    formulas): the per-locus mixture model, EM ctDNA-fraction estimator, generalised-likelihood-ratio
    detection statistic, IMAFv2, fragment-size weighting, KDE size profile, patient-specific outlier
    suppression, locus-noise filtering + control-derived background estimation (see the concept page
    for the exact equations).
  - **Silverman (1986)** *Density Estimation* (rank 1) + R `density`/`bw.nrd0` (rank 3): the Gaussian
    kernel estimator (eq. 2.2a), `∫K = 1` normalisation, and the rule-of-thumb bandwidth
    `0.9·min(σ̂, IQR/1.34)·n^(−1/5)` behind the opt-in KDE-smoothed size profile.
  - (Avanzini 2020 cited for the Poisson `p = 1 − e^(−λ)` shared with ONCO-CTDNA-001; Lanczos 1964 for
    the log-gamma binomial coefficient.)

- **Methods validated:**
  - **`DetectMRD`** — panel-level positivity: **MRD-positive iff ≥ 2 tracked variants detected**
    (threshold parameterized); reports `DetectedVariantCount` / `TrackedVariantCount`.
  - **`TrackVariantsOverTime`** — longitudinal per-timepoint MRD status + first-positive timepoint flag.
  - Panel **Poisson detection probability** `p = 1 − e^(−n·f·m)` (delegates to the ONCO-CTDNA-001
    primitive for panel size m).
  - **`IntegratedMutantAlleleFractionV2`** — per-context background-subtracted, depth-weighted mean
    `weighted.mean(max(0, MEAN_AF − BACKGROUND_AF), TOTAL_DP)`; a locus with VAF ≤ background contributes 0.
  - **`EstimateInvarSignal`** — INVAR generalised-likelihood-ratio detection: per-locus mixture
    `q = p·g + e·(1−p)`, EM estimate of ctDNA fraction `p`, statistic `LR = logL(p̂) − logL(0)`;
    signal-to-noise (AF-)weighted.
  - **`EstimateInvarSignalWithSize`** / **`FragmentSizeProfile`** (+ `FromKernelDensity`) /
    **`InvarMolecule`** — fragment-size-weighted GLRT (short tumour fragments up-weighted); discrete
    empirical profile default, opt-in Gaussian-KDE profile.
  - **`SuppressOutlierLoci`** (repolish) — one-sided binomial outlier test with a Bonferroni threshold.
  - **`EstimateLocusBackground`** / **`PassesBothStrandsFilter`** — control-derived background-error
    estimation + locus-noise + both-strands filters.

- **Documented corner cases / failure modes:** exactly 1 variant detected ⇒ MRD-negative (below the
  ≥2 threshold); < 8 markers ⇒ sensitivity at ≤0.1% VAF compromised (affects sensitivity, not the rule);
  empty panel / no informative locus (all tumour AF = 0) / empty control panel ⇒ undefined input;
  pure-background sample ⇒ EM `p̂ ≈ 0`, `LR ≈ 0` ⇒ not detected; flat size profile (P1 == P0) ⇒ the
  with-RL GLRT reduces exactly to the no-size GLRT; outlier `x ≤ 0` (no mutant reads) ⇒ binomial tail 1
  ⇒ never an outlier; both-strands `AF == 0` passes vacuously, single-strand-only fails.

- **Datasets (deterministic, computed independently of the C# implementation):**
  - **Signatera positivity:** 2/16 → positive, 1/16 → negative, 0/16 → negative, 3/16 → positive.
  - **Poisson panel:** n=1000, f=0.001, m=1 (λ=1) → 0.6321205588; m=16 (λ=16) → 0.9999998875.
  - **INVAR GLRT synthetic recovery:** inj=0 → p̂≈0, LR≈0 (not detected); inj=0.01 → p̂≈0.01002, LR≈4.06;
    inj=0.02 → LR≈11.81; inj=0.05 → LR≈44.14 (monotone in signal).
  - **AF-weighting:** AF-weighted LR ≈ 2.66 > flat-mean-AF LR ≈ 1.91 (same low-signal mixture).
  - **Size-weighting:** with-RL LR ≈ 0.1969 > no-size LR ≈ 0.1478 (short tumour fragments).
  - **Outlier suppression (repolish):** 9 clean loci (1 alt each, binomial tail 0.632, kept) + 1 planted
    (50 alt, tail 3.7e-66, removed); after removal residual = pure background ⇒ IMAFv2 = 0.
  - **Control-derived background / locus-noise:** pooled control allele fraction recovers the injected
    per-locus error (0.002); recurrent (5/20 controls) or high-background (0.0125) loci fail LOCUS_NOISE.
  - **KDE size profile:** single obs at x₀=100, h=0.5 → normalised P(100)=0.68454, P(101)=P(99)=0.15773
    (analytic Gaussian-bin integral via Φ, sums to 1).

- **Coverage recommendations (14 items):** MUST — `DetectMRD` positive iff ≥2 detected (0/1/2/3);
  correct Detected/Tracked counts; panel Poisson `p = 1 − e^(−nfm)`; `IMAF` depth-weighted mean VAF;
  `IMAFv2` background-subtracted depth-weighted (VAF ≤ background → 0); `EstimateInvarSignal` pure
  background → p̂≈0/LR≈0/not detected, and recovers injected fraction; AF-weighted LR ≥ unweighted; LR
  monotone in signal. SHOULD — custom positivity threshold shifts the call; `TrackVariantsOverTime`
  per-timepoint status + first-positive flag; detectionThreshold gates the call; out-of-range AF /
  background / empty panel throw. COULD — null/empty panel + invalid threshold raise documented
  exceptions.

## Deviations and assumptions

- **ASSUMPTION — per-variant "detected" criterion.** A tracked variant is counted *detected* in plasma
  at **≥ 1 supporting alt read (default, configurable)** — the minimal source-consistent presence rule.
  The cited sources define positivity at the **panel** level (≥2 variants) and require per-locus signal
  above background but publish **no** universal per-locus read-count cutoff (it is instrument/error-model
  specific, e.g. INVAR's trinucleotide GLRT). Correctness-affecting only for the per-variant flag, which
  is a tunable threshold; it does not change the panel-level ≥2 rule.
- **RESOLVED — discrete-vs-KDE size profile.** `FragmentSizeProfile(...)` keeps the **discrete** empirical
  `COUNT/TOTAL` per-length proportion (with uniform fall-back) as the **default**; `FromKernelDensity(...)`
  adds the **opt-in** Gaussian-KDE weight (Silverman 1986 eq. 2.2a, `bw.nrd0` bandwidth, analytic Φ bin
  integral) matching INVAR2's `estimate_real_length_probability`. All mixture/EM/GLRT equations unchanged.

No source contradictions: the Signatera ≥2/16 positivity rule (Reinert/PMC9265001), the Poisson panel
LoD (Natera/Avanzini, shared with ONCO-CTDNA-001), the INVAR integrate-across-loci principle (Wan 2020),
and the exact INVAR2 GLRT/EM/IMAFv2/size-weighting/outlier/background formulas each cover a disjoint stage
and reinforce one another. **Not for clinical or diagnostic use.**
