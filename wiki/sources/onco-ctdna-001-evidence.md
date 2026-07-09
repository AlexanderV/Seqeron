---
type: source
title: "Evidence: ONCO-CTDNA-001 (ctDNA detection — Poisson limit-of-detection, tumor-fraction, mean-VAF)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-CTDNA-001-Evidence.md
sources:
  - docs/Evidence/ONCO-CTDNA-001-Evidence.md
source_commit: d40f826d6627e4defc37d0248ca0911eec5bffdf
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-CTDNA-001

The validation-evidence artifact for test unit **ONCO-CTDNA-001** — **ctDNA analysis**
(Poisson limit-of-detection, tumor-fraction estimation, mean variant-allele-fraction
summarization). The **eleventh ingested unit of the Oncology family** and one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct
method is synthesized in [[ctdna-detection-and-tumor-fraction]]; [[test-unit-registry]] tracks the
unit. It is the **quantification / detection-limit** liquid-biopsy layer — the counterpart of the
biological-origin filter [[clonal-hematopoiesis-cfdna-filtering]] on the same cfDNA input.

## What this file records

- **Online sources (seven; mutually consistent, no contradictions):**
  - **Newman et al. (2014)** *Nat. Med.* 20(5):548–554, PMC4016134 (rank 1 — CAPP-Seq primary):
    detection is parameterized by a **mutant allele fraction** — CAPP-Seq "detected defined inputs
    … between 0.025% and 10%" with **96% specificity down to ~0.02%**; **background error floor**
    mean 0.006% / median 0.0003%; observed patient ctDNA fractions **~0.02%–3.2%, median ~0.1%**
    pre-treatment, summarized as a fraction **across SNV/indel reporters**.
  - **US Patent 11,085,084 B2** (rank 2 — formal spec of the detection model; restates the
    **Avanzini et al. 2020** *Sci. Adv.* shedding-model derivation, paywalled this session): the
    **Poisson detection model** — mean **λ = n·d** (n = sequenced **genome equivalents**, d =
    detection limit / mutant allele fraction); single-reporter **x = 1 − e^(−nd)**; **k-reporter
    p = 1 − e^(−ndk)**; **low-burden regime λ < 3** is Poisson-sampling-governed (results shift
    across the LoD from small input/recovery changes).
  - **Avanzini et al. (2020)** *Sci. Adv.* 6(50):eabc4308 (rank 1 — the shedding-model primary the
    patent restates; detection equations corroborated verbatim via the patent).
  - **Alcaide et al. (2020)** *Sci. Rep.* 10:12564, PMC7387491 (rank 1): **1 ng cfDNA ≈ 303 haploid
    genome equivalents** (mass→molecule conversion).
  - **Devonshire et al. (2014)** *Anal Bioanal Chem*, PMC4182654 (rank 1): **one haploid genome =
    3.3 pg** ⇒ 1000/3.3 ≈ 303 copies/ng (consistent with Alcaide).
  - **Pessoa et al. (2023)** review, PMC10314661 (rank 1): worked molecule count — at **VAF 0.1%
    over ~15,000 haploid GE ⇒ 15 tumor molecules** (λ = n·d = 15,000 × 0.001 = 15), corroborating
    λ = n·d.
  - **Antonello et al. (2024)** CNAqc, *Genome Biology* 25:38 (rank 1): the clonal-heterozygous
    diploid identity — expected VAF `v = m·π / [2(1−π) + π·n_tot]` reduces to **v = π/2** for m=1,
    n_tot=2, hence **tumor fraction = 2·VAF** (the same relation used by
    `OncologyAnalyzer.EstimatePurityFromVaf`).
  - (Wan et al. 2017 *Nat Rev Cancer* cited for context only — no value taken.)

- **Methods validated:**
  - **`DetectionProbability`** — the Poisson detection probability **p = 1 − e^(−n·d·k)** for n
    genome equivalents, mutant allele fraction d, k reporters.
  - **Detectability test** — deterministic boolean against a caller-supplied probability threshold
    (default **0.95**, the 95%-sensitivity assay convention) **plus** the requirement **λ = n·d·k ≥ 1**
    (at least one mutant molecule expected). Only the returned probability is non-assumption.
  - **`CalculateTumorFraction`** — **2 × mean clonal-heterozygous VAF** (copy-neutral diploid).
  - **`CalculateMeanVaf`** — mean of `altReads/totalReads` across reporters (Newman's across-reporter
    fraction summary).
  - **Genome-equivalents helper** — 1 ng ⇒ ≈ 303 haploid GE; 3.3 pg ⇒ 1 GE.

- **Documented corner cases / failure modes:** allele fractions at/under the analytic background
  (median 0.0003% / mean 0.006%) are indistinguishable from sequencing error (validated range starts
  0.025%); a fraction below the assay LoD is reported **not-detected, not zero ctDNA**; the
  **Poisson-limited low-input regime λ < 3** makes a true positive a stochastic miss (detection is a
  probability, not a guarantee); **λ = 0** (n=0 or d=0) ⇒ P = 1 − e⁰ = 0.

- **Datasets (deterministic, hand-derived):**
  - **CAPP-Seq range (Newman 2014):** detection 0.025%–10%; specificity floor ~0.02% @ 96%; median
    pre-treatment fraction ~0.1%; background 0.006% / 0.0003%.
  - **Poisson worked example (Pessoa 2023):** n=15,000, d=0.001 ⇒ λ=15, P(≥1) = 1 − e⁻¹⁵ ≈ 0.99999969.
  - **Mass→molecule (Devonshire/Alcaide):** 3.3 pg / haploid genome; ≈ 303 GE per ng cfDNA.

- **Coverage recommendations:** MUST — `DetectionProbability` = 1 − e^(−ndk) (n=15000,d=0.001,k=1);
  `CalculateTumorFraction` = 2 × mean clonal-het VAF; `CalculateMeanVaf` = mean altReads/totalReads;
  GE helper (303/ng, 3.3 pg/GE); λ=0 ⇒ p=0 not-detected; below-LoD (d<0.025%) detect=false when λ<1.
  SHOULD — k>1 raises probability; input validation (negative n, d∉[0,1], k<1, null set, VAF>0.5 for
  tumor fraction). COULD — monotonicity of p in n and d.

## Deviations and assumptions

- **ASSUMPTION — detection-decision rule.** The Poisson model gives a *probability* p = 1 − e^(−ndk);
  the literature fixes **no** universal probability threshold for *declaring* a variant detected. The
  implementation therefore returns the exact **p** (non-assumption) and offers a deterministic
  detectability test against a **caller-supplied threshold (default 0.95**, the same 95% convention
  as `CalculateVAFConfidenceInterval` in this class) plus **λ = n·d·k ≥ 1**. Changing the threshold
  changes only the boolean flag, never the returned probability.

No source contradictions among the seven references: the Poisson LoD model (patent / Avanzini), the
mass→molecule conversion (Devonshire / Alcaide, both giving 3.3 pg ⇒ 303/ng), the detection-range and
across-reporter fraction (Newman), the λ = n·d worked count (Pessoa), and the TF = 2·VAF diploid
relation (CNAqc / Antonello) each cover a disjoint stage and reinforce one another.
