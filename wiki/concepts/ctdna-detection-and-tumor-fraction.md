---
type: concept
title: "ctDNA detection (Poisson LoD) and tumor-fraction estimation"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-CTDNA-001-Evidence.md
source_commit: d40f826d6627e4defc37d0248ca0911eec5bffdf
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-ctdna-001-evidence
      evidence: "Test Unit ID: ONCO-CTDNA-001; Algorithm: ctDNA Analysis (Poisson limit-of-detection, tumor-fraction estimation, mean variant-allele-fraction summarization)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:clonal-hematopoiesis-cfdna-filtering
      source: onco-ctdna-001-evidence
      evidence: "Both operate on the same cell-free DNA / liquid-biopsy input (Newman 2014 CAPP-Seq plasma ctDNA); ONCO-CTDNA-001 quantifies detection probability + tumor fraction while ONCO-CHIP-001 removes clonal-hematopoiesis false positives — complementary stages of the cfDNA pipeline."
      confidence: high
      status: current
---

# ctDNA detection (Poisson LoD) and tumor-fraction estimation

The **quantification / limit-of-detection layer** of the Oncology family for **circulating tumor DNA
(ctDNA)** in a **liquid biopsy** (cell-free DNA / cfDNA from plasma). Where the sibling
[[clonal-hematopoiesis-cfdna-filtering]] *filters out* non-tumor cfDNA calls, this unit *quantifies*
the tumor signal: the **probability of detecting** ≥1 tumor molecule under Poisson sampling, the
**tumor fraction** from clonal VAF, and the **mean VAF** across reporters — the quantification
primitive beneath **minimal residual disease (MRD)** monitoring. The multi-variant MRD *calling*
layer that reuses this Poisson primitive is [[tumor-informed-mrd-detection]] (ONCO-MRD-001).
Validated under test unit **ONCO-CTDNA-001**; the
literature-traced record is [[onco-ctdna-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade
([[scientific-rigor|research-grade]]), **not for clinical or diagnostic use**.

## 1. The Poisson limit-of-detection model (Patent US 11,085,084 / Avanzini 2020)

ctDNA detection is fundamentally a **rare-molecule sampling** problem. With **n** sequenced **genome
equivalents** and a mutant allele fraction **d** (the fraction of cfDNA molecules carrying the
variant — the assay's detection limit), the number of observed tumor molecules is **Poisson with
mean λ = n·d**. Detection probabilities:

- **Single reporter:** `x = 1 − e^(−n·d)` — probability of observing ≥1 ctDNA molecule.
- **k independent reporters:** `p = 1 − e^(−n·d·k)` — integrating multiple tracked variants raises
  sensitivity (Newman 2014's multi-reporter design).

`DetectionProbability(n, d, k)` returns this exact **p**. The **low-burden regime λ < 3** is
Poisson-sampling-governed: small changes in plasma input or recovery shift results across the LoD, so
a true positive can be **missed purely from sampling** — detection is a *probability, not a guarantee*.

**Detectability decision.** The model yields a probability, not a verdict; the literature fixes **no**
universal threshold for *declaring* a variant detected. The unit therefore exposes a deterministic
detectability test: `p ≥ τ` (caller-supplied, **default 0.95** = the 95%-sensitivity assay
convention) **AND** `λ = n·d·k ≥ 1` (at least one mutant molecule expected). The returned probability
is exact; only the boolean flag depends on τ.

## 2. Tumor fraction = 2 × VAF (copy-neutral diploid; CNAqc / Antonello 2024)

`CalculateTumorFraction` = **2 × mean clonal heterozygous VAF**. For a clonal heterozygous SNV on a
copy-neutral diploid background the expected-VAF relation `v = m·π / [2(1−π) + π·n_tot]` (m=1, n_tot=2)
reduces to **v = π/2**, so the observed VAF is *half* the cellular prevalence and **tumor fraction =
2·VAF**. This is the same diploid-heterozygous identity that
[[allele-specific-copy-number-ascat]]/[[cancer-cell-fraction-clonal-clustering]] use for purity from
VAF (`EstimatePurityFromVaf`), specialized here to the copy-neutral liquid-biopsy case.

## 3. Mean VAF and genome-equivalents conversion

- **`CalculateMeanVaf`** = mean of `altReads / totalReads` across reporters — Newman 2014 summarizes a
  patient's ctDNA level as a fraction **across the SNV/indel reporter set**.
- **Genome-equivalents helper** — converts cfDNA mass to the molecule count **n** the Poisson model
  needs: **one haploid genome = 3.3 pg** (Devonshire 2014) ⇒ **≈ 303 haploid GE per ng cfDNA**
  (1000/3.3, corroborated by Alcaide 2020).

## 4. Worked oracles

| Quantity | Inputs | Result |
|----------|--------|--------|
| DetectionProbability | n=15000, d=0.001, k=1 (λ=15) | 1 − e⁻¹⁵ ≈ 0.99999969 |
| DetectionProbability | λ=0 (n=0 or d=0) | 1 − e⁰ = 0 → not detected |
| Below-LoD detect | d < 0.025%, λ < 1 | detect = false |
| k>1 sensitivity | fixed n,d, k=10 vs k=1 | p(k=10) > p(k=1) |
| TumorFraction | clonal het VAF 0.10 | 0.20 |
| GE conversion | 1 ng cfDNA | ≈ 303 haploid GE |
| GE conversion | 3.3 pg | 1 GE |

**Detection range (Newman 2014):** validated 0.025%–10% mutant allele fraction; specificity floor
~0.02% @ 96%; median pre-treatment ctDNA fraction ~0.1%; analytic **background error floor** mean
0.006% / median 0.0003% — fractions at/under this are indistinguishable from sequencing error. A
fraction below the assay LoD is reported **not-detected, not zero ctDNA**.

## 5. Relationship to the rest of the Oncology family

Operates on the **same cfDNA input** as [[clonal-hematopoiesis-cfdna-filtering]] but answers the
complementary question — *how much tumor is there / is it detectable?* rather than *is this call
tumor-derived?* Tumor fraction and CCF share the diploid-heterozygous VAF inversion with the
clonal-structure layers [[allele-specific-copy-number-ascat]] and
[[cancer-cell-fraction-clonal-clustering]] (this unit is the copy-neutral, whole-sample specialization).
Downstream detected variants feed the clinical-significance units
[[clinical-actionability-oncokb-levels]] and [[cancer-variant-tier-classification-amp-asco-cap]].

## 6. Assumption and scope

- **ASSUMPTION — detection-decision rule.** Only the Poisson probability `p = 1 − e^(−ndk)` is
  non-assumption. The declare-detected threshold (default 0.95) plus the `λ ≥ 1` guard is a flagged,
  source-consistent convention (the same 0.95 used by `CalculateVAFConfidenceInterval`); changing it
  moves only the boolean flag, never the returned probability.
- **Scope** — copy-neutral diploid tumor-fraction (`TF = 2·VAF`); input validation covers negative n,
  d ∉ [0,1], k < 1, null variant set, and VAF > 0.5 for tumor fraction.

Sources are mutually consistent — the Poisson LoD model (Patent/Avanzini), the mass→molecule
conversion (Devonshire/Alcaide, both 3.3 pg ⇒ 303/ng), the detection range and across-reporter
fraction (Newman), the λ = n·d worked count (Pessoa), and the TF = 2·VAF diploid relation
(CNAqc/Antonello) each cover a disjoint stage. **Not for clinical or diagnostic use.**
