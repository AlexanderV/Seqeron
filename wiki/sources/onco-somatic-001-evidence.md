---
type: source
title: "Evidence: ONCO-SOMATIC-001 (somatic mutation calling — tumor vs matched-normal VAF classification)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-SOMATIC-001-Evidence.md
sources:
  - docs/Evidence/ONCO-SOMATIC-001-Evidence.md
source_commit: cd2346b7cf5bc7b9f84f0e5cfa7716f539b2cbce
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-SOMATIC-001

The validation-evidence artifact for test unit **ONCO-SOMATIC-001** — **Somatic Mutation Calling**
(tumor vs matched-normal VAF classification): classify each candidate variant as **Somatic /
Germline / NotDetected** by comparing tumor allele frequency `f_t` against normal allele frequency
`f_n`, plus `FilterGermlineVariants` (return only somatic calls). This is the **thirty-third ingested
unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct rule-based caller is
synthesized in its own concept, [[somatic-variant-calling-tumor-normal]]; [[test-unit-registry]]
tracks the unit. This is the **upstream SNV caller** whose raw call set the two pre-interpretation QC
filters clean — the technical-artifact [[sequencing-artifact-detection]] and the biological-origin
[[clonal-hematopoiesis-cfdna-filtering]].

## What this file records

- **Online sources (four; peer-reviewed callers + a limit-of-detection study, mutually consistent):**
  - **Saunders et al. (2012)** "Strelka" — *Bioinformatics* 28(14):1811–1817 (rank 1, PMID 22581179).
    Defines the **somatic state `S = {(f_t, f_n): f_t ≠ f_n}`** (tumor and normal allele frequencies
    must differ), **restricted to a homozygous-reference normal genotype** (quality reflects the joint
    `P(S, G_n = ref/ref | D)` — the raw somatic probability alone over-calls in LOH / copy-number
    regions). Generative model: normal = diploid germline mixture + noise; tumor = normal mixture +
    somatic variation. Thresholds are **empirically chosen**, not absolute probability calibrations.
  - **Kim et al. (2018)** "Strelka2" — *Nature Methods* 15:591–594 (rank 1). Represents **continuous
    allele frequencies** for tumor and normal, leveraging the normal's expected genotype; classifies
    somatic by comparing tumor-vs-normal VAF at each site and requiring somatic LOD + somatic VAF to
    pass.
  - **Benjamin et al. / GATK Mutect2 (2019)** model documentation `mutect.tex` (rank 3, canonical
    reference impl, raw source read). **Germline filter** compares unnormalized probabilities of
    (germline het in normal) / (germline hom-alt) / (somatic, tumor-only). **"If we have no matched
    normal, ℓ_n = 1"** — a somatic variant shows strong tumor evidence (ℓ_t) while the normal shows
    essentially none (ℓ_n ≈ 1); given a matched normal, Mutect2 **skips variants clearly present in the
    germline**. `--af-of-alleles-not-in-resource` default **1e-6** tumor-normal / 5e-8 tumor-only /
    4e-3 mitochondrial. FilterMutectCalls uses the **diploid-het reference fraction 1/2** to separate
    germline from somatic in impure samples.
  - **Yan et al. (2021)** — *Scientific Reports* 11:11640 (rank 1). **Tumor limit of detection (LoD)
    = 5% VAF** (verbatim: WES has an LoD at VAF 5%; putative mutations at ≤ 5% VAF are frequently
    sequencing errors → high false-positive risk). Normal/baseline thresholds sit in the **sub-percent
    band (≈0.01%–0.5%)**; a **1% ceiling for "absent in normal"** is a conservative noise band
    consistent with the ref/ref restriction.

- **Documented corner cases / failure modes:**
  - **Tumor-only mode** (no matched normal): ℓ_n = 1; classification relies on the tumor alone.
  - **Low tumor purity / impure samples:** allele fractions drift from the diploid 1/2 →
    FilterMutectCalls separates somatic from germline by deviation from 1/2.
  - **LOH / copy-number regions:** raw somatic probability over-calls → ref/ref-normal restriction.
  - **Sub-5% VAF subclonal calls:** frequently sequencing errors → treated as NotDetected at the LoD.

- **Datasets (documented oracles — synthetic tumor/normal panel, `VAF = altReads/totalReads`,
  somatic score `= max(0, f_t − f_n)`):**

  | Variant | Tumor VAF f_t | Normal VAF f_n | Expected | Score |
  |---------|---------------|----------------|----------|-------|
  | A clear somatic | 0.25 (25/100) | 0.00 (0/100) | Somatic | 0.25 |
  | B germline het | 0.48 (48/100) | 0.50 (50/100) | Germline | 0.00 |
  | C sub-LoD | 0.02 (2/100) | 0.00 | NotDetected | 0.00 |
  | D CHIP-like in normal | 0.30 (30/100) | 0.03 (3/100) | Germline | 0.00 |
  | E tumor-only (no normal cov) | 0.20 (20/100) | 0.00 (0/0) | Somatic | 0.20 |
  | F at tumor threshold | 0.05 (5/100) | 0.00 | Somatic | 0.05 |
  | G at normal threshold | 0.30 (30/100) | 0.01 (1/100) | Somatic | 0.29 |

## Deviations and assumptions

- **ASSUMPTION 1 — normal "absent" ceiling = 1% VAF.** Strelka restricts to a ref/ref normal genotype
  but publishes no single VAF cutoff (its decision is Bayesian). A concrete rule-based cutoff is
  required; **1%** sits in the sub-percent normal baseline band (Yan 2021) and is exposed as a
  configurable `normalVafThreshold`, so the correctness-critical decision is **parameterized, not
  hard-wired to an invented value**.
- **ASSUMPTION 2 — somatic score `= max(0, f_t − f_n)`.** No authoritative source publishes this exact
  closed-form scalar (Mutect2/Strelka use Bayesian LOD scores not retrievable in full). The **monotone
  separation score** is a transparent, bounded-[0,1] surrogate for "present in tumor, absent in
  normal", documented as a simplification, **not presented as a caller LOD**.
- **Coverage recommendations:** MUST-test clear somatic (high f_t / zero f_n → Somatic, score = f_t −
  f_n); germline (comparable f_t/f_n → Germline); sub-5% tumor VAF → NotDetected; tumor-only (ℓ_n = 1)
  → Somatic; boundaries (f_t = 0.05 present, f_n = 0.01 absent); `FilterGermlineVariants` returns only
  somatic calls. SHOULD-test CHIP-like low-VAF normal contamination (f_n just above 1% → Germline).
  COULD-test custom thresholds change classification.

No source contradictions — Strelka (ref/ref normal, `f_t ≠ f_n`), Strelka2 (continuous VAF, somatic
LOD + VAF), Mutect2 (germline filter, ℓ_n = 1 tumor-only, 1/2 diploid-het fraction), and Yan 2021 (5%
tumor LoD, sub-percent normal band) each cover a disjoint facet of the same tumor-normal comparison and
reinforce one another.
