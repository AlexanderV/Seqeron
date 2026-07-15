---
type: concept
title: "Somatic variant calling (tumor vs matched-normal VAF classification)"
tags: [oncology, algorithm]
sources:
  - docs/algorithms/Oncology/Somatic_Mutation_Calling.md
  - docs/Evidence/ONCO-SOMATIC-001-Evidence.md
  - docs/Evidence/VARIANT-CALL-001-Evidence.md
source_commit: 77ba02823160a6580a9db64a5a86a1c49f733432
created: 2026-07-10
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-somatic-001-evidence
      evidence: "Test Unit ID: ONCO-SOMATIC-001 ... Algorithm: Somatic Mutation Calling (tumor vs matched normal classification)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:sequencing-artifact-detection
      source: onco-somatic-001-evidence
      evidence: "Yan 2021: putative mutations at <=5% VAF are frequently sequencing errors; the FilterArtifacts technical-artifact QC filter cleans the raw somatic call set this caller produces before clinical interpretation."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:clonal-hematopoiesis-cfdna-filtering
      source: onco-somatic-001-evidence
      evidence: "Dataset D (CHIP-like in normal, f_n=0.03) is classified Germline by the normal-VAF rule; the CHIP filter is the biological-origin QC filter downstream of this caller for cfDNA."
      confidence: medium
      status: current
---

# Somatic variant calling (tumor vs matched-normal VAF classification)

The **foundational somatic SNV caller** of the Oncology family — the **thirty-third ingested
ONCO-\* unit** — that classifies each candidate variant as **Somatic / Germline / NotDetected** by
comparing the **tumor allele frequency `f_t`** against the **matched-normal allele frequency `f_n`**,
plus `FilterGermlineVariants` (keep only the somatic calls). Every variant frequency is
`VAF = altReads / totalReads` — the empirical VAF primitive owned by
[[variant-allele-frequency-and-binomial-ci]] (ONCO-VAF-001), which also gives its Wilson binomial
confidence interval. This is the **upstream caller** that produces the raw somatic call set
which the two pre-interpretation QC filters clean — the technical-artifact
[[sequencing-artifact-detection]] (`FilterArtifacts`) and the biological-origin
[[clonal-hematopoiesis-cfdna-filtering]] (`FilterCHIP`) — before any variant reaches the
clinical-significance layers [[clinical-actionability-oncokb-levels]] and
[[cancer-variant-tier-classification-amp-asco-cap]]. Validated under test unit **ONCO-SOMATIC-001**;
the literature-traced record is [[onco-somatic-001-evidence]], [[test-unit-registry]] tracks the unit,
and [[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade
([[scientific-rigor|research-grade]]), **not for clinical or diagnostic use**.

## 1. The somatic state (Saunders 2012 / Strelka)

Strelka defines the somatic state formally as

```
S = { (f_t, f_n) : f_t ≠ f_n }
```

— the **tumor and normal allele frequencies must differ** — **restricted to a homozygous-reference
normal genotype** (quality reflects the joint `P(S, G_n = ref/ref | D)`). The ref/ref restriction is
load-bearing: the *raw* somatic probability alone over-calls in **LOH / copy-number regions**, so a
variant only counts as somatic when it is essentially **present in the tumor and absent in the
normal**. The successor **Strelka2** (Kim 2018) represents **continuous allele frequencies** for both
samples and classifies somatic by requiring a somatic LOD **and** a somatic VAF to pass — the same
tumor-vs-normal comparison this unit implements as an explicit rule.

## 2. Three decision thresholds

The classifier reduces the Bayesian callers' machinery to three transparent, **configurable** rules:

| Rule | Meaning | Value | Source |
|------|---------|-------|--------|
| **Tumor LoD** | below this f_t → **NotDetected** | **f_t ≥ 0.05** (5%) | Yan 2021 (WES LoD 5%; ≤5% frequently sequencing errors) |
| **Normal "absent" ceiling** | f_n at/below this = absent in normal | **f_n ≤ 0.01** (1%) — `normalVafThreshold` | ASSUMPTION: sub-percent normal band (Yan 2021), ref/ref (Saunders 2012) |
| **Germline separation** | f_n above the ceiling with substantial f_t → **Germline** | f_n > 0.01 | Mutect2 germline filter (Benjamin 2019) |

Decision logic per variant:

- **f_t < 0.05** → **NotDetected** (below the limit of detection; sub-clonal calls at ≤5% VAF are a
  documented high-false-positive band).
- **f_t ≥ 0.05 AND f_n ≤ 0.01** → **Somatic** (present in tumor, absent in normal — `f_t ≠ f_n` with a
  ref/ref-consistent normal).
- **f_t ≥ 0.05 AND f_n > 0.01** → **Germline** (the variant is also present in the normal — Mutect2
  skips variants clearly present in the germline).

## 3. Somatic score (documented simplification)

The reported per-call score is the **monotone separation**

```
somatic score = max(0, f_t − f_n)   ∈ [0, 1]
```

a transparent surrogate for "present in tumor, absent in normal". **ASSUMPTION:** no authoritative
source publishes this exact closed form — Mutect2/Strelka use Bayesian LOD scores not fully retrievable
— so this is documented **as a simplification, not presented as a caller LOD**. A germline variant
(`f_t ≈ f_n`) scores 0; a clear somatic variant (`f_n = 0`) scores `f_t`.

## 4. Tumor-only mode (Mutect2 ℓ_n = 1)

When there is **no matched normal** (no normal coverage), Mutect2's likelihood collapses: **"If we have
no matched normal, ℓ_n = 1"**. Classification then relies on the tumor evidence alone — a variant with
`f_t ≥ 0.05` and no normal to contradict it is called **Somatic** (dataset E: f_t = 0.20, normal 0/0 →
Somatic, score 0.20). The `--af-of-alleles-not-in-resource` prior differs by mode (1e-6 tumor-normal,
5e-8 tumor-only, 4e-3 mitochondrial), and FilterMutectCalls separates germline from somatic in impure
samples by **deviation from the diploid heterozygous fraction 1/2**.

## 5. Worked oracles

`VAF = altReads/totalReads`, score `= max(0, f_t − f_n)`:

| Variant | f_t | f_n | Status | Score | Why |
|---------|-----|-----|--------|-------|-----|
| A clear somatic | 0.25 | 0.00 | **Somatic** | 0.25 | tumor present, normal absent |
| B germline het | 0.48 | 0.50 | **Germline** | 0.00 | f_n > 1%, `f_t ≈ f_n` |
| C sub-LoD | 0.02 | 0.00 | **NotDetected** | 0.00 | f_t < 0.05 |
| D CHIP-like in normal | 0.30 | 0.03 | **Germline** | 0.00 | f_n = 3% > 1% ceiling |
| E tumor-only | 0.20 | 0.00 (0/0) | **Somatic** | 0.20 | no matched normal, ℓ_n = 1 |
| F at tumor threshold | 0.05 | 0.00 | **Somatic** | 0.05 | f_t = 0.05 boundary present (inclusive ≥) |
| G at normal threshold | 0.30 | 0.01 | **Somatic** | 0.29 | f_n = 0.01 boundary absent (inclusive ≤) |

Boundary conventions: **f_t = 0.05 is detected** (inclusive `≥`); **f_n = 0.01 is absent** (inclusive
`≤`). `FilterGermlineVariants` returns only the Somatic calls — its output is a **subset of the input**
(the germline/not-detected calls are removed), the same subset contract as the downstream QC filters.

The four entry points all live in `OncologyAnalyzer` (`Seqeron.Genomics.Oncology`):
`CallSomaticMutations(variants, τ_t, τ_n)` classifies every variant in one **O(n)** pass returning an
order-preserving list, `Classify` does a single variant, `FilterGermlineVariants` returns the somatic
subset, and `CalculateSomaticScore` computes the `[0, 1]` separation. Classification is **pure
in-memory** over `VariantObservation` records (tumor/normal alt+total read counts — no VCF parsing in
this unit, and no suffix-tree/substring machinery). An uncovered site (`totalReads = 0`) yields
`VAF = 0` (INV-05) — a tumor-only / uncovered-normal signal rather than an error — and thresholds
outside `[0, 1]`, null input, or `alt > total` are contract failures (`ArgumentOutOfRangeException` /
`ArgumentNullException`).

## 6. Relationship to the rest of the Oncology family

This caller sits at the **head of the somatic-SNV pipeline**. It is the tumor-vs-matched-normal VAF
counterpart to the **germline, reference↔query** caller [[germline-variant-calling-snp-indel]]
(VARIANT-CALL-001) — same goal (detect variants), different evidence — and both feed the same
downstream VEP-style annotator [[variant-effect-annotation-vep]]. Its raw call set flows into two
pre-interpretation QC filters that both output a subset:

- [[sequencing-artifact-detection]] removes **technical** false positives (OxoG / FFPE deamination /
  strand bias) — a sub-5% VAF somatic-looking call here is exactly what its GIV / FisherStrand rules
  guard against.
- [[clonal-hematopoiesis-cfdna-filtering]] removes **biological** (blood-clone / CHIP) false positives
  in cfDNA. Note dataset **D** (CHIP-like, f_n = 3%): this caller already classifies it **Germline** via
  the normal-VAF rule when a matched normal is present; the CHIP filter is the counterpart when the
  confounding clone is in matched WBC rather than the germline. On the same cfDNA input,
  [[ctdna-detection-and-tumor-fraction]] is the quantification counterpart.

Cleaned somatic variants then feed the clinical-significance units
[[clinical-actionability-oncokb-levels]] and [[cancer-variant-tier-classification-amp-asco-cap]], the
gene-level [[driver-gene-classification-20-20-rule]], the mutation-spectrum
[[sbs96-mutational-signature-catalog]], and every clonal-structure layer
([[allele-specific-copy-number-ascat]], [[cancer-cell-fraction-clonal-clustering]]) — a mis-classified
germline or artifact left in the call set would corrupt all of them.

## 7. Corner cases and scope

- **LOH / copy-number regions** — raw somatic probability over-calls; the ref/ref normal restriction
  (f_n ≤ 1% ceiling) mitigates.
- **Low tumor purity** — allele fractions drift from 1/2; germline separation degrades (Mutect2 uses
  deviation from 1/2). Out of scope for the rule-based classifier, which uses the normal VAF directly.
- **Sub-5% VAF subclonal** — NotDetected at the standard LoD (documented false-positive band).
- **Tumor-only mode** — no normal coverage → ℓ_n = 1, tumor-only classification.
- **Custom thresholds** — `normalVafThreshold` (and the tumor LoD) are caller knobs; changing them
  changes classification (a COULD-test).
- **Null / empty inputs** — documented API-contract failure modes.

**Two flagged assumptions**, both source-consistent and exposed as parameters rather than hard-wired:
the **1% normal "absent" ceiling** (Strelka's ref/ref restriction is Bayesian with no published
scalar; 1% sits in Yan 2021's sub-percent band) and the **`max(0, f_t − f_n)` somatic score** (a
bounded monotone surrogate, not a caller LOD). Sources are mutually consistent — Strelka / Strelka2
(the `f_t ≠ f_n` ref/ref definition and continuous-VAF comparison), Mutect2 (the germline filter and
tumor-only ℓ_n = 1), and Yan 2021 (the 5% tumor LoD and sub-percent normal band) each cover a disjoint
facet. **Not for clinical or diagnostic use.**
