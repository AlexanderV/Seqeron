---
type: concept
title: "Cancer cell fraction (CCF) estimation + 1D clonal/subclonal clustering"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-CCF-001-Evidence.md
source_commit: 280c1e8efb03103002e158f47f67a1e48267d3b0
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-ccf-001-evidence
      evidence: "Test Unit ID: ONCO-CCF-001 ... Algorithm: Cancer Cell Fraction (CCF) point estimation and 1D CCF clustering into clones/subclones"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:allele-specific-copy-number-ascat
      source: onco-ccf-001-evidence
      evidence: "CCF = f·(ρ·N_T + 2(1−ρ)) / (ρ·m) — the estimator consumes sample purity ρ, local total copy number N_T, and mutation multiplicity m, which are the outputs of the upstream ASCAT copy-number / purity fit."
      confidence: high
      status: current
---

# Cancer cell fraction (CCF) estimation + 1D clonal/subclonal clustering

The **downstream clonal-structure layer** of the Oncology family: given per-mutation VAF plus a
tumor's purity and local copy number, estimate each somatic mutation's **cancer cell fraction
(CCF)** — the fraction of tumor cells carrying it — then **cluster** those CCFs into
**clones/subclones** and label the top cluster **clonal**. Validated under test unit
**ONCO-CCF-001**; the literature-traced record is [[onco-ccf-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the evidence-artifact pattern.

This unit sits **on top of** the copy-number/purity substrate
[[allele-specific-copy-number-ascat]] — that page already derives the per-mutation CCF closed form
as its §4 (multiplicity / CCF VAF-inversion). This concept is the **standalone `EstimateCCF`
estimator** (with the reported-value cap) plus the genuinely distinct **`ClusterCCFValues` 1D
deconvolution** that ASCAT does not cover. It is a research-grade correctness reference —
[[scientific-rigor|research-grade]], **not for clinical or diagnostic use**.

## 1. CCF point estimation (McGranahan 2016 / PICTograph / PMC7867630)

Three primary sources give the **same** closed form (mutually corroborating — m and ρ both in the
denominator):

```
CCF = f · (ρ·N_T + 2(1−ρ)) / (ρ · m)          # f = VAF, ρ = purity, N_T = tumor total CN, m = multiplicity
```

derived by inverting the VAF generative model `VAF = (m·CCF·ρ) / (N_T·ρ + 2(1−ρ))` (Zheng 2022,
PICTograph). Equivalently CCF = observed mutation copy number `n_mut` / multiplicity, with
`n_mut = VAF·(1/ρ)·[ρ·CN_t + 2(1−ρ)]` (McGranahan 2016; normal CN = 2). The Box-1 definition
(Tarabichi 2021) states it plainly: **CCF = CP / purity** where CP is cellular prevalence — so the
cluster with the largest CCF is the clonal one (see §3).

**Multiplicity m** is itself `m = f·(ρ·N_T + 2(1−ρ)) / ρ`, rounded to the nearest non-zero integer
for clonal copy-number regions; valid range **1 ≤ m ≤ tumor copy number**.

**Domain guards:** purity ρ ∈ (0,1] (it divides), VAF ∈ [0,1], copy number ≥ 1, multiplicity ∈
[1, tumor CN] — invalid inputs throw.

## 2. The [0,1] reported cap (assumption, source-consistent)

The **raw** formula is not intrinsically bounded ≤ 1: with VAF sampled above its expectation the
CNAqc vignette shows **CCF = 1.06** (VAF 0.471, m=1). The registry invariant is 0 ≤ CCF ≤ 1, so the
unit reports **min(1, raw)** as the bounded CCF while **also exposing the uncapped raw value**.
Justification (an explicit ASSUMPTION, no source forbids it): the McGranahan clonal definition —
a mutation present in **all** cancer cells has CCF = 1 — so capping to 1 is the biologically correct
ceiling. Worked cap case: VAF 0.60, ρ 0.80, N_T 2, m 1 → raw 1.50 → **reported 1.0**.

## 3. 1D CCF clustering → clones/subclones (deterministic Lloyd k-means)

Many mutations' CCFs are **clustered into clones/subclones**. Sources name CCF clustering broadly
(Dirichlet process, variational beta mixtures), but this unit deliberately uses a **deterministic,
fully-specified 1D** method (ASSUMPTION — no fabricated Dirichlet process): **Lloyd's k-means**
(Lloyd 1982) with **quantile seeding**.

- **Assignment step:** each CCF → nearest centroid by least squared Euclidean distance.
- **Update step:** recompute each centroid as the mean of its assigned CCFs.
- **Objective:** minimise within-cluster sum of squares `Σ_i Σ_{x∈S_i} ‖x − μ_i‖²`.
- **Determinism (no RNG):** sort the values and seed the k centroids at evenly-spaced **quantiles** —
  identical output across repeated runs and across shuffled input grouped by value.
- **Clonal-cluster rule (Tarabichi 2021):** "the cluster with the highest CP can be deemed clonal"
  — i.e. the **highest-centroid** cluster is the clonal one; the rest are subclonal lineages.

**Oracle:** CCFs `{1.0, 0.98, 0.96, 0.50, 0.48, 0.52}`, k=2 → centroids `{0.50, 0.98}`, low cluster
= indices {3,4,5}, high = {0,1,2}, **clonal = the 0.98-centroid cluster**. Boundary: k=1 → one
cluster at the overall mean.

## Relationship to the ASCAT copy-number layer

Complementary, not overlapping: [[allele-specific-copy-number-ascat]] **produces** purity ρ, total
copy number N_T, and multiplicity m (its ASCAT/ASPCF/Battenberg stages) and states the CCF formula
as its terminal §4; this unit **consumes** those to estimate + cap the per-mutation CCF and then
**clusters** the resulting CCF vector into clonal/subclonal populations — the reconstruction step
ASCAT stops short of. The two are siblings under the clinical-interpretation ONCO units
[[clinical-actionability-oncokb-levels]], [[cancer-variant-tier-classification-amp-asco-cap]], and
the QC filter [[sequencing-artifact-detection]].

## Corner cases and failure modes

- **CCF > 1 from noise** — raw formula exceeds 1 (CNAqc 1.06); reported value capped, raw exposed.
- **Multi-copy loci (N_T > 2)** — multiplicity m ≠ 1; an integer m (≥1, ≤ tumor CN) must be
  supplied/estimated or CCF is ambiguous.
- **Unknown / invalid purity** — ρ divides, must be in (0,1].
- **Null/empty inputs** — documented failure modes on both `EstimateCCF` and `ClusterCCFValues`.

## Scope and limitations

Sources are mutually consistent: Tarabichi 2021 (*Nat. Methods*), Zheng 2022 (*Bioinformatics*,
PICTograph), McGranahan 2016 (*Science*), CNAqc (reference implementation), and Lloyd 1982 (k-means)
each cover a disjoint stage. **Two flagged assumptions**, both source-consistent rather than
source-mandated: (1) the [0,1] reported cap (invariant + McGranahan clonal definition), and (2)
deterministic Lloyd k-means with quantile seeding as the concrete 1D clustering method (sources name
the clustering only broadly). **Not for clinical or diagnostic use.**
