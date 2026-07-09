---
type: concept
title: "Allele-specific copy number + tumor purity/ploidy fit (ASCAT / ASPCF / subclonal / multiplicity)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-ASCAT-001-Evidence.md
source_commit: ee6df31c1d54006d2fcc189bb794df7d37c89cf0
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-ascat-001-evidence
      evidence: "Test Unit ID: ONCO-ASCAT-001 ... Algorithm: Upstream derivation of allele-specific copy-number segments, joint purity/ploidy fit (ASCAT), and mutation multiplicity"
      confidence: high
      status: current
---

# Allele-specific copy number + tumor purity/ploidy fit (ASCAT / ASPCF / subclonal / multiplicity)

The **upstream tumor copy-number layer** of the Oncology family: from per-locus **logR** (total
intensity) and **BAF** (B-allele frequency) tracks, derive **allele-specific integer copy number**
(nA/nB per segment) by jointly fitting tumor **purity ρ** and **ploidy ψ** — the **ASCAT** method
(Van Loo et al. 2010) — plus its segmentation front-end (**ASPCF**), a **subclonal** two-state
extension (Battenberg), and the downstream **mutation multiplicity / cancer-cell-fraction (CCF)**
inversion. Validated under test unit **ONCO-ASCAT-001**; the literature-traced record is
[[onco-ascat-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

This is the tumor **allele-specific / clonal** layer. It sits below the clinical-interpretation ONCO
siblings — the therapeutic ranking [[clinical-actionability-oncokb-levels]], the four-tier
[[cancer-variant-tier-classification-amp-asco-cap]], and the QC filter [[sequencing-artifact-detection]]
— and is distinct from the chromosome-arm-scale total-copy-number [[aneuploidy-detection]] (a log2
depth-ratio caller with no allelic contrast or purity fit). The same logR + BAF allelic-contrast idea,
restricted to the HLA genes, is the **HLA-locus specialization**
[[hla-nomenclature-and-allele-specific-loh]] (LOHHLA — allele-specific HLA loss-of-heterozygosity).
The downstream **genomic-scar** unit [[homologous-recombination-deficiency-score]] (ONCO-HRD-001) reads
LOH / TAI / LST scar counts straight off these allele-specific major/minor CN segments and sums them
into the HRD score; its LOH term is the standalone genome-wide LOH caller
[[loss-of-heterozygosity-detection]] (ONCO-LOH-001), which counts minor-CN-0 / major-CN-≠0 segments
> 15 Mb over the same segments.

## 1. ASCAT joint (ρ, ψ) fit — the core

Each segment carries one fitted **logR `r`** and one/two **BAF `b`** values. Allele-specific copy
numbers follow the ASCAT inversion (reference `ascat.runAscat.R`, γ = platform parameter):

```
nA = (rho-1 - (b-1)*2^(r/gamma) * ((1-rho)*2 + rho*psi)) / rho
nB = (rho-1 +  b   *2^(r/gamma) * ((1-rho)*2 + rho*psi)) / rho
```

with `rho` = aberrant-cell fraction (purity ρ), `psi` = tumor ploidy ψ, `gamma` = platform gain
(**γ = 1 for massively-parallel sequencing — WGS/WES/TS**; ≈0.55 for SNP arrays only).

**Grid-search objective ("sunrise" plot).** ρ and ψ are chosen over a grid to minimise the
segment-length-weighted **squared distance of the minor allele to the nearest non-negative integer**
(germline-heterozygous SNPs *should* land on whole numbers):

```
d = Σ  |nMinor − max(round(nMinor),0)|²  · length · (b==0.5 ? 0.05 : 1)
TheoretMaxdist = Σ  0.25 · length · (b==0.5 ? 0.05 : 1)        # 0.25 = (½)² worst case
goodnessOfFit  = (1 − d/TheoretMaxdist) · 100                   # reported as a %
```

Integer assignment rounds and clamps: `nA = max(round(nAfull),0)`, `nB = max(round(nBfull),0)`; the
**major allele is the larger** of the two. **BAF = 0.5 (balanced) segments are down-weighted ×0.05** —
they carry little allele-specific information (cannot tell 1+1 from 2+2 except via logR).

**ψ scale:** ploidy = DNA amount relative to a haploid genome; pure diploid → ψ = 2, ">2.7n" marks
aneuploidy. The *fitted* ψ here is the joint grid-search output; the **post-hoc** summary of ψ as a
length-weighted mean over already-called segments — plus the whole-genome-doubling flag — is the
distinct downstream unit [[tumor-ploidy-estimation-and-whole-genome-doubling]] (ONCO-PLOIDY-001).

## 2. ASPCF segmentation front-end (Nilsen 2012 / Ross 2021)

Before fitting, the raw tracks are segmented by **penalised least squares (PCF)** — within-segment SSE
plus a per-segment penalty γ that trades fit against parsimony:

```
L(S | y, γ) = Σ_{I∈S} Σ_{j∈I} (y_j − ȳ_I)² + γ·|S|
```

solved to the **global optimum** by an O(n²) dynamic program `e_k = min_j ( d_{jk} + e_{j−1} + γ )`.
**ASPCF** (allele-specific PCF) segments logR and BAF **jointly with common change points** but keeps
**separate per-track means**: `L(S|y₁,y₂,γ) = L(S|y₁,γ) + L(S|y₂,γ)`, γ charged once per segment. BAFs
are **mirrored** (folded to their distance from 0.5) in regions of allelic imbalance before joint
segmentation, so a copy-neutral-LOH segment separates from a balanced segment that share the same logR.
Large γ collapses to one segment; small γ recovers every level. (copynumber default γ = 40; ASCAT
later 70 — probe-scale-specific, so the repository exposes γ as a required caller parameter.)

## 3. Subclonal two-state extension (Battenberg / Nik-Zainal 2012)

A segment is **clonal** (one state, all tumor cells) or **subclonal** (two cell populations, two integer
states). An observed real-valued allele copy number is a fraction-weighted mix of the **two bracketing
integers**:

```
n_obs = f·n₁ + (1−f)·n₂ ,   n₂ = ⌊n_obs⌋ , n₁ = ⌈n_obs⌉ , f = (n_obs − n₂)/(n₁ − n₂) ∈ [0,1]
```

An integer `n_obs` collapses to a single clonal state (f ≈ 0 or 1). Three-or-more-population and
non-adjacent-state mixtures are out of scope (documented limitation).

## 4. Mutation multiplicity and CCF (McGranahan 2016 / PICTograph / DeCiFering)

Once purity and allele-specific CN are known, a somatic mutation's **multiplicity m** (mutated-allele
copies) and **cancer-cell fraction (CCF)** follow from the VAF generative model:

```
VAF = (m · CCF · p) / (c · p + 2·(1−p))                       # PICTograph
m   = VAF·(c·p + 2(1−p)) / p   at clonal CCF = 1              # invert, round, clamp to [1, major-CN]
CCF = (F · v) / (ρ · M) ,  F = ρ·N_tot + 2(1−ρ)               # DeCiFering closed form; v = VAF, M = multiplicity
```

with p/ρ = purity, c/N_tot = tumor total CN, M = multiplicity, 2 = normal CN. **Multiplicity clamp:** a
raw value < 0.5 would round to 0 (no mutated copy) — non-physical for an observed variant — so it is
clamped to ≥ 1 (and ≤ major-allele CN).

This same CCF closed form is the standalone estimator of the **downstream** clonal-structure unit
[[cancer-cell-fraction-clonal-clustering]] (ONCO-CCF-001), which adds the reported-value **[0,1] cap**
(exposing the uncapped raw) and a **deterministic 1D Lloyd k-means** that deconvolutes the per-mutation
CCF vector into clones/subclones — the reconstruction step this copy-number layer stops short of.

## Corner cases and failure modes

- **Non-identifiability / multiple optima:** the sunrise plot can show several local minima (e.g. a
  2n vs 4n solution); ASCAT takes the **global minimum** over the grid. A planted single-solution genome
  must have a unique minimum within tolerance.
- **γ must match the platform** (γ = 1 sequencing, ≈0.55 arrays); a wrong γ rescales logR and biases CN.
- **Balanced (BAF = 0.5) segments** are down-weighted ×0.05 in the goodness-of-fit.
- **ASPCF ≤ greedy:** the DP penalised cost is ≤ any greedy segmentation (it returns the global optimum).

## Planted-truth test strategy

Tests **invert the ASCAT forward model** (γ = 1) to synthesise deterministic inputs from known
(ρ₀, ψ₀, integer nA/nB), then assert recovery:

```
b = (ρ·nB + (1−ρ)·1) / (ρ·(nA+nB) + 2(1−ρ))                    # BAF at a germline-het SNP
r = log2( (ρ·n + 2(1−ρ)) / (ρ·ψ + 2(1−ρ)) )                    # logR, genome-average segment → r = 0
```

Canonical oracles: ρ₀ = 0.80, ψ₀ ∈ {2, 3}; segments 1+1 (balanced, b=0.5, r=0), 2+0 (copy-neutral LOH),
2+1 (gain); a clonal m=1 mutation on a 1+1 (CN=2) segment → CCF ≈ 1.0. ASPCF: a two-level logR track
(0.0 for loci 0–9, 1.0 for 10–19) with γ = 0.5 recovers one breakpoint at index 10. Subclonal:
nA_obs = 1.4, nB_obs = 0.6 → states (2,0)/(1,1), f ≈ 0.4; a pure-clonal 2+1 collapses to one state.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for the allele-specific copy-number
derivation. The forward-model BAF/logR equations are used **only to synthesise planted-truth test
inputs** — the production derivation consumes **measured** logR/BAF. **Not for clinical or diagnostic
use.** No source contradictions: Van Loo 2010 + the ASCAT R reference (core fit), Nilsen 2012 + Ross
2021 (ASPCF segmentation), Nik-Zainal 2012 / Battenberg (subclonal), and McGranahan 2016 / PICTograph /
DeCiFering (multiplicity/CCF) each cover a disjoint stage and are mutually consistent.
