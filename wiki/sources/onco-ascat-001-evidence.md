---
type: source
title: "Evidence: ONCO-ASCAT-001 (allele-specific copy number + purity/ploidy fit — ASCAT/ASPCF/subclonal/multiplicity)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-ASCAT-001-Evidence.md
sources:
  - docs/Evidence/ONCO-ASCAT-001-Evidence.md
source_commit: ee6df31c1d54006d2fcc189bb794df7d37c89cf0
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-ASCAT-001

The validation-evidence artifact for test unit **ONCO-ASCAT-001** — **upstream derivation of
allele-specific copy-number segments, joint purity/ploidy fit (ASCAT), and mutation multiplicity**. The
**fourth ingested unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its
own concept, [[allele-specific-copy-number-ascat]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (all mutually consistent, spanning four algorithm stages, no contradictions):**
  - **Van Loo et al. (2010)** ASCAT, PNAS 107(39):16910 (rank 1, primary) — two input tracks **logR**
    (total intensity `r`) + **BAF** (allelic contrast `b`); a **grid search over ploidy ψ and aberrant
    fraction ρ** ("sunrise" goodness-of-fit plot); the optimal solution minimises distance of
    allele-specific copy numbers to **non-negative integers** for germline-het SNPs; ψ = DNA relative
    to a haploid genome (pure diploid → 2; ">2.7n" = aneuploidy).
  - **ASCAT R reference `ascat.runAscat.R`** (rank 3, original authors) — the verbatim nA/nB inversion
    `nA=(rho-1-(b-1)*2^(r/gamma)*((1-rho)*2+rho*psi))/rho` and mirror for nB; the length-weighted
    squared minor-allele integer-distance `d` with **BAF=0.5 down-weighted ×0.05**;
    `TheoretMaxdist` (0.25 = (½)² worst case) and `goodnessOfFit = (1 − m/TheoretMaxdist)·100`;
    integer assignment `round`+clamp-0, **major = larger** of {nA,nB}.
  - **ASCAT README — γ:** for HTS (WGS/WES/TS) **gamma must be 1**; default 0.55 is SNP-array only.
  - **Nilsen et al. (2012)** copynumber PCF/ASPCF, BMC Genomics 13:591 (rank 1) — the penalised
    least-squares criterion `L(S|y,γ)=Σ_{I}Σ_{j}(y_j−ȳ_I)²+γ|S|`; the O(n²) DP recurrence
    `e_k=min_j(d_{jk}+e_{j−1}+γ)` (global optimum); multi-track joint cost = sum of per-track SSE, γ
    charged once per segment; conservative default γ = 40.
  - **Ross et al. (2021)** allele-specific multi-sample segmentation, Bioinformatics 37:1909 (rank 1) —
    joint logR+BAF objective with **common change points, separate per-track means**; **BAF mirroring**
    to a single track in regions of allelic imbalance.
  - **Nik-Zainal et al. (2012) / Battenberg** (rank 1/3) — two-state clonal/subclonal model
    (`nMaj1/nMin1/frac1` + `nMaj2/nMin2/frac2`, `frac1+frac2=1`); decomposition
    `n_obs = f·n₁ + (1−f)·n₂` with the two bracketing integers `n₂=⌊n_obs⌋, n₁=⌈n_obs⌉`,
    `f=(n_obs−n₂)/(n₁−n₂)`; integer `n_obs` ⇒ single clonal state.
  - **McGranahan et al. (2016)** clonal neoantigens, Science (rank 1) — expected AF / mutation copy
    number and CCF closed form (already cited at `OncologyAnalyzer.cs` line 6965).
  - **Zheng et al. (2022) PICTograph** (rank 1) — VAF generative model
    `VAF=(m·CCF·p)/(c·p+2(1−p))`; multiplicity by inversion at CCF=1, round + clamp to [1, major-CN].
  - **DeCiFering (2021, PMC8542635)** (rank 1) — CCF closed form `c=(F·v)/(ρ·M)`, `F=ρ·N_tot+2(1−ρ)`
    (confirms the implemented `EstimateCcf`).

- **Documented corner cases / failure modes:** balanced BAF=0.5 segments carry little allele-specific
  information (×0.05 weight); **non-identifiability** — the sunrise plot can show several local minima
  (2n vs 4n), ASCAT takes the global minimum; **γ must match the platform** or logR/CN are biased; the
  rounded **multiplicity is clamped to [1, major-allele CN]** (a raw <0.5 would round to 0 = non-physical
  for an observed variant).

- **Datasets (planted-truth, deterministic — synthesised by inverting the ASCAT forward model, γ=1):**
  - **Core ASCAT:** ρ₀ = 0.80, ψ₀ ∈ {2, 3}; Segment A nA=1/nB=1 (balanced, b=0.5, r=0), Segment B
    nA=2/nB=0 (copy-neutral LOH), Segment C nA=2/nB=1 (gain); a clonal m=1 mutation on a CN=2 (1+1)
    segment → CCF ≈ 1.0. Forward model
    `b=(ρ·nB+(1−ρ))/(ρ·n+2(1−ρ))`, `r=log2((ρ·n+2(1−ρ))/(ρ·ψ+2(1−ρ)))`.
  - **ASPCF breakpoint:** two-level logR track (0.0 for loci 0–9, 1.0 for 10–19, BAF 0.5), γ = 0.5 →
    2 segments, one breakpoint between index 9 and 10, means 0.0 / 1.0.
  - **Subclonal mixture (Battenberg):** ρ=1.0, ψ=2.0, γ=1; nA_obs = 0.4·2+0.6·1 = 1.4,
    nB_obs = 0.6 → states (2,0)/(1,1), fraction f ≈ 0.4; pure-clonal control nA_obs=2/nB_obs=1 →
    single state (f ≈ 0 or 1).

- **Coverage recommendations (13 items):** MUST-test allele-specific segmentation recovers planted
  breakpoints; the joint (ρ,ψ) grid recovers ρ₀=0.80 / ψ₀∈{2,3} + integer nA/nB; derived multiplicity
  = planted m; end-to-end CCF ≈ 1.0 on a clonal mutation; ASPCF recovers the planted breakpoint and its
  DP cost ≤ greedy; mirrored-BAF joint cost separates copy-neutral-LOH from balanced; large-γ→one
  segment / small-γ→each level; subclonal fit recovers f₀=0.4 with states (2,0)/(1,1); a pure-clonal
  integer segment collapses to one state. SHOULD-test that GoF at the true (ρ,ψ) beats a wrong one, and
  null/invalid-grid-bounds throw. COULD-test a balanced-only genome (all 1+1, BAF≈0.5 down-weighted).

## Deviations and assumptions

- **ASSUMPTION — germline-het-SNP BAF forward model** used **only to synthesise planted-truth inputs**
  (`b = (ρ·nB + (1−ρ)) / (ρ·(nA+nB) + 2(1−ρ))`); it is the exact algebraic inverse of the two cited
  ASCAT equations. Production consumes **measured** logR/BAF.
- **ASSUMPTION — logR normalisation reference = average sample ploidy** (`ρ·ψ + 2(1−ρ)`), so a
  genome-average segment has r = 0 (matches ASCAT's tumor-baseline-corrected logR). Synthesis-only.
- **ASSUMPTION (ASPCF) — γ is a required exposed parameter, not hard-coded.** The penalty form `+γ|S|`
  and DP recurrence are sourced verbatim; the numeric default (copynumber 40; ASCAT later 70) is
  probe-scale-specific, so tests use a γ derived from each dataset's ΔSSE (optimum provable).
- **ASSUMPTION (subclonal) — the two-state mixture uses the two bracketing integers** `⌊n_obs⌋, ⌈n_obs⌉`;
  three-or-more-population and non-adjacent-state mixtures are out of scope (documented limitation).

No source contradictions — the four algorithm stages (ASCAT core fit, ASPCF segmentation, Battenberg
subclonal, multiplicity/CCF) draw on disjoint primary literature that is mutually consistent.
</content>
