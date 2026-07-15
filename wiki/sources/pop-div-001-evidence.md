---
type: source
title: "Evidence: POP-DIV-001 (Diversity statistics — nucleotide diversity π, Watterson's θ, Tajima's D, heterozygosity)"
tags: [validation, population-genetics]
doc_path: docs/Evidence/POP-DIV-001-Evidence.md
sources:
  - docs/Evidence/POP-DIV-001-Evidence.md
source_commit: 909848bd266e323a5133b9158dd0cd3bf668ef43
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: POP-DIV-001

The validation-evidence artifact for test unit **POP-DIV-001** — the population-genetics
**diversity-statistics** panel: nucleotide diversity **π**, Watterson's **θ_W**, **Tajima's D**,
and **heterozygosity** (observed/expected gene diversity). It is one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the formulae,
worked oracle, invariants, and corner cases are synthesized in the dedicated concept
[[genetic-diversity-statistics]]. This is the **second POP-\* (population-genetics) unit** in
the wiki, sibling of [[ancestry-estimation-admixture]] (POP-ANCESTRY-001). See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (all "exact match", no deviations):**
  - **Nei & Li (1979)**, *PNAS* 76(10):5269–73 (DOI 10.1073/pnas.76.10.5269) — **nucleotide
    diversity** `π = Σ_{i<j} d_ij / (C(n,2)·L)`: mean pairwise nucleotide differences per site over
    all `C(n,2) = n(n−1)/2` sequence pairs (`d_ij` = differences between sequences i,j; L =
    length).
  - **Watterson (1975)**, *Theoretical Population Biology* 7(2):256–276 — the **Watterson estimator**
    `θ̂_W = K / a_n` where K = number of segregating sites and `a_n = Σ_{i=1}^{n−1} 1/i` (the
    (n−1)-th harmonic number). Per-site form used in the implementation: `θ_W = S / (a_n·L)`.
  - **Tajima (1989)**, *Genetics* 123(3):585–95 — **Tajima's D** =
    `(k̂ − S/a_1) / √(e_1·S + e_2·S(S−1))`, the scaled difference between the two diversity
    estimators (π-based `k̂` vs Watterson-based `S/a_1`), expected equal under neutrality/constant
    size. **k̂ is the average pairwise difference count, NOT per-site.** Variance components:
    `a_1 = Σ 1/i`, `a_2 = Σ 1/i²`, `b_1 = (n+1)/(3(n−1))`, `b_2 = 2(n²+n+3)/(9n(n−1))`,
    `c_1 = b_1 − 1/a_1`, `c_2 = b_2 − (n+2)/(a_1·n) + a_2/a_1²`, `e_1 = c_1/a_1`,
    `e_2 = c_2/(a_1²+a_2)`. Sign interpretation: **D≈0** neutral/mutation-drift equilibrium;
    **D<0** excess rare alleles → selective sweep or population expansion; **D>0** deficit of rare
    alleles → balancing selection or population contraction.
  - **Nei (1978)**, *Genetics* 89(3):583–590 + **Wikipedia "Zygosity"** — heterozygosity. Diploid
    observed heterozygosity `H_o = (# heterozygous individuals)/n` needs genotypes; for **haploid
    sequence** data the unit uses Nei's **unbiased gene diversity** as the H_o analogue,
    `H_obs = n/(n−1)·(1/L)·Σ_pos(1 − Σ_i p_i²)`, which is mathematically **π for haploid data**.
    Expected heterozygosity = **basic gene diversity** `H_exp = (1/L)·Σ_pos(1 − Σ_i p_i²)`
    (`H_e = 1 − Σ f_i²` per site). Relationship: `H_obs = n/(n−1)·H_exp` (Nei 1978 bias
    correction).

- **Datasets (documented oracle — Wikipedia Tajima's D worked example):** n=5, L=20, five
  0/1 sequences (Y,A,B,C,D). S = 4 segregating sites (positions 3,7,13,19);
  `a_1 = 1+1/2+1/3+1/4 = 25/12 ≈ 2.0833`; pairwise-difference total = 20 over 10 comparisons →
  **k̂ = 2.0**, π = k̂/L = 0.1, `S/a_1 ≈ 1.92`, θ_W(per-site) `= 4/(2.0833·20) ≈ 0.096`,
  numerator `d = k̂ − S/a_1 = 0.08`; with `a_2 ≈ 1.4236`, `b_1 = 0.5`, `b_2 ≈ 0.3667`,
  `c_1 = 0.02`, `c_2 ≈ 0.0227`, `e_1 ≈ 0.0096`, `e_2 ≈ 0.00393`, Var `≈ 0.0856` →
  **D ≈ 0.273** (tests TD-C01/TD-C02).

## Deviations and assumptions

**None** — every formula is an exact match to its published source (π/Nei-Li 1979, θ_W/Watterson
1975, D/Tajima 1989, gene diversity + Nei-1978 bias correction). Notable API/contract points:
the method signature `CalculateTajimasD(averagePairwiseDifferences, segregatingSites, sampleSize)`
takes **k̂ (unnormalized pairwise-difference average), NOT per-site π** as its first argument, and
computes the Watterson estimate `S/a_1` internally. Edge cases: n=0 → zeros; n=1 → π=0, θ
undefined→0; n=2 → π defined but Tajima's D undefined→0 (Tajima 1989 requires n≥3); S=0 /
monomorphic / all-identical → all metrics 0; Var≤0 → D=0 (numerical guard). No source
contradictions.
