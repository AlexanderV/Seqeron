---
type: concept
title: "Population differentiation (Fst, F-statistics, pairwise Fst)"
tags: [population-genetics, algorithm]
mcp_tools:
  - fst
  - pairwise_fst
sources:
  - docs/Evidence/POP-FST-001-Evidence.md
source_commit: 6a9852103155b627075f1a105de26fac5b97f70a
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pop-fst-001-evidence
      evidence: "Test Unit ID: POP-FST-001 ... Algorithm: Fst (Fixation Index), F-Statistics (Fis, Fit, Fst)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:allele-genotype-frequencies
      source: pop-fst-001-evidence
      evidence: "CalculateFst consumes per-population allele frequencies p1, p2 (with sizes n1, n2): pBar = (n1*p1 + n2*p2)/(n1+n2); variance = (n1*(p1-pBar)^2 + n2*(p2-pBar)^2)/(n1+n2). These frequencies are the output of the allele/genotype-frequency primitive."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:genetic-diversity-statistics
      source: pop-fst-001-evidence
      evidence: "F-statistics use expected heterozygosity: Fis = 1 - Hi/Hs, Fit = 1 - Hi/Ht, Fst = 1 - Hs/Ht; the Wright denominator pBar(1-pBar) is expected heterozygosity under Hardy-Weinberg — the same gene-diversity quantity POP-DIV-001 computes as H = 1 - Sum p_i^2."
      confidence: high
      status: current
---

# Population differentiation (Fst, F-statistics, pairwise Fst)

Quantify **how genetically differentiated populations are** — the fraction of total genetic variation
that is due to *between-population* structure rather than *within-population* diversity. This is a
population-genetics `POP-*` unit (**POP-FST-001**) in the family anchored by
[[ancestry-estimation-admixture]]. It is genuinely distinct from its POP siblings: it does not
*count* alleles ([[allele-genotype-frequencies]]) nor summarise the variation *inside one sample*
([[genetic-diversity-statistics]]) — it compares **two or more populations' allele frequencies** and
returns a differentiation scalar. It **consumes** the per-population allele frequencies produced by
the [[allele-genotype-frequencies]] primitive. Validated under test unit **POP-FST-001**; the
literature-traced record is [[pop-fst-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## Fst — Wright's variance-based fixation index

The implementation uses **Wright's (1965) variance definition directly** (Wikipedia *Fixation
index*):

```
Fst = sigma_S^2 / (pBar * (1 - pBar))
```

where `sigma_S^2` is the variance in allele frequency **among** subpopulations (weighted by
population size) and `pBar*(1-pBar)` is the expected heterozygosity of the total population under
Hardy–Weinberg. For two populations with sizes `n1, n2` and allele frequencies `p1, p2`:

```
pBar     = (n1*p1 + n2*p2) / (n1 + n2)                          # size-weighted mean
sigma_S^2 = (n1*(p1-pBar)^2 + n2*(p2-pBar)^2) / (n1 + n2)
het       = pBar * (1 - pBar)
Fst       = sigma_S^2 / het
```

**Multi-locus aggregation is ratio-of-sums**, not a mean of per-locus ratios:
`Fst = (Sum_l sigma_S,l^2) / (Sum_l pBar_l(1-pBar_l))`.

**Important distinction:** this computes the **population parameter** from *known* allele
frequencies. It is **not** the Weir & Cockerham (1984) θ estimator, which uses ANOVA variance
components with finite-sample (unequal-sample-size) bias correction — the implementation applies
**no** Weir–Cockerham correction. This is a documented modelling choice, not a bug.

## F-statistics — the heterozygosity partition

`CalculateFStatistics` uses the heterozygosity-based definitions (Wikipedia *F-statistics*), given
observed individual heterozygosity `Hi`, expected within-subpopulation `Hs`, and expected total
`Ht`:

```
Fis = 1 - Hi/Hs     # inbreeding of Individual within Subpopulation
Fit = 1 - Hi/Ht     # inbreeding of Individual within Total
Fst = 1 - Hs/Ht     # effect of Subpopulation structure within Total
```

The **partition identity holds exactly** (algebraic, not approximate): `(1 - Fit) = (1 - Fis)(1 -
Fst)`, because `(Hi/Hs)(Hs/Ht) = Hi/Ht`.

## Invariants and value ranges

- **Fst range `0 ≤ Fst ≤ 1`.** `Fst = 0` = complete panmixia (no differentiation, freely
  interbreeding); `Fst = 1` = complete differentiation (fixed differences, no shared diversity).
- **Fis, Fit range `-1 … 1`** — Fis can be **negative** under an excess of heterozygotes.
- **Pairwise Fst matrix:** diagonal `Fst(i,i) = 0`, symmetric `Fst(i,j) = Fst(j,i)`, non-negative
  for distinct populations.
- **Not a metric:** Fst does **not** satisfy the triangle inequality, so it is not a distance in the
  mathematical sense (Wikipedia).

## Interpretation scale and reference oracles

Wright / Hartl & Clark qualitative bands: `0–0.05` little, `0.05–0.15` moderate, `0.15–0.25` great,
`> 0.25` very great differentiation. Literature reference values used as oracles: human population
pairs from Cavalli-Sforza (1994) — Danish–English 0.0021, Lapps–Sardinians 0.0667, Mbuti–Papua New
Guineans 0.4573, 42-population mean 0.1338; continental HapMap pairs from Elhaik (2012) — Europe vs
Sub-Saharan Africa 0.153, Europe vs East Asia 0.111, within-continent `< 0.01`.

## Edge cases

- **Identical populations / all-monomorphic sites → Fst = 0** (no variance).
- **Fixed differences `p1=1, p2=0` → Fst = 1.0 exactly** (pBar=0.5, variance=0.25, het=0.25,
  ratio=1). Verifiable by hand.
- **Denominator zero** (`pBar = 0` or `1`, e.g. both populations fixed for the same allele, or empty
  populations) → **return 0** (design decision for the 0/0 case, not undefined behaviour).
- **Unequal sample sizes** are handled by size-weighting `c_i = n_i/N` (Wright 1965), *not* by
  Weir–Cockerham correction.

## Scope

Faithful implementation of Wright's variance-based Fst and the heterozygosity-based F-statistics
partition (exact match to Wikipedia *Fixation index* / *F-statistics*, Wright 1965). It computes the
**population parameter from known allele frequencies**; it does **not** implement the Weir &
Cockerham (1984) θ estimator, per-locus significance testing, or bootstrap confidence intervals. No
source contradictions; the Evidence file records no deviations (Open Questions: none).
