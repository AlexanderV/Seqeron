---
type: concept
title: "Effective Number of Codons (ENC / Nc)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/CODON-ENC-001-Evidence.md
source_commit: 5bc4ea5003342f5c5c657d68183edfae40fba29a
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: codon-enc-001-evidence
      evidence: "Test Unit ID: CODON-ENC-001 ... Algorithm: Effective Number of Codons (ENC / Nc), Wright 1990"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:codon-adaptation-index
      source: codon-enc-001-evidence
      evidence: "ENC (Wright 1990) summarizes a gene's overall codon bias as a reference-free number in [20,61]; CAI (Sharp & Li 1987) summarizes it as a reference-based score in [0,1] — two different single-number codon-bias measures"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:relative-synonymous-codon-usage
      source: codon-enc-001-evidence
      evidence: "Codon homozygosity F̂ = (n·Σp_i² − 1)/(n − 1) is built from the same synonymous-codon frequencies p_i = n_i/n that RSCU formalizes"
      confidence: high
      status: current
---

# Effective Number of Codons (ENC / Nc)

A single **reference-free codon-usage-bias measure** for a gene, defined by **Wright (1990)** as
an analogue of the population-genetics "effective number of alleles". Where
[[codon-adaptation-index|CAI]] scores a gene against a reference set, **Nc quantifies how many
codons are *effectively* in use** across all amino acids — from **20** (each amino acid uses a
single codon: maximal bias) up to **61** (every synonymous codon used equally: no bias). Lower Nc
= stronger bias. Validated as [[codon-enc-001-evidence|CODON-ENC-001]] against Fuglsang (2004,
2006), which reproduce Wright's original equations verbatim; see [[test-unit-registry]] for how
the unit is tracked.

## The measure

Built up from **codon homozygosity** `F̂` per amino acid, then aggregated by degeneracy class.

**Per-amino-acid homozygosity (Eq. 1):** for an amino acid with `n` total codons in the gene and
synonymous-codon frequencies `p_i = n_i/n`:

    F̂ = ( n · Σ_i p_i² − 1 ) / ( n − 1 )

This is the sampling-**without-replacement** estimator (Fuglsang 2006 concludes it is the superior
one). The **effective codons for that amino acid** is its reciprocal (Eq. 2): `N̂c(aa) = 1 / F̂`.

**Gene-level aggregation (Eq. 3)** averages `F̂` within each **degeneracy class** of the standard
genetic code and sums the reciprocals, weighted by how many amino acids fall in each class:

    N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆

The constants encode the standard-code degeneracy partition of the 20 amino acids: **2** singlets
(Met/ATG, Trp/TGG — always exactly one codon), **9** doublets (2-fold), **1** triplet (isoleucine,
3-fold), **5** quartets (4-fold), **3** sextets (6-fold). Stop codons are excluded.

**Unbiased limit** (each synonymous codon equally likely): `F̂₂ = 0.5`, `F̂₃ = 1/3`, `F̂₄ = 0.25`,
`F̂₆ = 1/6`, giving `2 + 9/0.5 + 1/(1/3) + 5/0.25 + 3/(1/6) = 2 + 18 + 3 + 20 + 18 = 61`.

## Corner-case rules (from Wright / Fuglsang)

- **Amino acid with n ≤ 1.** `F̂` is undefined (the `n−1` denominator is 0). Nc requires ≥ 2 codons
  for each represented amino acid; a class member with no estimable `F̂` is dropped and the class
  average is taken over the remaining members (Eq. 4 within-class averaging, e.g.
  `F̂₄ = (F̂_pro + F̂_gly + F̂_ala + F̂_val)/4` when threonine is absent).
- **Empty 3-fold class.** Isoleucine is the *only* 3-fold amino acid, so if it is absent the class
  average is undefined; Wright's explicit fallback is **`F̂₃ = (F̂₂ + F̂₄)/2`** (Eq. 5a).
- **Upper-bound re-adjustment.** Eq. 3 can yield `N̂c > 61` for nearly-uniform short genes; Wright
  prescribes **re-adjusting the result down to 61**.
- **Per-aa overshoot.** For small `n`, an evenly-used amino acid can give `N̂c(aa)` above its
  degeneracy (e.g. an evenly-split 2-fold codon → 3 > 2); the gene-level value is still capped at 61.
- **Calculability.** Nc is defined only when at least one amino acid per degeneracy class has ≥ 2
  codons; Wright used *E. coli* K12 as the reference organism in the original paper.

## Deviation / assumption (one, from the artifact)

- **Lower clamp at 20.** Wright/Fuglsang state Nc *approaches* 20 under extreme bias and explicitly
  prescribe re-adjusting **down to 61** at the top — they give no hard clamp at 20. 20 is the
  structural minimum (every class collapses to one codon ⇒ each `N̂c(aa)=1`), so a defensive
  `Math.Max(20, …)` cannot lower a legitimately-computed value. Treated as a defensive bound, not a
  Wright-prescribed parameter.

## Worked oracles

- **Fully unbiased gene** (equal codon counts throughout, large `n`) → Nc **61**.
- **Maximally biased gene** (one codon per amino acid) → Nc **20**.
- **Two-fold hand example (Eq. 1):** Phe TTT×3, TTC×1 (n=4, p=¾,¼): `Σp² = 10/16`,
  `F̂ = (4·0.625 − 1)/3 = 0.5`, `Nc(Phe) = 2`. Even split TTT×2/TTC×2: `F̂ = 1/3`, `Nc(Phe) = 3`
  (illustrates the per-aa overshoot).
- **Fuglsang "no-bias-discrepancy" simulation** → Nc **40.5** (per-aa Nc: 2-fold 1.5 ×9, Ile 2.0,
  4-fold 2.5 ×5, 6-fold 3.5 ×3, singlets 1.0 ×2).
- Invariant **20 ≤ Nc ≤ 61** for arbitrary inputs; isoleucine-absent gene uses the `(F̂₂+F̂₄)/2`
  fallback; null → ArgumentNullException, empty/whitespace → 0.

## Place in the codon-usage family

Nc is the **reference-free** counterpart of the reference-based [[codon-adaptation-index|CAI]]: both
reduce a whole gene's codon bias to one number, but Nc needs no reference gene set (it measures
evenness of synonymous-codon usage intrinsically, 20–61) while CAI measures adaptation toward a
reference set (0–1). Both are built from the same synonymous-family codon frequencies that
[[relative-synonymous-codon-usage|RSCU]] formalizes (`F̂` is a function of the `p_i`). Other
siblings still in `docs/Evidence/` include rare-codon analysis, codon optimization, and raw
codon-frequency/usage tables. See [[algorithm-validation-evidence]] for the shared evidence pattern.
</content>
</invoke>
