---
type: concept
title: "Linkage disequilibrium (D, D', r², haplotype blocks)"
tags: [population-genetics, algorithm]
sources:
  - docs/Evidence/POP-LD-001.md
source_commit: fadbea3029500764efb2211347df8b83ad90d190
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pop-ld-001-evidence
      evidence: "Test Unit ID: POP-LD-001 ... Methods: CalculateLD, FindHaplotypeBlocks (population-genetics unit)."
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:allele-genotype-frequencies
      source: pop-ld-001-evidence
      evidence: "D = p_AB − p_A·p_B and r² = D²/(p_A·q_A·p_B·q_B) are built from the per-allele frequencies p_A, p_B (q = 1 − p) and the two-locus haplotype frequency p_AB — the allele/haplotype frequencies produced by the frequency primitive."
      confidence: high
      status: current
---

# Linkage disequilibrium (D, D', r², haplotype blocks)

Quantify the **non-random association between alleles at two different loci** — the signature that
two genetic sites are inherited together more (or less) often than chance predicts. This is a
population-genetics `POP-*` unit (**POP-LD-001**) in the family anchored by
[[ancestry-estimation-admixture]]. It is genuinely distinct from its POP siblings: it is a **pairwise
inter-locus** quantity, not a per-locus count ([[allele-genotype-frequencies]]), a within-sample
diversity summary ([[genetic-diversity-statistics]]), a between-population differentiation scalar
([[population-differentiation-fst]]), or a single-locus goodness-of-fit test
([[hardy-weinberg-equilibrium-test]]). It **consumes** the allele/haplotype frequencies produced by
the [[allele-genotype-frequencies]] primitive. Validated under test unit **POP-LD-001**; the
literature-traced record is [[pop-ld-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The three LD measures

For two biallelic loci with alleles A/a (frequencies `p_A`, `q_A = 1 − p_A`) and B/b (`p_B`,
`q_B = 1 − p_B`), and observed two-locus haplotype frequency `p_AB` (Wikipedia *Linkage
disequilibrium*):

```
D  = p_AB − p_A · p_B                        # raw coefficient (Robbins 1918)
```

`D` is the excess of the AB haplotype over the independence expectation `p_A·p_B`. Its magnitude is
hard to compare across loci because its range depends on the allele frequencies, so two
normalizations are used:

```
D' = |D| / D_max,  clamped to [0, 1]         # Lewontin 1964
  D_max = min(p_A·q_B, q_A·p_B)   when D < 0
  D_max = min(p_A·p_B, q_A·q_B)   when D ≥ 0

r² = D² / (p_A · q_A · p_B · q_B)            # Hill & Robertson 1968
```

`D'` measures how close the pair is to the theoretical maximum LD given the marginal frequencies
(range 0–1 as `|D'|`); `r²` is the squared correlation between the two loci (range 0–1) and is the
measure that drives GWAS tagging and block detection.

## r² without phase — the diploid-covariance shortcut

The implementation of `CalculateLD` does **not** require phased haplotypes. Using the diploid-frequency
result (Wright 1933, Hill & Robertson 1968) that the diploid correlation `R_AB` equals the
haplotype-level `r_AB`, it computes `r²` as the **squared Pearson correlation of the 0/1/2 genotype
dosage vectors** `X₁, X₂` (0/1/2 = count of alternate alleles):

```
r² = Cov(X₁, X₂)² / (Var(X₁) · Var(X₂))
Cov = Σ (X₁ᵢ − μ₁)(X₂ᵢ − μ₂) / n            Var = Σ (Xᵢ − μ)² / n
```

`D` is then recovered from the same covariance, because in the 0/1/2 encoding
`Cov_diploid(X₁, X₂) = 2D`, so **`D = Cov(X₁, X₂) / 2`**; `D'` follows the Lewontin branch above.
Consequences: identical genotype vectors → `r² = 1` (perfect LD); balanced/independent vectors →
`r² = 0`; **no haplotype phasing / EM is needed**.

## Haplotype-block detection

`FindHaplotypeBlocks` uses a **simplified adjacent-pair** version of Gabriel et al. (2002): scan
consecutive variants, and where the pairwise `r²` between neighbours is `≥ threshold` (default
**0.7**; strong-LD conventions use 0.7–0.8) they extend a block. A block requires **≥ 2 variants**.

## Invariants and value ranges

- **`0 ≤ r² ≤ 1`** and **`0 ≤ |D'| ≤ 1`** (clamped).
- **Empty input → `r² = 0`, `D' = 0`** (no data).
- **Monomorphic locus** (zero variance ⇒ zero denominator) → **`r² = 0`**, division guarded — not
  `NaN`/undefined. All-identical genotypes likewise → `r² = 0` (no polymorphism to correlate).
- **Distance and variant IDs are preserved** from the input record.
- **Block invariants:** `Start ≤ End`; ≥ 2 variants per block; blocks are **non-overlapping** and
  **ordered by position**; adjacent variants within a block satisfy `r² ≥ threshold`.

## Edge cases

- Single genotype pair → computable but statistically unreliable (`n = 1`).
- Perfect correlation → `r² = 1`, `D' = 1`; **perfect anti-correlation also → `r² = 1`, `D' = 1`**
  (`r²` and `|D'|` are sign-blind).
- Single variant → no block possible; all variants in strong LD → one spanning block; no pair above
  threshold → no blocks; strong pairs separated by a weak gap → multiple separate blocks.

## Scope

Faithful implementation of the textbook two-locus LD measures and adjacent-pair block detection
(exact match to Wikipedia *Linkage disequilibrium* / *Haplotype block*, Lewontin 1964, Hill &
Robertson 1968, Gabriel et al. 2002). It computes `r²`/`D'` between **two biallelic loci** and
detects blocks along a variant run; it does **not** build a full pairwise LD matrix, phase
haplotypes (EM), fit an LD-decay curve, or reproduce the exact Gabriel confidence-interval bounds. No
source contradictions; the Evidence file records no deviations (Open Questions: none).
