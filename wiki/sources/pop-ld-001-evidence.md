---
type: source
title: "Evidence: POP-LD-001 (Linkage disequilibrium — D, D', r², haplotype blocks)"
tags: [validation, population-genetics]
doc_path: docs/Evidence/POP-LD-001.md
sources:
  - docs/Evidence/POP-LD-001.md
source_commit: fadbea3029500764efb2211347df8b83ad90d190
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: POP-LD-001

The validation-evidence artifact for test unit **POP-LD-001** — **linkage disequilibrium** between
two loci (`CalculateLD`) and **haplotype-block detection** (`FindHaplotypeBlocks`). It is one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
the formulae, worked oracles, invariants, and corner cases are synthesized in the dedicated concept
[[linkage-disequilibrium]]. This is a population-genetics `POP-*` unit that **consumes** the
allele/haplotype frequencies produced by [[allele-genotype-frequencies]] (POP-FREQ-001) and sits in
the family anchored by [[ancestry-estimation-admixture]]. See [[test-unit-registry]] for how units
are tracked.

## What this file records

- **Online sources (all "exact match", no deviations):**
  - **Wikipedia "Linkage disequilibrium"** — formal `D = p_AB − p_A·p_B`, Lewontin's `D'`
    normalization, Hill & Robertson `r²`, decay of LD with recombination, and the diploid-frequency
    result that the diploid correlation `R_AB` equals the haplotype-level `r_AB` (Wright 1933, Hill &
    Robertson 1968) — the basis for computing `r²` **without haplotype phase**.
  - **Wikipedia "Haplotype block"** — block definition on a high-LD threshold; Gabriel et al. (2002)
    and Patil et al. (2001) block-definition methods.
  - **Lewontin, R.C. (1964)** *Genetics* 49(1):49–67 — introduced the `D'` normalization.
  - **Hill, W.G. & Robertson, A. (1968)** *TAG* 38(6):226–231 — introduced the `r²` correlation
    measure.
  - **Gabriel, S.B. et al. (2002)** *Science* 296(5576):2225–2229 — haplotype-block detection from
    LD thresholds.

- **Implementation choices (documented, not deviations):**
  - `CalculateLD` computes **`r²` as the squared Pearson correlation of the 0/1/2 genotype dosage
    vectors** — `r² = Cov(X₁,X₂)² / (Var(X₁)·Var(X₂))` — which is mathematically equivalent to the
    haplotype `r²` and requires **no phase information**. `D` is recovered from the diploid
    covariance: `Cov_diploid = 2D` in the 0/1/2 encoding, so `D = Cov(X₁,X₂)/2`; `D'` follows
    Lewontin (1964) with `D_max` branched on the sign of `D`, then `|D'|` clamped to `[0,1]`.
  - `FindHaplotypeBlocks` uses a **simplified adjacent-pair** Gabriel rule: consecutive variants with
    pairwise `r² ≥ threshold` (default **0.7**) form a block; a block needs **≥ 2 variants**.

- **Oracles & datasets:** perfect LD (identical genotype vectors, e.g. (0,0),(1,1),(2,2) repeated) →
  `r² ≈ 1`, high `D'`; no LD (independent/balanced genotypes) → `r² ≈ 0`, `D' ≈ 0`; perfect
  anti-correlation → `r² = 1`, `D' = 1`. Block oracles: single variant → no block; two variants in
  high LD → one block; two variants in low LD → no block; all-strong-LD run → one spanning block;
  non-contiguous strong pairs → multiple separate blocks.

## Deviations and assumptions

**None** as an algorithm deviation — every formula is an exact match to Wikipedia (Linkage
disequilibrium / Haplotype block), Lewontin (1964), Hill & Robertson (1968), and Gabriel et al.
(2002). Contract behaviours (design decisions, not undefined behaviour): **empty genotypes →
`r² = 0`, `D' = 0`**; **monomorphic locus** (zero variance, denominator 0) → **`r² = 0`** (division
guarded); all-identical genotypes → `r² = 0` (no polymorphism to correlate); distance and variant
IDs are **preserved** from input. Block invariants: `Start ≤ End`, ≥ 2 variants, non-overlapping,
position-ordered. Scope is `r²`/`D'` between two biallelic loci + adjacent-pair block detection only
— no full pairwise LD matrix, no haplotype phasing/EM, no LD-decay curve fitting, no exact Gabriel
confidence-interval bounds. No source contradictions; Open Questions: none.
