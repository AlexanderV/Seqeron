---
type: source
title: "Evidence: META-ALPHA-001 (alpha diversity — Shannon/Simpson/invSimpson/Chao1/Pielou/richness)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-ALPHA-001-Evidence.md
sources:
  - docs/Evidence/META-ALPHA-001-Evidence.md
source_commit: 88b3a1e12a0b76ef17934a3d6d3c12f96a1fe058
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-ALPHA-001

The validation-evidence artifact for test unit **META-ALPHA-001** — **alpha diversity**, the
within-sample diversity summary computed by `MetagenomicsAnalyzer.CalculateAlphaDiversity`. This is
the **first ingested unit of the Metagenomics family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The method is synthesized in its own
concept, [[alpha-diversity]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (all mutually consistent, no contradictions):**
  - **Wikipedia — Diversity index / Alpha diversity / Species richness / Species evenness**
    (accessed 2026-02-02) — the Shannon, Simpson, inverse-Simpson, Pielou and richness formulas, and
    the definition of alpha diversity as **mean species diversity at the local (single-sample) scale**.
  - **Shannon (1948)** *A Mathematical Theory of Communication* — Shannon entropy `H' = −Σ pᵢ ln(pᵢ)`.
  - **Simpson (1949)** *Measurement of Diversity* — Simpson concentration `λ = Σ pᵢ²`.
  - **Chao (1984)** *Non-parametric estimation of the number of classes in a population* — the Chao1
    richness estimator (needs singleton/doubleton counts).
  - **Hill (1973)** *Diversity and evenness: a unifying notation* — inverse Simpson `1/λ` = true
    diversity of order 2 (effective number of species).
  - **Pielou (1966)** — evenness `J = H / ln(S)`.
  - *(the algorithm doc adds **Whittaker 1960** for the α/β/γ diversity framework that places alpha as
    within-sample diversity.)*

- **Extracted formulas & properties:** `H = 0` for one species / `H = ln(R)` at perfect evenness (nats,
  natural log); `λ = 1` for one species / `λ = 1/R` at perfect evenness; `1/λ` = effective species
  count; `J ∈ [0,1]`, `J = 1` at perfect evenness, undefined (→0 by convention) when `S ≤ 1`;
  `Chao1 = S_obs + f₁²/(2f₂)`, bias-corrected `S_obs + f₁(f₁−1)/2` when `f₂ = 0`, `= S_obs` when
  `f₁ = 0`, and `= S_obs` for proportional (non-integer) data where singletons/doubletons are undefined.

- **Documented edge / corner cases:** empty input → all metrics `0`; single species → `H=0, λ=1,
  J` undefined/0; perfect evenness → `H=ln(S), λ=1/S, J=1`; zero abundances filtered out (ln(0)
  undefined); unnormalized abundances normalized to sum 1 before calculation; all-zero / null input →
  zero result; negative abundances = invalid/undefined.

- **Datasets (documented oracles):**
  - Single species (p=1.0) → H 0, λ 1.0.
  - Two equal species (0.5, 0.5) → H ln(2) ≈ 0.693, λ 0.5, 1/λ 2, J 1.0.
  - Four equal species (0.25 each) → H ln(4) ≈ 1.386, λ 0.25, 1/λ 4, J 1.0.
  - Uneven (0.9, 0.1) → H ≈ 0.325, J ≈ 0.469 (dominance lowers evenness).

## Implementation notes (from the Evidence file)

Filters non-positive abundances, normalizes to proportions summing to 1, computes Shannon with the
natural logarithm (`Math.Log`), Simpson `Σpᵢ²`, inverse Simpson `1/λ`, Pielou `H/ln(S)` for `S>1`
(else `0`, the standard ecological convention since `ln(1)=0`), and Chao1 with the bias-corrected
`f₂=0` branch for integer count data / `S_obs` fallback for proportional data. Null/empty → all metrics 0.

## Deviations and assumptions

- **Evidence file:** *"None. All formulas match external sources exactly."*
- **Nuance vs the algorithm doc** (`docs/algorithms/Metagenomics/Alpha_Diversity.md`, §5.4): that doc
  records **one accepted deviation** — Chao1 falls back to `ObservedSpecies` for **non-integer /
  proportional** abundance input (so callers passing relative abundances receive observed richness, not
  a singleton/doubleton unseen-richness estimate). This is a data-type gate, not a formula change; the
  Evidence file's "formulas match exactly" and the doc's Chao1-fallback note are consistent. Captured on
  [[alpha-diversity]].

No source contradictions — the encyclopedic and primary-literature formulas (Shannon 1948, Simpson
1949, Hill 1973, Chao 1984, Pielou 1966) are the standard, mutually consistent definitions.
</content>
</invoke>
