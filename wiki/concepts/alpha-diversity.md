---
type: concept
title: "Alpha diversity (within-sample diversity indices)"
tags: [metagenomics, algorithm]
sources:
  - docs/Evidence/META-ALPHA-001-Evidence.md
  - docs/algorithms/Metagenomics/Alpha_Diversity.md
source_commit: 88b3a1e12a0b76ef17934a3d6d3c12f96a1fe058
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: meta-alpha-001-evidence
      evidence: "Test Unit ID: META-ALPHA-001 ... Area: Metagenomics ... Canonical Methods: MetagenomicsAnalyzer.CalculateAlphaDiversity"
      confidence: high
      status: current
---

# Alpha diversity (within-sample diversity indices)

**Alpha diversity** summarizes the diversity *within a single sample or site* by combining **richness**
(how many taxa are present) and **evenness** (how evenly their abundances are distributed) into one set
of scalar statistics. It is the local term of Whittaker's α/β/γ framework (Whittaker 1960): alpha =
within-sample, beta = between-sample turnover, gamma = regional. This is the **first ingested unit of
the Metagenomics family** and the **anchor** for the metagenomics diversity topic; its between-sample
sibling [[beta-diversity]] (Bray-Curtis / Jaccard turnover) is a future unit. Validated under test unit
**META-ALPHA-001**; the record is [[meta-alpha-001-evidence]], [[test-unit-registry]] tracks the unit,
and [[algorithm-validation-evidence]] describes the artifact pattern.

The single entry point `MetagenomicsAnalyzer.CalculateAlphaDiversity(IReadOnlyDictionary<string,double>)`
takes one **taxon → abundance** map and returns an `AlphaDiversity` record carrying **six** metrics:
`ObservedSpecies`, `ShannonIndex`, `SimpsonIndex`, `InverseSimpson`, `PielouEvenness`, `Chao1Estimate`.
Abundances may be **counts or proportions** — the method internally normalizes the positive values to
proportions summing to 1 before computing Shannon/Simpson, and filters out any abundance `≤ 0` (ln(0) is
undefined). `O(n)` time/space in the number of taxa; deterministic.

## The six metrics

For proportions `pᵢ` (positive abundances normalized to sum 1) over `S` observed taxa:

```
ObservedSpecies  S_obs = |{ i : pᵢ > 0 }|                 (richness)
ShannonIndex     H     = −Σ pᵢ ln(pᵢ)                     (natural log ⇒ nats; Shannon 1948)
SimpsonIndex     λ     = Σ pᵢ²                            (concentration; Simpson 1949)
InverseSimpson   1/λ                                      (effective species = Hill order 2; Hill 1973)
PielouEvenness   J     = H / ln(S)   for S > 1, else 0    (Pielou 1966)
Chao1Estimate    Ŝ     = S_obs + f₁²/(2·f₂)               (unseen-richness estimator; Chao 1984)
```

- **Shannon `H`** is the entropy of a randomly drawn individual's identity: `0` for a single species (no
  uncertainty), maximal `ln(S)` at perfect evenness. The repository uses `Math.Log`, so `H` is reported
  in **nats** (natural log), not bits or `log₁₀`.
- **Simpson `λ`** is the probability that two randomly drawn individuals are the *same* species: `1` for
  one species, `1/S` at perfect evenness; it *rises* with dominance (so it is a concentration, not a
  diversity — higher λ = less diverse).
- **Inverse Simpson `1/λ`** flips λ into an effective-number-of-species reading (the true diversity of
  order 2): `2` for two equal species, `S` at perfect evenness.
- **Pielou `J`** rescales `H` by its maximum `ln(S)` so `J ∈ [0,1]`, `= 1` at perfect evenness, `→ 0`
  under dominance. It is **undefined when `S ≤ 1`** (`ln(1) = 0` ⇒ division by zero); the implementation
  returns `0` in that branch, the standard ecological convention.

## Chao1 (the one data-type-dependent branch)

`Chao1` estimates *unseen* richness from rare taxa:

```
f₂ > 0 :  Ŝ = S_obs + f₁² / (2·f₂)              (standard Chao1)
f₂ = 0 :  Ŝ = S_obs + f₁·(f₁ − 1) / 2           (bias-corrected form)
f₁ = 0 :  Ŝ = S_obs                             (no singletons ⇒ nothing to extrapolate)
```

where `f₁` = number of **singleton** taxa (abundance exactly 1) and `f₂` = **doubleton** taxa
(abundance exactly 2). Singletons/doubletons only exist for **integer count data**. The implementation
therefore gates Chao1 on the input *looking like counts*: if every positive abundance is effectively
integer-valued it applies the estimator (with the `f₂ = 0` bias-corrected branch); if any positive
abundance is **non-integer / proportional**, it returns `ObservedSpecies` for `Chao1Estimate` instead.
This is the single accepted **deviation** (algorithm doc §5.4): relative-abundance callers get observed
richness, not an extrapolated estimate — the output shape is preserved but no unseen-richness correction
is applied. The other five metrics are computed identically for counts and proportions (they are
normalized either way).

## Invariants and edge cases

- **INV-01:** `ObservedSpecies` = count of taxa with strictly positive abundance (non-positive filtered).
- **INV-02:** single species → `H = 0`.
- **INV-03:** `λ = 1` for a single species, decreasing as abundance spreads.
- **INV-04:** `InverseSimpson = 1/λ` when `λ > 0`.
- **INV-05:** `PielouEvenness = 0` when `ObservedSpecies ≤ 1` (undefined branch → 0).
- **Empty / null input → every field `0`** (all-zero `AlphaDiversity` record); likewise when no positive
  abundances remain after filtering.
- **Zero abundances** are dropped before any calculation; **negative** abundances are invalid/undefined.

Worked oracles (from [[meta-alpha-001-evidence]]): single species → `H 0 / λ 1`; `(0.5, 0.5)` →
`H ln2 ≈ 0.693 / λ 0.5 / 1/λ 2 / J 1`; `(0.25, 0.25, 0.25, 0.25)` → `H ln4 ≈ 1.386 / λ 0.25 / 1/λ 4 /
J 1`; uneven `(0.9, 0.1)` → `H ≈ 0.325 / J ≈ 0.469` (dominance lowers evenness).

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the six standard alpha-diversity
indices over one abundance profile. All formulas match their primary sources exactly (Shannon 1948,
Simpson 1949, Hill 1973, Chao 1984, Pielou 1966); the only implementation-specific behaviour is the
Chao1 count-vs-proportion gate and the `S ≤ 1 → J = 0` convention. It computes summary statistics over a
*supplied* abundance map — taxonomic classification / profiling that produces that map, and the
between-sample turnover metrics, are separate units ([[beta-diversity]] and the rest of the metagenomics
family). No rarefaction, no confidence intervals, no Hill-number series beyond order 2. No source
contradictions.
</content>
