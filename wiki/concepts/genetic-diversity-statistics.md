---
type: concept
title: "Genetic diversity statistics (π, Watterson's θ, Tajima's D, heterozygosity)"
tags: [population-genetics, algorithm]
mcp_tools:
  - diversity_statistics
  - tajimas_d
sources:
  - docs/Evidence/POP-DIV-001-Evidence.md
  - docs/algorithms/Population_Genetics/Diversity_Statistics.md
source_commit: e1e652ea6b5955cb98edcb4a870d0f0c4b24a90b
created: 2026-07-10
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pop-div-001-evidence
      evidence: "Test Unit ID: POP-DIV-001 ... Algorithm: Nucleotide Diversity, Watterson's Theta, Tajima's D"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:ancestry-estimation-admixture
      source: pop-div-001-evidence
      evidence: "POP-DIV-001 is the second population-genetics POP-* unit, sibling of POP-ANCESTRY-001; both operate on population sequence/genotype data but measure different quantities (within-population diversity vs ancestry proportions)."
      confidence: high
      status: current
---

# Genetic diversity statistics (π, Watterson's θ, Tajima's D, heterozygosity)

Summarize the **amount and pattern of sequence variation within a single population sample** as four
classical population-genetics scalars. This is the **second population-genetics `POP-*` unit**
(POP-DIV-001), sibling of the family anchor [[ancestry-estimation-admixture]] (POP-ANCESTRY-001).
It is genuinely distinct from ancestry estimation: ancestry decomposes *one admixed individual* into
a mixture of pre-defined ancestral sources (an EM mixture-weight problem), whereas these statistics
describe *how variable a whole sample is* and whether that variation departs from neutral
expectation. Validated under test unit **POP-DIV-001**; the literature-traced record is
[[pop-div-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The four statistics

Given `n` aligned sequences of length `L`, `S` = number of segregating (polymorphic) sites, `d_ij` =
pairwise differences between sequences i and j, and `p_i` = allele frequencies at a site:

- **Nucleotide diversity π** (Nei & Li 1979) — mean pairwise differences **per site**:
  `π = Σ_{i<j} d_ij / (C(n,2)·L)`, with `C(n,2) = n(n−1)/2` pairwise comparisons.
- **Watterson's θ_W** (Watterson 1975) — diversity estimated from the **count of segregating
  sites**: `θ̂_W = S / a_n`, harmonic-number normalizer `a_n = Σ_{i=1}^{n−1} 1/i`; per-site form
  `θ_W = S / (a_n·L)`. π and θ_W estimate the same population parameter θ = 4Nμ under neutrality.
- **Tajima's D** (Tajima 1989) — the **normalized difference** between the two estimators, testing
  the neutral hypothesis:
  `D = (k̂ − S/a_1) / √(e_1·S + e_2·S(S−1))`.
  `k̂` is the **average pairwise-difference count (NOT per-site π)**; `S/a_1` is the Watterson count
  estimate. Variance terms `a_1 = Σ 1/i`, `a_2 = Σ 1/i²`, `b_1 = (n+1)/(3(n−1))`,
  `b_2 = 2(n²+n+3)/(9n(n−1))`, `c_1 = b_1 − 1/a_1`, `c_2 = b_2 − (n+2)/(a_1·n) + a_2/a_1²`,
  `e_1 = c_1/a_1`, `e_2 = c_2/(a_1²+a_2)`.
- **Heterozygosity / gene diversity** (Nei 1978) — from the per-site allele frequencies `p_i`
  produced by the [[allele-genotype-frequencies]] primitive (POP-FREQ-001), expected
  `H_exp = (1/L)·Σ_pos(1 − Σ_i p_i²)`
  (per-site `1 − Σ f_i²`); unbiased observed analogue for **haploid** data
  `H_obs = n/(n−1)·H_exp`, which equals π for haploid sequences. (True diploid `H_o` = fraction of
  heterozygous individuals; unavailable without genotypes.)

## Tajima's D sign interpretation

| D | Cause (allele-frequency spectrum) | Biological reading |
|---|---|---|
| ≈ 0 | `k̂ ≈ S/a_1` | neutral, mutation–drift equilibrium |
| < 0 | excess of rare alleles (`k̂ < S/a_1`) | selective sweep / population expansion |
| > 0 | deficit of rare alleles (`k̂ > S/a_1`) | balancing selection / population contraction |

## Worked oracle (Wikipedia Tajima's D example)

n=5, L=20, sequences Y/A/B/C/D; **S = 4** (positions 3,7,13,19); `a_1 = 25/12 ≈ 2.0833`. Pairwise
differences total 20 over 10 comparisons → **k̂ = 2.0**, π = 0.1, `S/a_1 ≈ 1.92`, numerator
`d = 0.08`, θ_W(per-site) `≈ 0.096`. With `a_2 ≈ 1.4236`, `b_1 = 0.5`, `b_2 ≈ 0.3667`, `c_1 = 0.02`,
`c_2 ≈ 0.0227`, `e_1 ≈ 0.0096`, `e_2 ≈ 0.00393`, Var `≈ 0.0856` → **D ≈ 0.273** (tests TD-C01/TD-C02).

## Invariants and edge cases

- **Sample-size guards:** `n=0` → all zeros; `n=1` → π=0, θ undefined→0; `n=2` → π defined but
  Tajima's D **undefined→0** (Tajima 1989 requires n≥3).
- **No polymorphism:** `S=0` / monomorphic / all-identical sequences → π=θ=D=H=0.
- **Numerical guard:** variance ≤ 0 → D=0.
- **API-shape note:** `CalculateTajimasD(averagePairwiseDifferences, segregatingSites, sampleSize)`
  expects **k̂**, not per-site π; it derives `S/a_1` internally. Passing π (per-site) here would be a
  units error.

## Implementation shape (POP-DIV-001)

`PopulationGeneticsAnalyzer` exposes four public entry points plus two private helpers; the
integrated `CalculateDiversityStatistics(IEnumerable<IReadOnlyList<char>>)` orchestrates the rest and
returns the `DiversityStatistics` record.

- **π** via `CalculateNucleotideDiversity` — explicit `O(n²·m)` pairwise mismatch counting (`n`
  sequences, `m` sites), `O(n)` space. **θ_W** via `CalculateWattersonTheta(S, n, L)` — `O(n)`
  harmonic-number loop. **Tajima's D** via `CalculateTajimasD(averagePairwiseDifferences, S, n)` —
  `O(n)`.
- **Segregating-site counting is anchored to the first sequence:** a site is called polymorphic when
  any sequence's character at that position differs from the *first* sequence's character there (not
  a full all-pairs distinctness test). Assumes equal-length aligned input.
- **Integrated k̂ conversion:** `CalculateDiversityStatistics` computes per-site π first, then feeds
  Tajima's D the unnormalized `k̂ = π · L` — deliberately reconciling the per-site-π vs count-k̂
  confusion (INV/units note above).
- **API-name gotcha:** the record field is spelled **`SegregratingSites`** (source-code typo) and
  that spelling is part of the public API, preserved in source and tests for compatibility. The two
  heterozygosity fields carry repository-specific names: `HeterozygosityExpected` = basic
  `1 − Σ p_i²`, `HeterozygosityObserved` = Nei's `n/(n−1)` bias-corrected per-site gene diversity
  (equals π for haploid data) — **not** a diploid observed-heterozygosity count.
- **No explicit equal-length validation pass:** the implementation indexes every sequence at the
  first sequence's positions, so callers must pre-align/pre-normalize; shorter sequences yield invalid
  indexing or uninterpretable comparisons (accepted deviation). No p-value/significance layer for
  Tajima's D is provided.

## Scope

Faithful, deviation-free implementation of the four classical estimators (all "exact match" to Nei
& Li 1979, Watterson 1975, Tajima 1989, Nei 1978). It computes **summary statistics** from an
existing alignment/allele-count table; it does **not** call variants, phase haplotypes, correct for
recombination, or run other neutrality tests (Fu & Li's D, Fay & Wu's H). No source contradictions.
