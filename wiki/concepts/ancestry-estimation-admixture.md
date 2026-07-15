---
type: concept
title: "Ancestry estimation (supervised/projection ADMIXTURE)"
tags: [population-genetics, algorithm]
mcp_tools:
  - estimate_ancestry
sources:
  - docs/Evidence/POP-ANCESTRY-001-Evidence.md
  - docs/algorithms/Population_Genetics/Ancestry_Estimation.md
source_commit: c691bc46fb22a8ae15f49a0fcac5c96548a50b4e
created: 2026-07-10
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pop-ancestry-001-evidence
      evidence: "Test Unit ID: POP-ANCESTRY-001 ... Algorithm: Ancestry Estimation (supervised / projection ADMIXTURE)"
      confidence: high
      status: current
---

# Ancestry estimation (supervised/projection ADMIXTURE)

Estimate an individual's **ancestry proportions** — what fraction of their genome derives from each
of K reference populations — from unlinked SNP genotypes. This is the **first population-genetics
`POP-*` unit** in the wiki and the anchor for that family. It implements the **supervised /
projection** mode of **ADMIXTURE** (Alexander, Novembre & Lange 2009; Alexander & Lange 2011): the
reference-panel allele frequencies **F** are treated as **fixed, known constants** and only the
ancestry vector **Q** for each query individual is estimated by EM. Validated under test unit
**POP-ANCESTRY-001**; the literature-traced record is [[pop-ancestry-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
artifact pattern.

It sits alongside — but is distinct from — the wiki's **phylogenetics** `PHYLO-*` family. Those
methods place *populations/taxa on a tree* from distances ([[evolutionary-distance-matrix]],
[[distance-based-tree-construction]]); ancestry estimation instead decomposes *one admixed
individual* into a **mixture of pre-defined ancestral sources** — a likelihood/EM mixture-weight
problem, not a tree.

## Model (Alexander et al. 2009)

Genotypes are the **allele-1 copy counts** `g_ij ∈ {0, 1, 2}` for individual `i` at SNP `j` (2 =
homozygous 1/1, 1 = heterozygous 1/2, 0 = homozygous 2/2). Two parameter matrices:

- **Q = {q_ik}** (I×K) — fraction of individual `i`'s genome from population `k`.
- **F = {f_kj}** (K×J) — frequency of allele 1 at SNP `j` in population `k`. These per-population
  allele frequencies are exactly the output of the foundational
  [[allele-genotype-frequencies]] primitive (POP-FREQ-001), here taken as fixed/known input.

Assuming Hardy–Weinberg within populations and independence across SNPs, the **log-likelihood
(Eq. 2)** is

```
L(Q,F) = Σ_i Σ_j { g_ij·ln(Σ_k q_ik f_kj) + (2−g_ij)·ln(Σ_k q_ik (1−f_kj)) }
```

subject to `0 ≤ f_kj ≤ 1`, `q_ik ≥ 0`, and **`Σ_k q_ik = 1`** (proportions sum to one). In this
unit **F is fixed** (learned from labelled reference panels, each representing a 100%-single-ancestry
source), so the optimization reduces to estimating Q per individual.

## EM update (FRAPPE, Eq. 4)

Ancestry proportions are refined by the FRAPPE EM iteration:

```
q_ik^{n+1} = (1/2J) · Σ_j [ g_ij · a^n_ijk + (2 − g_ij) · b^n_ijk ]
a^n_ijk = q^n_ik f^n_kj      / (Σ_m q^n_im f^n_mj)          # posterior that an allele-1 copy came from k
b^n_ijk = q^n_ik (1 − f^n_kj) / (Σ_m q^n_im (1 − f^n_mj))   # posterior that an allele-2 copy came from k
```

Each iteration is an **EM ascent step**: the log-likelihood is non-decreasing, and the estimate is
automatically renormalized (the `1/2J` averaging over `2J` allele copies keeps `Σ_k q_ik = 1`).
**Convergence (Eq. 5):** stop when `L(Q^{n+1}) − L(Q^n) < ε`; ADMIXTURE's default is `ε = 10⁻⁴`.
The unstructured joint (Q+F) EM is **O(IJK²)** per iteration [Alexander et al. 2009]; here **F is fixed**,
so there is no per-population panel-update term and the per-individual cost drops to **O(iterations·J·K)**
(one pass over SNPs × populations, `O(K)` working space), giving **O(I·iterations·J·K)** for the full call.

## Worked oracles

- **Symmetric two-population panel, diagnostic individual.** K = 2 (A, B), J = 2,
  `f_A = [0.8, 0.2]`, `f_B = [0.2, 0.8]`; genotype `g = [2, 0]`; start `q⁰ = (0.5, 0.5)`. One EM
  iteration gives **`q = (0.8, 0.2)` exactly** (L = −1.5426…); then (0.9412, 0.0588), (0.9846,
  0.0154), … converging to **`(1.0, 0.0)`**. L strictly increases each step; q sums to 1 every
  iteration.
- **Single-SNP closed forms** (J = 1, K = 2, `q⁰ = (0.5,0.5)`, `f = [[0.9],[0.1]]`):
  `g = 2 → (0.9, 0.1)`; **`g = 1 → (0.5, 0.5)`** (a heterozygote on a symmetric panel is
  uninformative — the a and b contributions cancel); `g = 0 → (0.1, 0.9)`.
- **Uninformative identical panels** (`f_A = f_B` at every SNP): a uniform q is a **fixed point** —
  `(0.5,0.5) → (0.5,0.5)` for any genotype, because the a/b numerators and denominators are
  proportional.

## Invariants and edge cases

- **INV — simplex:** every returned individual satisfies `Σ_k q_ik = 1`, `q_ik ≥ 0`.
- **INV — EM monotone ascent:** the Eq. 2 log-likelihood is non-decreasing across iterations (basis
  of the Eq. 5 stopping rule); the diagnostic individual converges to its true source population.
- **Uninformative fixed point:** identical reference panels leave a uniform individual uniform.
- **Missing / malformed genotypes:** values outside {0,1,2} are treated as **missing** and skipped
  for that SNP (no Eq. 2 term). An individual whose genotype vector length ≠ J is skipped. An
  individual with **all** SNPs missing is returned at the **uniform prior** `q = (1/K, …)` (no
  informative term ever updates q).
- **Empty inputs:** empty individuals or empty reference panels → empty result.

## Implementation (Seqeron)

`PopulationGeneticsAnalyzer.EstimateAncestry(individuals, referencePops, maxIterations)`
(`Seqeron.Genomics.Population`) is the public entry point; it iterates individuals and yields one
`AncestryProportion { IndividualId, Proportions }`, where `Proportions` is a dictionary keyed by
reference-population id that sums to 1 (INV-01). Per-individual work is delegated to the private
`EstimateIndividualAncestry` (FRAPPE EM, Eq. 4, F fixed) with the Eq. 2 objective computed by
`AncestryLogLikelihood` for the Eq. 5 stopping test.

- **Inputs:** `individuals` as `(IndividualId, IReadOnlyList<int> Genotypes)` with each genotype in
  `{0,1,2}` and length = panel SNP count; `referencePops` as `(PopulationId, IReadOnlyList<double>
  AlleleFrequencies)` with each frequency in `[0,1]` and the **same SNP order** across all panels.
- **`maxIterations` default 100** — an EM budget layered on top of the ε = 10⁻⁴ early stop; both
  reach the same maximum. The convergence tolerance is the constant `AncestryLogLikelihoodTolerance`.
- **Validation is total, not throwing:** empty inputs → empty result; length-mismatch individual
  skipped; out-of-range genotype treated as missing. No exceptions for these documented input classes.
- The dividing-by-`2J` step (J = *informative* SNP count for that individual) keeps `Σ_k q_ik = 1`
  without a separate renormalization. Not a substring search, so the repository suffix tree is N/A.

Related population-genetics unit: [[population-differentiation-fst]] (F_Statistics).

## Identifiability and scope

The likelihood is **invariant under permutations of the K population labels** (Eq. 2), so an
*unsupervised* fit has at least **K! equivalent global maxima**. This unit sidesteps that: with
**fixed, labelled** reference panels the labels are **pinned**, so the returned proportions are
directly interpretable — but only **relative to the supplied reference labels**. Two documented
API-shape [[research-grade-limitations|research-grade]] assumptions: (1) a fixed `maxIterations`
budget augments the ε = 10⁻⁴ convergence rule (the EM reaches the same maximum either way); (2)
missing-genotype handling is skip-the-site. This is the **supervised/projection** ("estimate Q given
fixed F", ADMIXTURE manual §2.10/§2.14) task only — it does **not** jointly learn F, run the
unsupervised block-relaxation/quasi-Newton ADMIXTURE optimizer, choose K, or model linkage. No
source contradictions: Alexander et al. (2009), Alexander & Lange (2011), and the ADMIXTURE 1.4
manual are mutually consistent, the manual's supervised/projection modes being the applied
specialization of the 2009 model.
