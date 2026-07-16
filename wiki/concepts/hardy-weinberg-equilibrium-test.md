---
type: concept
title: "Hardy-Weinberg equilibrium test (chi-square)"
tags: [population-genetics, algorithm]
mcp_tools:
  - hardy_weinberg_test
sources:
  - docs/Evidence/POP-HW-001-Evidence.md
  - docs/algorithms/Population_Genetics/Hardy_Weinberg_Test.md
source_commit: 758c875b14f9b85d6a9290da5a0fb1c56f8d6478
created: 2026-07-10
updated: 2026-07-16
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pop-hw-001-evidence
      evidence: "Test Unit ID: POP-HW-001 ... Algorithm: Hardy-Weinberg Equilibrium Chi-Square Test."
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:allele-genotype-frequencies
      source: pop-hw-001-evidence
      evidence: "The test first derives allele frequencies from observed genotype counts (p = (2*n_AA + n_Aa)/(2n), q = 1 - p) and expected counts E = {p^2*n, 2pq*n, q^2*n} before computing chi-square; it consumes the allele/genotype-frequency primitive (POP-FREQ-001, whose scope explicitly excludes the HWE test)."
      confidence: high
      status: current
---

# Hardy-Weinberg equilibrium test (chi-square)

Given a **biallelic** locus's observed diploid genotype counts, decide whether the population is in
**Hardy-Weinberg equilibrium (HWE)** — i.e. whether the genotypes occur in the `p²/2pq/q²`
proportions expected under random mating with no other evolutionary forces. This is a
population-genetics `POP-*` unit (**POP-HW-001**) in the family anchored by
[[ancestry-estimation-admixture]]. It is genuinely distinct from its POP siblings: it does not merely
*count* alleles — it **consumes** those counts from [[allele-genotype-frequencies]]
(POP-FREQ-001, whose scope explicitly leaves the HWE test to this unit) and turns them into a
**statistical hypothesis test** with a p-value. Validated under test unit **POP-HW-001**; the
literature-traced record is [[pop-hw-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The principle

For a biallelic locus with major-allele frequency `p` and minor-allele frequency `q = 1 − p`, the
**Hardy-Weinberg principle** (Hardy 1908, Weinberg 1908) predicts constant genotype frequencies
across generations, absent selection/mutation/migration/drift/non-random mating:

```
f(AA) = p²        f(Aa) = 2pq        f(aa) = q²
```

## The chi-square goodness-of-fit test

1. **Estimate allele frequencies** from the observed genotype counts `n_AA`, `n_Aa`, `n_aa`
   (`n = n_AA + n_Aa + n_aa`):

   ```
   p = (2·n_AA + n_Aa) / (2n)        q = 1 − p
   ```

2. **Expected counts** under HWE:

   ```
   E(AA) = p²·n     E(Aa) = 2pq·n     E(aa) = q²·n
   ```

3. **Pearson chi-square** statistic over the three genotype classes:

   ```
   χ² = Σ (O − E)² / E
      = (n_AA − E_AA)²/E_AA + (n_Aa − E_Aa)²/E_Aa + (n_aa − E_aa)²/E_aa
   ```

4. **Degrees of freedom = 1** (`#genotypes − #alleles = 3 − 2` for the biallelic case, because the
   allele frequency `p` was estimated from the same data). The **p-value** comes from the chi-square
   CDF at df = 1 (Seqeron uses a **lower-incomplete-gamma approximation**). At the default
   significance **α = 0.05** the critical value is **3.841**: `χ² ≥ 3.841 ⇒ reject HWE`.

`TestHardyWeinberg` returns `InEquilibrium` (the accept/reject decision), `ChiSquare`, and `PValue`.

## Critical values and decision

| df | α = 0.10 | α = 0.05 | α = 0.01 |
|----|----------|----------|----------|
| 1  | 2.706    | 3.841    | 6.635    |

- `p ≥ α` → **fail to reject H₀**: consistent with HWE (`InEquilibrium = true`).
- `p < α` → **reject H₀**: significant deviation from HWE (`InEquilibrium = false`).

## Worked oracles

- **Ford's scarlet tiger moth** (1469, 138, 5): `p ≈ 0.954`, expected `(1467.4, 141.2, 3.4)`,
  `χ² ≈ 0.83` → p > 0.05, **in equilibrium**.
- **Perfect HWE** (25, 50, 25 with `p = 0.5`): observed = expected → `χ² = 0`, in equilibrium.
- **Excess heterozygotes** (10, 80, 10): `p = 0.5`, expected `(25, 50, 25)` →
  `χ² = 9 + 18 + 9 = 36 ≫ 3.84`, **out of equilibrium** (signals heterozygote advantage /
  assortative mating — impossible under random mating).

## Invariants and edge cases

- **Zero sample size** (`n = 0`): returns `InEquilibrium = true, χ² = 0, PValue = 1`. No data means
  no evidence against H₀; PValue = 1 is maximal compatibility with the null — a direct consequence
  of the hypothesis-testing framework, not an ad-hoc convention.
- **Fixed / monomorphic allele** (e.g. all AA): `p = 1, q = 0`, expected = observed, `χ² = 0` →
  in equilibrium.
- **Division-by-zero protection**: a genotype class with **expected count 0** has its chi-square
  term **skipped** (consistent with the conditional Wikipedia formula), so a monomorphic or
  single-individual sample does not blow up.
- **Deviation causes** (why a locus fails the test): non-random mating / inbreeding (excess
  homozygotes), selection, mutation, migration/gene flow, small-population drift — and, in real
  genotype data, **genotyping error** (a common bioinformatics QC use of the test).

## Implementation shape (POP-HW-001)

`PopulationGeneticsAnalyzer.TestHardyWeinberg(string variantId, int observedAA, int observedAa,
int observedaa, double significanceLevel = 0.05)` is the single entry point. It returns a result
record echoing the inputs (`VariantId`, `ObservedAA/Aa/aa`) alongside the derived
`ExpectedAA/Aa/aa` (= `p²n`, `2pqn`, `q²n`), `ChiSquare`, `PValue`, and the boolean
`InEquilibrium` (strictly `PValue ≥ significanceLevel`). The p-value is `1 − ChiSquareCDF(χ², 1)`,
where the private `ChiSquareCDF(double, int)` routes through a `RegularizedGammaP(double, double)`
regularized lower-incomplete-gamma helper (this is the "lower-incomplete-gamma approximation"
named above). Cost is **O(1) time / O(1) space** — constant-time arithmetic on one genotype triple.

**Accepted deviation (§5.4):** the method takes `int` counts directly and does **not** validate that
genotype counts are non-negative, nor that `significanceLevel` lies in `[0, 1]`. Negative counts
produce uninterpretable frequencies/expectations; the guard is left to the caller. This is a
documented, accepted simplification, not a bug.

## Scope

Faithful implementation of the biallelic **chi-square goodness-of-fit** HWE test (exact match to the
Hardy-Weinberg / Pearson definitions, Ford 1971 oracle). It does **not** implement the **exact test**
(Wigginton et al. 2005, the small-sample / rare-allele alternative) nor multiallelic (>2-allele)
loci. No source contradictions — the algorithm is fully determined by the sources (Open Questions:
none).
