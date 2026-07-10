---
type: source
title: "Evidence: POP-HW-001 (Hardy-Weinberg equilibrium chi-square test)"
tags: [validation, population-genetics]
doc_path: docs/Evidence/POP-HW-001-Evidence.md
sources:
  - docs/Evidence/POP-HW-001-Evidence.md
source_commit: 0ddacceea3a398cc68e027111682aa4ce726fcb7
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: POP-HW-001

The validation-evidence artifact for test unit **POP-HW-001** — the **Hardy-Weinberg equilibrium
(HWE) chi-square test**: does a biallelic locus's observed genotype counts match the
`p²/2pq/q²` proportions expected under random mating? It is one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the formulae, worked
oracles, invariants, and corner cases are synthesized in the dedicated concept
[[hardy-weinberg-equilibrium-test]]. This is a population-genetics `POP-*` unit that **consumes**
the allele/genotype frequencies produced by [[allele-genotype-frequencies]] (POP-FREQ-001 —
which explicitly leaves the HWE test to this unit) and sits in the family anchored by
[[ancestry-estimation-admixture]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (all "exact match", no algorithm deviations):**
  - **Wikipedia "Hardy-Weinberg principle"** (after Hardy 1908 / Weinberg 1908) — the constancy
    of allele/genotype frequencies across generations absent evolutionary forces, and the
    expected genotype frequencies `f(AA)=p²`, `f(Aa)=2pq`, `f(aa)=q²`; the Ford scarlet-tiger-moth
    worked example; deviation causes; and the note that in real genotype data a deviation may
    signal genotyping error.
  - **Wikipedia "Chi-squared test"** — Pearson goodness-of-fit `χ² = Σ (O−E)²/E`.
  - **Hardy (1908) Science** and **Weinberg (1908)** — the original principle.
  - **Ford (1971) Ecological Genetics** — the moth dataset (via Wikipedia).
  - **Emigh (1980) Biometrics** and **Wigginton et al. (2005) AJHG** — chi-square-test comparison
    and the **exact test** method (cited as the alternative to the chi-square approximation).

- **Formulae recorded:** allele frequencies from genotype counts `p = (2·n_AA + n_Aa)/(2n)`,
  `q = 1 − p`; expected counts `E = {p²n, 2pqn, q²n}`; chi-square `Σ (n − E)²/E` over the three
  genotype classes; **1 degree of freedom** (`#genotypes − #alleles = 3 − 2`); the p-value from the
  chi-square CDF (df=1); default significance `α = 0.05` (critical value 3.841).

- **Implementation notes:** `TestHardyWeinberg` in `PopulationGeneticsAnalyzer.cs` computes allele
  frequencies, expected counts, and the chi-square statistic, then derives the p-value via a
  **chi-square CDF using a lower-incomplete-gamma approximation**. Returns `InEquilibrium`,
  `ChiSquare`, `PValue`.

- **Oracles & datasets:** Ford's moth (1469, 138, 5) → `p≈0.954`, `χ²≈0.83`, InEquilibrium = true
  (HW-E01); perfect HWE (25, 50, 25) → `χ²=0` (HW-E02); excess heterozygotes (10, 80, 10) →
  `χ²=36 ≫ 3.84`, InEquilibrium = false (HW-E03); zero samples (0,0,0) → InEquilibrium = true,
  PValue = 1 (HW-E04); fixed/monomorphic allele (100, 0, 0) → InEquilibrium = true (HW-E05);
  all heterozygotes (0, 100, 0) → InEquilibrium = false (HW-E06).

## Deviations and assumptions

**None** as an algorithm deviation — every formula is an exact match to the Hardy-Weinberg / Pearson
chi-square definitions. Documented contract behaviours: **zero sample size** returns
`InEquilibrium = true, ChiSquare = 0, PValue = 1` (no evidence against H₀ → maximal compatibility
with the null, a direct consequence of the hypothesis-testing framework, not an ad-hoc assumption);
a genotype class with **expected count 0** has its chi-square term **skipped** (division-by-zero
protection, consistent with the conditional Wikipedia formula). Scope is the **biallelic chi-square
goodness-of-fit** test only — the **exact test** (Wigginton 2005) and multiallelic loci are noted as
out of scope. No source contradictions; Open Questions: none.
