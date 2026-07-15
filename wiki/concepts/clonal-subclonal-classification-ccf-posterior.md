---
type: concept
title: "Clonal vs subclonal classification (CCF posterior, Landau/ABSOLUTE)"
tags: [oncology, algorithm]
sources:
  - docs/algorithms/Oncology/Clonal_Subclonal_Classification.md
  - docs/Evidence/ONCO-CLONAL-001-Evidence.md
source_commit: 9a7b5ef76a7b1587ce00068d78fadb60ab086bab
created: 2026-07-09
updated: 2026-07-14
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-clonal-001-evidence
      evidence: "Test Unit ID: ONCO-CLONAL-001 ... Algorithm: Clonal vs Subclonal Mutation Classification (cancer cell fraction posterior)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:cancer-cell-fraction-clonal-clustering
      source: onco-clonal-001-evidence
      evidence: "Both classify mutations clonal vs subclonal from CCF, but ONCO-CLONAL-001 uses a Bayesian per-mutation CCF posterior P(c) ∝ Binom(a|N,f(c)) with the rule clonal iff P(CCF>0.95)>0.5 (Landau 2013), whereas ONCO-CCF-001 uses point-estimate + Lloyd k-means clustering."
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:allele-specific-copy-number-ascat
      source: onco-clonal-001-evidence
      evidence: "The expected allele fraction f(c) = αMc/(2(1−α)+αq) consumes sample purity α, per-locus absolute copy number q, and multiplicity M — the outputs of the upstream ASCAT copy-number / purity fit."
      confidence: high
      status: current
---

# Clonal vs subclonal classification (CCF posterior, Landau/ABSOLUTE)

The **probabilistic clonal-structure classifier** of the Oncology family: given a somatic mutation
observed in `a` of `N` reads at a locus of absolute copy number `q` and multiplicity `M` in a sample
of purity `α`, build a **Bayesian posterior over the cancer cell fraction (CCF)** and label the
mutation **clonal** iff the posterior mass above CCF 0.95 exceeds 0.5 — the ABSOLUTE-style rule of
Landau et al. (2013). Validated under test unit **ONCO-CLONAL-001**; the literature-traced record is
[[onco-clonal-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade
correctness reference — [[scientific-rigor|research-grade]], **not for clinical or diagnostic use**.

## Relationship to CCF estimation + clustering

This unit is the **posterior-probability** approach to the same clonal/subclonal question that
[[cancer-cell-fraction-clonal-clustering]] (ONCO-CCF-001) answers with a **point estimate + Lloyd
k-means**. They are genuinely distinct algorithms — an `alternative_to` pair, not overlapping:

- **ONCO-CCF-001** — deterministic CCF point estimate `CCF = f·(ρ·N_T+2(1−ρ))/(ρ·m)`, then 1D
  k-means deconvolution; **highest-centroid cluster = clonal**.
- **ONCO-CLONAL-001** (this page) — a full **Binomial CCF posterior** per mutation, then a
  probabilistic call from its mass above 0.95. No clustering; each mutation is classified on its own
  read-count uncertainty.

Both consume the purity / copy-number / multiplicity substrate produced upstream by
[[allele-specific-copy-number-ascat]].

## 1. Expected allele fraction and the CCF posterior (Landau 2013)

Expected VAF of a mutation present in a fraction `c` of cancer cells (Landau's M=1 form):

```
f(c) = α·c / (2(1−α) + α·q)          # α = purity, c = CCF, q = per-locus absolute copy number
```

generalised by Satas et al. (2021, DeCiFering) Eq. 1 to arbitrary **multiplicity M** (SNV copies in
mutated cells):

```
f(c) = α·M·c / (2(1−α) + α·q)
```

The posterior over `c` is a Binomial likelihood with a **uniform prior**, evaluated on a **regular
100-point grid** `c ∈ [0.01, 1]` and normalised to sum to 1:

```
P(c) ∝ Binom(a | N, f(c)),    c on 100-point grid over [0.01, 1]
```

**Classification rule (verbatim, Landau 2013):** clonal iff `P(CCF > 0.95) > 0.5`, subclonal
otherwise. Note the two thresholds are distinct: **0.95** is the CCF cut, **0.5** is the posterior
probability cut.

## 2. Why probabilistic, not a point estimate

Classification is by the **posterior mass** above 0.95, not the point estimate alone: a mutation with
a point CCF near 1 but wide uncertainty (shallow coverage) can still be **subclonal**. Contrast the
worked cases — A1 (`a=N=300`, deep, P=1.0 → clonal) vs a hypothetical shallow read count with the
same mean but P(CCF>0.95) ≤ 0.5 → subclonal.

## 3. Multiplicity raises CCF for the same VAF (Satas 2021)

For `M > 1` (e.g. a mutation on both copies after copy-neutral LOH or duplication), the same VAF maps
to a **lower** CCF; ignoring M **overestimates** CCF. Worked oracle **case E**: `a=100, N=100, q=2,
M=2, ρ=1.0` → CCF mean 0.994, `P(CCF>0.95)=0.998` → clonal — the M=2 lift that the M=1 model would
place lower.

## 4. Point-estimate variant: `IdentifyClonalMutations`

A companion primitive classifies pre-computed CCF point values by the **strict** threshold **CCF >
0.95** (boundary 0.95 **excluded**). Oracle: inputs `{0.96, 0.95, 1.00, 0.50, 0.951}` →
clonal indices **{0, 2, 4}** (0.95 is not strictly greater; 0.50 is subclonal).

## Worked oracles (grid evaluation, computed independently of the implementation)

| Case | a | N | q | M | α | CCF mean | P(CCF>0.95) | Status |
|------|---|---|---|---|---|----------|-------------|--------|
| A1 | 300 | 300 | 2 | 1 | 1.0 | 0.9995 | 1.0000 | Clonal |
| B2 | 400 | 1000 | 2 | 1 | 0.8 | 0.9725 | 0.8642 | Clonal |
| C1 | 240 | 1000 | 2 | 1 | 0.8 | 0.6013 | 0.0000 | Subclonal |
| D  | 200 | 1000 | 2 | 1 | 1.0 | 0.4012 | ≈0 | Subclonal |
| E  | 100 | 100 | 2 | 2 | 1.0 | 0.9943 | 0.9980 | Clonal |

Registry invariants: **ClonalCount + SubclonalCount = total**; **ClonalFraction =
ClonalCount / total**.

## Corner cases and failure modes

- **CCF lower bound > 0** — the grid is `c ∈ [0.01, 1]`; a *detected* mutation is present in ≥1
  cancer cell, so CCF is never exactly 0.
- **Multiplicity matters** — supplying the wrong M (or defaulting to 1 on a multi-copy locus)
  overestimates CCF and can flip a subclonal call to clonal.
- **Domain guards** — purity α ∈ (0,1], valid read counts / copy number / multiplicity, CCF ∈ [0,1];
  invalid inputs throw. Empty variant set → empty calls, counts 0, ClonalFraction 0.

## Assumption: per-variant local copy number, not a genome-wide ploidy scalar

The registry stub signature was `ClassifyClonality(variants, purity, ploidy)`, but Landau's model
uses the **per-locus** absolute somatic copy number `q`, not a single genome-wide ploidy scalar. The
canonical method therefore carries `q` per variant (`ClonalityVariant.LocalCopyNumber`) and takes
`(variants, purity)`. This mirrors the prior ONCO-WGD decision (a registry scalar superseded by
per-segment data to match the authoritative definition). **Non-correctness-affecting** — API shape
only; the numerical rule and outputs are exactly Landau's.

## Implementation (ONCO-CLONAL-001 spec)

Entry points in `OncologyAnalyzer.cs`:

- `OncologyAnalyzer.ClassifyClonality(variants, purity)` — posterior-grid clonal/subclonal
  classification returning per-variant `ClonalityCall` (`Ccf`, `ProbabilityClonal`, `Status`) plus
  `ClonalCount` / `SubclonalCount` / `ClonalFraction`.
- `OncologyAnalyzer.IdentifyClonalMutations(ccfValues)` — point-estimate selection (returns 0-based
  indices with CCF > 0.95).

**Numerics:** the Binomial likelihood is evaluated in **log-space** and the constant binomial
coefficient `C(N,a)` is **omitted** — it cancels under grid normalisation. For the degenerate
all-zero posterior (e.g. `f ≈ 0` with `a > 0`), a **flat posterior** over the grid is substituted so
the result stays well-defined (classified subclonal). Purely numerical (no string search / suffix
tree). **Complexity:** `ClassifyClonality` is `O(n·G)` time, `O(n)` space (n variants, `G = 100`
grid points, constant); `IdentifyClonalMutations` is `O(m)` over m CCF values.

## Scope and limitations

Two mutually-consistent primary sources (Landau 2013, *Cell*; Satas 2021, *Cell Systems* DeCiFering)
— no contradictions. Landau supplies the posterior + 0.95/0.5 classification rule; Satas generalises
the expected-allele-fraction relation to arbitrary multiplicity M. The single flagged assumption
(per-variant `q` over a ploidy scalar) is API-shape only. **Not for clinical or diagnostic use.**
