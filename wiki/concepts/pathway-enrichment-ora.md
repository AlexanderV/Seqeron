---
type: concept
title: "Pathway enrichment / over-representation analysis (ORA) — hypergeometric right-tail test"
tags: [metagenomics, algorithm]
sources:
  - docs/Evidence/META-PATHWAY-001-Evidence.md
  - docs/algorithms/Metagenomics/Pathway_Enrichment_ORA.md
source_commit: 14005a6134e0cf637135d54041289269cd71c467
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: meta-pathway-001-evidence
      evidence: "Test Unit ID: META-PATHWAY-001, Algorithm: Metabolic Pathway Enrichment (ORA via the hypergeometric test), Method FindPathwayEnrichment / HypergeometricUpperTail"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:functional-prediction
      source: meta-pathway-001-evidence
      evidence: "FindPathwayEnrichment / HypergeometricUpperTail is validated standalone here (META-PATHWAY-001) and also exercised as component B of functional prediction (META-FUNC-001); same shared machinery"
      confidence: high
      status: current
---

# Pathway enrichment / over-representation analysis (ORA)

**Over-representation analysis (ORA)** answers "given a query gene list, which metabolic **pathways** (or
GO terms / gene sets) are represented among them **more than chance would predict**?" Seqeron implements
it as a single specification-driven, deterministic method — `MetagenomicsAnalyzer.FindPathwayEnrichment`,
scoring each candidate pathway with the **hypergeometric right-tail** (upper-tail) test
`HypergeometricUpperTail`. This is the classic GO::TermFinder / clusterProfiler enrichment statistic.

This concept **owns** the ORA / hypergeometric machinery, validated under test unit
**META-PATHWAY-001** ([[meta-pathway-001-evidence]]). The *same* method also appears as **component B** of
[[functional-prediction]] (META-FUNC-001), where it scores pathways over a metagenome's homology-transferred
functions — that page links here rather than re-deriving the statistic. [[test-unit-registry]] tracks the
unit; [[algorithm-validation-evidence]] describes the artifact pattern. The numerical core is **exact**
with respect to the cited formulas.

## The hypergeometric right-tail test

Draw a **query** of `n` genes without replacement from a **background universe** of `N` genes that
contains `M` members of a pathway; the overlap `X` (query genes in the pathway) is hypergeometric. The
over-representation p-value is the probability of observing `x` **or more** overlaps (Boyle et al. 2004
GO::TermFinder; PNNL ORA §8.2):

```
P(X ≥ x) = 1 − Σ_{i=0}^{x−1}  C(M,i) · C(N−M, n−i) / C(N, n)
```

- **Symbols (verbatim from the sources):** `N` = total background genes; `M` = pathway/term size (genes
  annotated to the node); `n` = query size; `x` (`= k`) = the overlap. The reference R computation is
  `phyper(q = x−1, m = M, n = N−M, k = n, lower.tail = FALSE)` — the `q = x−1` shift is exactly what makes
  `lower.tail = FALSE` yield P(X ≥ x) rather than P(X > x).
- **One-sided / over-representation only:** this is the **upper tail** ("x or more"); under-representation
  is a separate lower-tail test. The **hypergeometric** (sampling *without* replacement) is used rather
  than the binomial because it is exact for a finite universe.
- **Numerics:** the tail is summed in **log-space** via a Lanczos **log-Gamma** to avoid factorial
  overflow (validated to `N = 8000`); the final p is **clamped to `[0,1]`**; pathways are returned
  **sorted ascending by p-value** (Boyle's P-value term ranking).

## Background / query handling

- An explicit `backgroundGenes` argument sets `N`. When it is null/empty, `N` **defaults to the union of
  all pathway members** (the documented, caller-overridable convenience — see ASM below).
- Pathway members are **intersected with the background** before counting, so `M` and `x` are measured
  against the actually-sampled universe; the **query is always unioned into the background**.
- **Not implemented:** multiple-testing correction (BH / Bonferroni / FDR) — apply your own over the
  returned p-values.

## Invariants and edge cases

- **p ∈ [0,1]** always.
- **`p = 1`** when `x = 0` (empty upper-sum), or on degenerate margins `M = 0` or `n = 0` (no success can
  be drawn).
- **Symmetry:** the hypergeometric is invariant under swapping the two marked groups, so `(M, n)` and
  `(n, M)` give the **same** p-value — e.g. `(400, 100)` and `(100, 400)` both yield 7.884747×10⁻⁸.
- **Without-replacement feasibility:** terms with `i > M` or `(n − i) > (N − M)` are infeasible
  (`C(·,·) = 0`) and drop out of the sum.
- Results are sorted ascending by p; the reported overlap/size fields echo the inputs.

Worked oracles (from [[meta-pathway-001-evidence]]): PNNL §8.2 `N=8000, M=100, n=400, x=20` →
`P(X ≥ 20) = 7.88 × 10⁻⁸`; exact rational cases `N10/M5/n5/x5 → 1/252`, `N4/M2/n2/x1 → 5/6`,
`N10/M5/n5/x0 → 1`, `N10/M5/n5/x1 → 251/252`.

## Scope and limitations

A [[research-grade-limitations|research-grade]] enrichment statistic: the hypergeometric right-tail
p-value is **exact** and source-backed, but ORA is a *threshold*-based method — it needs a pre-selected
query list and treats every gene inside it equally, ignoring effect-size ranking (the domain of GSEA) and
inter-gene dependence within a pathway. No FDR/multiple-testing correction is applied, and no gene-set
database is bundled (pathways are caller-supplied). The single accepted deviation (background defaulting
to the pathway-member union when no universe is given) affects `N` on the no-background path only, not the
formula. Though registered under metagenomics via `FindPathwayEnrichment`, the statistic is generic (the
sources are GO::TermFinder and proteomics ORA). No source contradictions.
</content>
