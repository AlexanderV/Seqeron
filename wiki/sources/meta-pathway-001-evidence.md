---
type: source
title: "Evidence: META-PATHWAY-001 (pathway enrichment / over-representation analysis via the hypergeometric right-tail test)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-PATHWAY-001-Evidence.md
sources:
  - docs/Evidence/META-PATHWAY-001-Evidence.md
source_commit: 14005a6134e0cf637135d54041289269cd71c467
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-PATHWAY-001

The validation-evidence artifact for test unit **META-PATHWAY-001** — **metabolic pathway enrichment /
over-representation analysis (ORA)** by the **hypergeometric right-tail test**
(`MetagenomicsAnalyzer.FindPathwayEnrichment` / `HypergeometricUpperTail`). This is the **dedicated unit
for the ORA/hypergeometric machinery** that also appears as **component B** of
[[functional-prediction]] (META-FUNC-001). The method is synthesized in its own concept,
[[pathway-enrichment-ora]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Boyle et al. 2004 — GO::TermFinder** (*Bioinformatics* 20(18):3710–3715, PMC3037731; authority rank
    1 paper / 3 reference implementation) — the over-representation p-value as the **right tail of the
    hypergeometric distribution** `P = 1 − Σ_{i=0}^{k−1} C(M,i)·C(N−M,n−i)/C(N,n)`, with the verbatim
    symbol definitions (*N* = background gene count, *M* = pathway/term size, *n* = query size, *k* =
    overlap), the one-sided **P(X ≥ k)** ("k or more") semantics, and the rationale that the
    hypergeometric (sampling **without** replacement) is more accurate than the binomial.
  - **PNNL — Proteomics Data Analysis in R/Bioconductor §8.2 (ORA)** (rank 3; documents the canonical
    R/Bioconductor `phyper` reference) — the identical formula
    `P(X ≥ x) = 1 − P(X ≤ x−1) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M,n−i)/C(N,n)`, and the reference call
    `phyper(q = x−1, m = M, n = N−M, k = n, lower.tail = FALSE)` (the `q = x−1` shift is what makes
    `lower.tail = FALSE` compute P(X ≥ x) rather than P(X > x)).

- **Documented corner / failure cases:**
  - **x = 0 (no overlap):** empty upper-sum ⇒ `P = 1 − 0 = 1` (largest p-value, no over-representation
    possible).
  - **Degenerate margins (N = 0, M = 0, or n = 0):** no successes can be drawn ⇒ `P(X ≥ x) = 1`.
  - **Without-replacement constraint:** terms with `i > M` or `(n − i) > (N − M)` are infeasible
    (`C(·,·) = 0`) and contribute 0.

- **Datasets (documented oracles):**
  - **PNNL §8.2 worked example:** `N = 8000, M = 100, n = 400, x = 20` → **P(X ≥ 20) = 7.88 × 10⁻⁸**
    (7.884747×10⁻⁸); the hypergeometric is **symmetric** under swapping the two marked groups, so
    `(M, n) = (400, 100)` and `(100, 400)` give the same p-value (both orientations independently computed).
  - **Exact hand-derived rational cases** (exact binomial coefficients): all-query-in-pathway
    `N10/M5/n5/x5` → `1/252 = 0.003968…`; partial overlap `N4/M2/n2/x1` → `5/6 = 0.8333…`; no-overlap
    `N10/M5/n5/x0` → `1`; at-least-one `N10/M5/n5/x1` → `251/252 = 0.99603…`.

## Recommended test coverage (from the Evidence file)

MUST: PNNL 7.88×10⁻⁸; exact 1/252 (all overlap) and 5/6 (partial overlap); x = 0 → p = 1; result fields
match inputs and results sorted **ascending by p-value** (Boyle's P-value ranking of terms). SHOULD:
**symmetry** — (M, n) and (n, M) give the same p; caller-supplied vs defaulted background changes N and
hence p. COULD: many pathways processed independently in one call.

## Deviations and assumptions

- **ASSUMPTION — background defaulting.** When the caller supplies no background universe, the
  implementation uses the **union of all pathway members (plus the query)** as the background. The cited
  sources require a *caller-defined* background and do not prescribe a default, so this is a documented,
  caller-overridable convenience. It does **not** change the formula when a background is supplied;
  it is correctness-affecting only for the no-background path, which is tested against a hand-derived
  expectation. No source contradictions.

## Relationship to META-FUNC-001

`FindPathwayEnrichment` / `HypergeometricUpperTail` is the same machinery that META-FUNC-001 exercises as
**component B** of [[functional-prediction]]. META-PATHWAY-001 is its **dedicated** unit: it frames ORA
standalone with its own GO::TermFinder + PNNL sources, adds the **M↔n symmetry** invariant and the exact
rational oracles, and validates `docs/algorithms/Metagenomics/Pathway_Enrichment_ORA.md`. The
[[pathway-enrichment-ora]] concept **owns** the method; [[functional-prediction]] links to it.
