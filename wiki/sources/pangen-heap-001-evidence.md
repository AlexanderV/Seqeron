---
type: source
title: "Evidence: PANGEN-HEAP-001 (Pan-genome growth model — Heaps'-law fit of the new-gene curve)"
tags: [validation, comparative-genomics, pan-genome]
doc_path: docs/Evidence/PANGEN-HEAP-001-Evidence.md
sources:
  - docs/Evidence/PANGEN-HEAP-001-Evidence.md
source_commit: 7ed06ad09b44f5c210d8477a3ffa537c7e9920f4
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PANGEN-HEAP-001

The validation-evidence artifact for test unit **PANGEN-HEAP-001** — the **pan-genome growth model**:
count how many **new** gene clusters each added genome contributes from a **presence/absence matrix**,
fit the **power law** `y = K·x^(-alpha)` to that new-gene curve by least squares, and classify the
pan-genome **open vs closed** by the fitted decay exponent **alpha**. This is a **pan-genome family**
(`PANGEN-*`) Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm, its formulas, invariants,
worked oracles, and assumptions are summarized in the dedicated concept
[[pan-genome-heaps-law-fit]]. It is the fitting engine underneath the open/closed facet of the
occupancy partition [[pan-genome-core-accessory-partition]] (PANGEN-CORE-001) and consumes the gene
families from the clustering step [[pan-genome-gene-clustering]] (PANGEN-CLUSTER-001). See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **micropan `heaps()`** (`R/powerlaw.R` + CRAN refman, authority 3, implements Tettelin 2008) — the
    reference implementation. Model `y.hat <- p[1]·x^(-p[2])` (p[1] = Intercept/K, p[2] = decay
    exponent alpha); **binarization** (counts > 0 → 1); the **first-appearance** new-cluster rule
    (`(cm==1)[i] & (cm==0)[i-1]` — a cluster is new at genome i when its cumulative presence is 1 at
    row i and 0 at row i−1); **genome index starts at 2** (first genome has no predecessor); pooling
    over `n.perm` random orderings (default 100, "certainly a minimum"); least-squares objective
    `J = sqrt(Σ(y−K·x^(-alpha))²)/|x|` minimized by `optim(method="L-BFGS-B")` over **K ∈ [0,10000],
    alpha ∈ [0,2]** from start `(mean(y at x==2), 1)`; verbatim open/closed rule **"if alpha<1.0 the
    pan-genome is open, if alpha>1.0 it is closed."**
  - **Tettelin, Riley, Cattuto, Medini 2008** (Curr Opin Microbiol 11(5):472–477, authority 1) — the
    new-gene-discovery curve modelled by a **power law** (not only exponential decay); openness judged
    by the power-law decay exponent; pan-genome = core (all strains) + dispensable (subset-shared +
    strain-specific); **open** = keeps gaining new genes, **closed** = bounded gene pool.
  - **Tettelin et al. 2005** (PNAS 102(39):13950–13955, authority 1, coined "pan-genome") — qualitative
    anchor for *S. agalactiae*: core-decay `F_c = κc·exp[−n/τc] + Ω` (Ω = 1806 extrapolated core,
    r²=0.990) and new-gene `F_s = κs·exp[−n/τs] + tg(θ)` (tg(θ) = **33 ± 3.5** asymptotic new genes per
    genome, r²=0.995); 2nd genome added **161** new genes, 5th added **54**, core ≈ **80%** of a single
    genome; asymptote nonzero (p < 6×10⁻⁴) → **open** pan-genome.
- **Corner cases / failure modes:** < 2 genomes → new-gene curve `x = 2:ng` empty, no fit defined
  (contract returns degenerate Intercept = 0, not an exception); genome index starts at 2 (first genome
  no predecessor); **binary presence only** (copy-number > 1 / duplicate ids collapse to 1 before
  counting); alpha clamped to [0,2] with alpha = 1 the open/closed boundary (strict `<` → boundary is
  not-open/closed).
- **Datasets (documented oracles):**
  - *Closed-form power curve (exact, derived)* — x = [2,3], y = [8,4]: solving `8/4 = (3/2)^alpha` →
    **alpha = ln 2/ln(3/2) ≈ 1.70951**, **K = 8·2^alpha ≈ 26.16400**, J = 0; alpha > 1 → **closed**.
  - *Constant new-gene curve (open boundary, exact)* — x = [2,3], y = [1,1]: best power fit `K·x^(-alpha)
    ≡ 1` → **alpha = 0, K = 1**, J = 0; alpha < 1 → **open**.
  - *Tettelin 2005 qualitative anchor* — 161 / 54 new genes at genomes 2 / 5, tg(θ) = 33 ± 3.5, ≈ 80%
    core, openness = open.
- **Coverage recommendations:** MUST-test the first-appearance new-gene-count extraction, exact
  (K, alpha) recovery on the closed curve (1.7095 / 26.1640 / not-open), the constant curve (alpha 0 /
  K 1 / open), the strict alpha < 1 / alpha > 1 classification, the < 2-genome degenerate fit (no
  exception), and binarization (a cluster present multiple times counts once); SHOULD-test that the
  dictionary/convenience overload clusters-then-fits and agrees with the matrix path; COULD-test that
  the predictor is monotone decreasing in N for alpha > 0.

## Deviations and assumptions

Two assumptions, both source-backed, no deviations from the literature. (1) **Optimizer method
(assumption).** micropan calls `optim(method="L-BFGS-B")`; this implementation minimizes the *identical*
objective over the *identical* box constraints (K ∈ [0,10000], alpha ∈ [0,2]) from the *identical* start
point using deterministic bounded coordinate descent. The method is non-correctness-affecting — the
objective, bounds, and start (which determine the result) are copied verbatim, and for exact
power-curve data within bounds the global minimum is unique and both optimizers reach it (recovered
(K, alpha) matches the analytic solution to < 1e-9). (2) **Permutation RNG (assumption).** micropan uses
R's `sample()` (permutation-dependent pooled curve); this implementation fixes the seed and uses natural
input order for the first permutation, so single-permutation fixed-order fits are exactly reproducible.
The pool-over-orderings averaging principle matches the source; the specific RNG stream does not affect
correctness for the deterministic test matrices. No source contradictions.
