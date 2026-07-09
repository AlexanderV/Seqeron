---
type: source
title: "Evidence: ONCO-SIG-003 (Signature exposure bootstrap confidence intervals)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-SIG-003-Evidence.md
sources:
  - docs/Evidence/ONCO-SIG-003-Evidence.md
source_commit: 2c404cc9eb23b82d7b7e6aeb757a40051a8f84fd
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-SIG-003

The validation-evidence artifact for test unit **ONCO-SIG-003** — **bootstrap confidence intervals for
signature exposures**, the uncertainty-quantification layer that wraps the ONCO-SIG-002 NNLS refit
([[mutational-signature-fitting-and-extraction]]). It is the **thirty-first ingested unit of the Oncology
family** and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern. The distinct method is synthesized in its own concept,
[[signature-exposure-bootstrap-confidence-intervals]]; [[test-unit-registry]] tracks the unit.

This unit answers "**how confident is the fit?**" — not a new decomposition but a resample-and-refit loop that
turns the single point-estimate exposure vector from ONCO-SIG-002 into a per-signature confidence interval.

## What this file records

The Evidence file was built in two dated layers (see its Change History):

- **2026-06-14 — the multinomial bootstrap + percentile CI:**
  - **Senkin 2021 — MSA** (*BMC Bioinformatics* 22:540; rank 1) — parametric bootstrap over the mutation
    catalog; per bootstrap sample, **NNLS attribution is re-applied** to derive the vector of signature
    activities; **95% CIs = the [2.5%, 97.5%] percentiles** of the bootstrap activities; **1000 replicates**
    is the standard count (Figure 2).
  - **sigminer `sig_fit_bootstrap`** (Wang S. et al.; rank 3, citing Huang 2018 as the primary method) — the
    reference implementation's resampling step is `sample(seq(K), total_count, replace = TRUE, prob =
    catalog/sum(catalog))` then `table(...)`, which is **exactly a multinomial(N, p)** draw with
    `N = Σ catalog`, `pₖ = catalogₖ / N` (fixed total burden). Re-verified 2026-06-23: sigminer is
    **multinomial-only** (fixed N). `sig_fit()` refit per replicate; **≥ 100 replicates recommended**.
  - **Huang, Wojtowicz & Przytycka 2018** (*Bioinformatics* 34(2):330–337; rank 1) — the primary method:
    bootstrap to measure **confidence/stability of the decomposition**; implemented by `SignatureEstimation`
    (`bootstrapSigExposures`).
  - **Efron 1979** (*Annals of Statistics* 7(1):1–26; rank 1, via MIT 18.05 notes) — the **percentile
    bootstrap CI**: the 2.5/97.5 percentiles of the bootstrap distribution as the 95% CI; generally
    `[100·½(1−c), 100·(1−½(1−c))]` percentiles for confidence `c`.
  - **Hyndman & Fan 1996 — type-7 quantile** (*The American Statistician* 50(4):361–365; via Wikipedia,
    rank 4) — the R/NumPy default sample-quantile: `h = p·(n−1)`, `Q(p) = x₍⌊h⌋₎ + (h−⌊h⌋)·(x₍⌊h⌋₊₁₎ −
    x₍⌊h⌋₎)`. Realizes the percentile interval on the finite set of bootstrap exposures.

- **2026-06-23 — the Poisson resampling variant (limitation fix):**
  - **Senkin 2021 — Poisson model** (re-retrieved PMC article; rank 1) — MSA's alternative: each mutation
    category count is an **independent Poisson draw**, so the **total burden N is NOT fixed** (contrast the
    fixed-N multinomial). The verbatim justification: *"the conditional distribution of a vector of
    independent Poisson variables is equivalent to multinomial distribution"* — the two schemes differ only
    by whether N is conditioned. Derived construction: channel `k` resampled as **Poisson(observedₖ)**
    independently, NNLS refit per replicate, CIs = the same [2.5%, 97.5%] type-7 percentiles.
  - **Knuth 1997** (*TAOCP* Vol. 2 §3.4.1; rank 5) — the multiplication-of-uniforms Poisson-deviate
    generator used to reproduce the per-replicate Poisson draws in the differential test.

## Corner cases / failure modes

- **Zero-mutation catalog (N = 0):** nothing to resample → every replicate exposure 0 → interval collapses to
  `[0, 0]`.
- **Degenerate single-channel catalog (multinomial):** a multinomial draw over one non-zero channel is
  **deterministic** (all N fall in it) → resample = observed exactly → every replicate exposure = the point
  estimate (`[10]` over `[[1.0]]` → point/mean/lower/upper all 10).
- **Poisson zero-count channel:** `Poisson(0) = 0` every replicate — that channel (and any signature isolated
  to it) is exactly 0 always.
- **Poisson single non-zero channel is NOT degenerate:** `Poisson(λ>0)` fluctuates (**variance = mean = λ**),
  so the single-channel interval has **positive width** — the observable distinction between the Poisson and
  multinomial schemes (`[12]` over `[[1.0]]` → each replicate an independent Poisson(12)).
- **Single replicate (R = 1):** the percentile of a one-element distribution is that value → lower = upper =
  mean = the single replicate's exposure.

## Datasets (deterministic worked oracles)

- **Type-7 percentile worked values** (hand-derived from the Hyndman & Fan formula): sorted `[0,1,2,3,4]`
  (n=5) → p=0.5→2.0, p=0.025→0.1, p=0.975→3.9; sorted `[2,4,6,8]` (n=4) → p=0.025→2.15, p=0.975→7.85,
  p=0.5→5.0.
- **Deterministic single-channel bootstrap (multinomial collapse):** catalog `[10]`, signatures `[[1.0]]` →
  every resampled catalog `[10]`, NNLS exposure 10 each replicate → point/mean/lower/upper = 10.
- **Poisson single-channel bootstrap:** catalog `[12]`, `[[1.0]]` → each replicate exposure an independent
  Poisson(12) draw; mean and [2.5%, 97.5%] type-7 bounds equal those of the same seeded RNG order verified
  against an independent Knuth Poisson generator; a zero-count channel → 0 every replicate.
- **Zero-mutation catalog:** `[0,0,0]` over `[[1,0,0],[0,1,0]]` → all-zero intervals per signature.

## Deviations and assumptions

- **ASSUMPTION — type-7 quantile convention.** Sources fix the percentile *method* ([2.5%, 97.5%]) but not the
  finite-sample interpolation rule. The **type-7** (R/NumPy default; Hyndman & Fan 1996) linear-interpolation
  estimator is adopted — correctness-affecting (it sets the exact bound values) but the documented default of
  the reference tooling (R `quantile`, used by sigminer), so source-aligned rather than invented. The dominant
  degenerate/deterministic tests are convention-independent.
- **ASSUMPTION — fixed RNG seed (42).** The bootstrap is randomized; no source prescribes a seed. A fixed
  default seed makes results reproducible (matching the repository's Phylogenetics bootstrap convention).
  Non-correctness-affecting for the deterministic cases (single-channel collapse, N = 0), whose outcomes are
  seed-independent.
- **Backward-compatibility requirement:** the **multinomial** scheme (sigminer's fixed-N default) remains the
  default — calling without `resampling` must be byte-for-byte equal to explicit `Multinomial`; the `Poisson`
  path produces a genuinely different bootstrap distribution.

No source contradictions — Senkin (MSA), sigminer, Huang 2018, Efron, and Hyndman & Fan are mutually consistent
on the resample-refit-percentile pipeline; the only cross-source nuance is multinomial (fixed N, sigminer) vs
Poisson (unfixed N, Senkin), which the Poisson↔multinomial conditional-equivalence resolves as two variants of
one method.
