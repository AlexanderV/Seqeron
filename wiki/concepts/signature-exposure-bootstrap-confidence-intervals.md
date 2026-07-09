---
type: concept
title: "Signature exposure bootstrap confidence intervals (resample-refit-percentile)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-SIG-003-Evidence.md
source_commit: 2c404cc9eb23b82d7b7e6aeb757a40051a8f84fd
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-sig-003-evidence
      evidence: "Test Unit ID: ONCO-SIG-003 ... Algorithm: Signature Exposure Estimation — Bootstrap Confidence Intervals"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:mutational-signature-fitting-and-extraction
      source: onco-sig-003-evidence
      evidence: "For each bootstrap sample, NNLS attribution is applied to derive the vector of signature activities — the same NNLS refit (FitSignatures) of ONCO-SIG-002 is re-run per replicate."
      confidence: high
      status: current
---

# Signature exposure bootstrap confidence intervals

The Oncology family's **signature-exposure uncertainty** unit (**ONCO-SIG-003**): a resample-and-refit loop
that turns the single point-estimate exposure vector produced by the NNLS refit
([[mutational-signature-fitting-and-extraction]], ONCO-SIG-002) into a **per-signature confidence interval**.
It adds no new decomposition — it quantifies how stable the existing fit is under resampling of the mutation
catalog. The literature-traced record is [[onco-sig-003-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## The pipeline: resample → refit → percentile

Three steps, each pinned to a primary source:

1. **Resample** the observed 96-channel catalog `R` times (default 1000; sigminer recommends ≥ 100).
2. **Refit** each resampled catalog with the same **NNLS** attribution used by ONCO-SIG-002
   (`FitSignatures`), producing one exposure vector per replicate — an `R × signatures` bootstrap matrix.
3. **Percentile CI:** for each signature, the confidence interval is the **[2.5%, 97.5%] percentiles** of its
   `R` bootstrap exposures (Efron 1979 percentile method; for confidence `c`, the
   `[100·½(1−c), 100·(1−½(1−c))]` percentiles). One interval per signature, in signature order; the reported
   **point estimate** is `FitSignatures(observed).Exposures` (contract consistency with ONCO-SIG-002).

Percentiles are read off the finite bootstrap sample by the **type-7** sample-quantile (R/NumPy default;
Hyndman & Fan 1996): `h = p·(n−1)`, `Q(p) = x₍⌊h⌋₎ + (h−⌊h⌋)·(x₍⌊h⌋₊₁₎ − x₍⌊h⌋₎)`.

**Type-7 worked oracles:** sorted `[0,1,2,3,4]` (n=5) → p=0.5 → 2.0, p=0.025 → 0.1, p=0.975 → 3.9; sorted
`[2,4,6,8]` (n=4) → p=0.025 → 2.15, p=0.5 → 5.0, p=0.975 → 7.85.

## Two resampling schemes: multinomial vs Poisson

The only place the two source lineages differ is **how a replicate catalog is drawn** — and specifically
**whether the total mutational burden N is held fixed**:

- **Multinomial (default; sigminer, fixed N).** Draw `N = Σ catalog` category indices with replacement,
  weighted by `pₖ = catalogₖ / N`, and tabulate → `multinomial(N, p)`. The total burden is conserved every
  replicate. This is the byte-for-byte default; calling without a `resampling` argument equals explicit
  `Multinomial`.
- **Poisson (Senkin 2021 variant, unfixed N).** Draw each channel independently as **Poisson(observedₖ)**;
  N is **not** fixed. Senkin's justification: the conditional distribution of a vector of independent Poisson
  variables **equals** the multinomial — so the two schemes are the same method with/without conditioning on N.
  The observable consequence is **variance = mean** per channel.

Both feed the identical refit + percentile steps; the `Poisson` path produces a genuinely different bootstrap
distribution from `Multinomial`.

## Corner cases (the discriminating tests)

| Case | Multinomial | Poisson |
|------|-------------|---------|
| Zero-mutation catalog (N = 0) | every replicate 0 → interval `[0,0]` | same — Poisson(0) = 0 |
| Single non-zero channel `[10]`/`[12]` over `[[1.0]]` | **deterministic collapse**: draw always `[N]` → exposure = point estimate every replicate → point/mean/lower/upper equal | **NOT degenerate**: Poisson(λ>0) fluctuates (var = mean) → **positive interval width** |
| Zero-count channel among others | resamples away as usual | exactly 0 every replicate (Poisson(0)=0) |
| Single replicate (R = 1) | lower = upper = mean = that one value | same |

The single-channel row is the **defining distinction**: the multinomial draw over one outcome is deterministic
(so its interval collapses), whereas Poisson(λ) has strictly positive width. This is exactly why the Poisson
variant was added (limitation fix) — the fixed-N multinomial understates uncertainty when the catalog is
concentrated.

## Invariants and validation handles

- **Ordering:** `Lower ≤ point-estimate region ≤ Upper` and `Lower ≤ Mean ≤ Upper`; the 97.5 bound ≥ the 2.5
  bound (Efron percentile ordering).
- **Non-negativity:** every bootstrap exposure ≥ 0 (inherited from NNLS, `x ≥ 0`; Lawson & Hanson 1974).
- **Determinism:** two calls with the same seed return identical intervals (fixed default seed 42).
- **Input validation** (mirrored from `FitSignatures`): null catalog/signatures → `ArgumentNullException`;
  ragged/empty signatures, catalog-length mismatch, negative count → `ArgumentException`; replicates < 1,
  confidence ∉ (0,1), or an undefined `resampling` enum → `ArgumentOutOfRangeException`.

## Relation to the oncology family

This is the **uncertainty layer directly above** [[mutational-signature-fitting-and-extraction]]
(ONCO-SIG-002): it `depends_on` that unit's NNLS refit (`FitSignatures`) and re-runs it per replicate, and
sits two steps above the [[sbs96-mutational-signature-catalog]] (ONCO-SIG-001) that builds the 96-channel
spectrum being resampled. A confidence interval on a signature exposure is what makes the fitted
somatic-mutation-process biomarker actionable — a signature whose interval excludes 0 is a *present* process,
distinguishing it from the copy-number-scar [[homologous-recombination-deficiency-score]] and the
mismatch-repair [[microsatellite-instability-detection]], and feeding the clinical-interpretation units
([[cancer-variant-tier-classification-amp-asco-cap]], [[clinical-actionability-oncokb-levels]]) with a
confidence qualifier rather than a bare point estimate.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for the signature-exposure bootstrap. The models
are MSA / Senkin 2021 (parametric bootstrap + Poisson variant), sigminer `sig_fit_bootstrap` (fixed-N
multinomial resample, the default), Huang, Wojtowicz & Przytycka 2018 (bootstrap confidence for decomposition),
Efron 1979 (percentile CI), and Hyndman & Fan 1996 (type-7 quantile). The **type-7** interpolation and a
**fixed seed (42)** are the two source-aligned assumptions; the **multinomial** scheme is the backward-
compatible default and the **Poisson** scheme its unfixed-N alternative. **Not for clinical or diagnostic
use.** No source contradictions.
