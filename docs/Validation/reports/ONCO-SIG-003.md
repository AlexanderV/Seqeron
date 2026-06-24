# Validation Report: ONCO-SIG-003 — Signature Exposure Estimation: Bootstrap Confidence Intervals

- **Validated:** 2026-06-24   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates, confidence, seed, resampling)`; enum `BootstrapResampling` (Multinomial / Poisson); internals `MultinomialResample`, `PoissonResample`, `SamplePoisson` (Knuth), `SampleBinomial`, `Percentile` (type-7), `Mean`; point estimate via `FitSignatures` (NNLS / Lawson-Hanson, ONCO-SIG-002).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

> Re-validation trigger: commit `4b5c1160` added the Senkin (2021) MSA **Poisson** resampling
> variant behind a `BootstrapResampling` enum that defaults to `Multinomial`. Prior cross-checklist
> verification was reset (☑→☐). This session re-validates both schemes from external sources.

## Stage A — Description

### Sources opened & what they confirm
- **Senkin S. (2021), MSA, BMC Bioinformatics 22:540** — re-fetched PMC8567580 this session.
  Confirmed verbatim:
  - *"For each bootstrap sample, NNLS attribution is applied to derive the vector of signature activities."*
  - *"95% confidence intervals are then derived for each signature attribution by taking [2.5%, 97.5%] percentiles of the resulting bootstrap activities."*
  - Poisson variant: *"Since the mutational burden of each bootstrap sample is fixed and equal to M, we slightly modify the method by drawing counts from independent binomial distributions, so that the total mutational burden is no longer fixed"*, and *"for any given mutation category … the distribution of bootstrapped mutation counts follows a Poisson distribution."*
  These match the Evidence artifact's quotes exactly. The Poisson construction = per-channel
  Poisson(observedₖ), N not fixed, NNLS refit per replicate, [2.5%, 97.5%] percentile CI.
- **sigminer `sig_fit_bootstrap`** (per Evidence, re-confirmed 2026-06-23 in the commit) — fixed-N
  **multinomial** resample `sample(K, ΣN, replace=TRUE, prob=catalog/Σcatalog)`; this is the default
  scheme. Cites Huang 2018 as the underlying method.
- **Efron (1979) percentile method** — 95% CI = [2.5%, 97.5%] percentiles; generally
  [½(1−c), 1−½(1−c)].
- **Hyndman & Fan (1996) type-7 quantile** — R/NumPy default linear interpolation; 0-based rank
  h = p·(n−1), Q(p) = x₍⌊h⌋₎ + (h−⌊h⌋)(x₍⌈h⌉₎ − x₍⌊h⌋₎).
- **Knuth, TAOCP Vol. 2 §3.4.1** (Poisson deviate) — re-confirmed via web search: L = e^(−λ),
  k = 0, p = 1, multiply uniforms until p ≤ L, return k−1. This is exactly the production
  `SamplePoisson`.

### Formula check
- Point estimate = NNLS exposures of the **observed** (un-resampled) catalog — matches Senkin/Huang.
- Multinomial resample: Multinomial(N, p), N = Σ catalog, pₖ = catalogₖ/N — matches sigminer/Huang.
- Poisson resample: each channel k ~ Poisson(observedₖ), independent, N free — matches Senkin Poisson variant.
- Percentile interval: Lower = Q((1−c)/2), Upper = Q(1−(1−c)/2) → 2.5/97.5 at c=0.95 — matches Efron/Senkin.
- Quantile: 0-based h = p·(n−1) linear interpolation = type-7 — matches Hyndman & Fan.

### Edge-case semantics check
- **N = 0** → every multinomial resample all-zero → all replicate exposures 0 → [0,0]. ✅
- **Single non-zero channel (multinomial)** → deterministic collapse → every replicate = point estimate. ✅
- **R = 1** → percentile of one element = that element → lower = upper = mean. ✅
- **Poisson zero-count channel** → Poisson(0) = 0 every replicate → that signature ≡ 0. ✅
- **Poisson single non-zero channel** → variance = mean > 0 → positive interval width (NOT collapse). ✅

### Independent cross-check (numbers, recomputed this session)
- **NNLS planted-exposure recovery (scipy `nnls`):** 5-channel × 3-signature toy matrix, planted
  exposures [40, 25, 10]; built catalogue = W·planted; refit → recovered **[40, 25, 10]** exactly.
  Confirms the refitting recovers known exposures.
- **Type-7 percentiles vs `numpy.quantile(method='linear')`:** [0,1,2,3,4] at p=0.025/0.5/0.975 →
  **0.1 / 2.0 / 3.9**; [2,4,6,8] → **2.15 / 5.0 / 7.85**. Matches the Evidence table and the code's `Percentile`.
- **Poisson bootstrap CI brackets the point estimate:** 20 000 Poisson(12) draws → 95% type-7 CI
  **[6, 19]**, brackets the point estimate 12. ✅

### Findings / divergences
- **Resolved divergence (prior report's PASS-WITH-NOTES note):** the prior report flagged that only
  the sigminer fixed-N multinomial was implemented while Senkin uses an unconditioned Poisson resample.
  Commit `4b5c1160` adds the Poisson scheme behind an opt-in enum (default still multinomial). Both
  published schemes are now available and externally grounded; the divergence is resolved → **PASS**.

## Stage B — Implementation

### Code path reviewed
- `BootstrapExposures` — `OncologyAnalyzer.cs:4774–4872` (validation, point estimate, replicate loop, percentile CI)
- `MultinomialResample` (sequential conditional-binomial) — `:4882–4918`
- `PoissonResample` — `:4929–4936`; `SamplePoisson` (Knuth) — `:4944–4961`
- `SampleBinomial` (sum of Bernoulli) — `:4968–4990`; `Mean` — `:4993–5002`; `Percentile` (type-7) — `:5012–5034`
- `BootstrapResampling` enum — `:4681`

### Formula realised correctly? (evidence)
- **Point estimate** = `FitSignatures(observed, signatures).Exposures` (un-resampled). ✅ (Senkin/Huang; INV-5)
- **Multinomial** via conditional binomials: channel k draws Binomial(remaining, pₖ/Σ_{i≥k}pᵢ); total
  conserved; N=0 → all zeros; last channel with mass takes the remainder. Standard exact construction. ✅
- **Poisson**: independent `SamplePoisson(observedₖ)` per channel; N not conserved; λ≤0 → 0. ✅
- **SamplePoisson** byte-matches Knuth (L=e^(−λ), multiply uniforms, return #draws−1). ✅
- **Percentile**: 0-based h = p·(n−1), linear interpolation; n=1 returns the element. Type-7. ✅
- **CI probabilities** `(1−c)/2` and `1−(1−c)/2`. ✅

### Cross-verification table recomputed vs code
| Check | External value | Code result |
|-------|----------------|-------------|
| NNLS recovery [40,25,10] | [40,25,10] (scipy nnls) | M4 / FitSignatures path (exact) |
| type-7 [0..4] @0.025/0.5/0.975 | 0.1 / 2.0 / 3.9 (numpy) | `Percentile` test exact (1e-12) |
| type-7 [2,4,6,8] @0.025/0.5/0.975 | 2.15 / 5.0 / 7.85 (numpy) | `Percentile` test exact (1e-12) |
| Poisson(12) 95% CI brackets 12 | [6,19] ⊃ 12 | P1/P3 positive-width, brackets point |
| Poisson(0)=0 channel | 0 always | P2 exact 0 |

### Variant/delegate consistency
- Default (no `resampling`) is byte-for-byte equal to explicit `Multinomial` (P5, INV-7). ✅
- Poisson path produces a genuinely different bootstrap distribution (P5). ✅
- No `*Fast`/delegate variants. Defaults: replicates 1000, confidence 0.95, seed 42 — source-attributed.

### Test quality audit
- **Sourced, exact assertions** (`Within(1e-10/1e-12)`) for all deterministic cases; ordering/non-negativity
  invariants on non-degenerate catalogs; determinism (M6, P4); failure modes (8 guards incl. undefined-scheme).
- **P1 is a real distribution lock, not a tautology:** the test re-derives the per-replicate exposures with
  an *independent* Knuth Poisson reference consuming `Random` in the same seeded order, and asserts the
  bootstrap mean and **both** type-7 percentile bounds match exactly. It would fail against any non-Poisson
  resample or any RNG-order mismatch.
- **Type-7 interpolation locked on non-constant samples** (M8b) via reflection on the two NumPy-confirmed
  datasets — fails against nearest-rank / 1-based / non-type-7 rules. (This was the prior session's A46 fix.)
- **Honest green:** ONCO-SIG-003 class = **25 passed, 0 failed**. Build: 0 warnings, 0 errors. No code
  changed this session; full suite not required, but the class is green against a freshly built test DLL.

### Findings / defects
- None. No implementation, description, or test-quality defect found this session.

## Verdict & follow-ups
- **Stage A: PASS** — multinomial and Poisson schemes, NNLS point estimate, percentile method, and type-7
  quantile are all externally sourced and correct; the prior multinomial-vs-Poisson note is resolved by the
  added opt-in Poisson variant.
- **Stage B: PASS** — code faithfully realises both schemes (Poisson = Knuth, multinomial = conditional
  binomial, type-7 = NumPy-confirmed); 25/25 tests green; no defect.
- **End-state: ✅ CLEAN** — algorithm fully functional; no code changed.
