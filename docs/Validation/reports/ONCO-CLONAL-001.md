# Validation Report: ONCO-CLONAL-001 — Clonal vs Subclonal Mutation Classification

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.ClassifyClonality(variants, purity)`, `OncologyAnalyzer.IdentifyClonalMutations(ccfValues)` (+ internal `ClassifyOne`, `BinomialLikelihoodKernel`; records `ClonalityVariant`, `ClonalityCall`, `ClonalityResult`, enum `ClonalityStatus`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm

1. **Landau DA et al. (2013), *Cell* 152(4):714–726** — primary source, retrieved this session via WebFetch of PMC3575604 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3575604/). Confirmed verbatim:
   - Expected allele fraction **f(c) = αc / (2(1−α) + αq)** (α = purity, c = CCF, q = absolute somatic copy number).
   - Posterior **P(c) ∝ Binom(a | N, f(c))**, **uniform prior on c**, over a **regular grid of 100 c values** spanning **c ∈ [0.01, 1]**, normalised by its sum.
   - Classification rule: a mutation is **clonal** if the **posterior probability that CCF > 0.95 is > 0.5**, **subclonal** otherwise.
   - Sanity reference: subclonal sSNV **median CCF = 0.49** (range 0.11–0.89) — confirms subclonal mutations sit well below 1.

2. **Satas G et al. (2021), *Cell Systems* 12(10):1004–1018 (DeCiFering)** — retrieved via WebFetch of PMC8542635. Confirmed:
   - **Eq. 1:** `c ≈ (1/ρ)·(ρ·N_tot + (1−ρ)·2)/M · v̂` (ρ = purity, N_tot = total CN in cancer cells, M = mutation multiplicity, v̂ = VAF).
   - Algebraically inverting Eq. 1 gives `v̂ = f(c) = ρ·M·c / (ρ·N_tot + 2(1−ρ))`, i.e. with N_tot = q this is exactly the multiplicity-general form `f(c) = ρ·M·c / (2(1−ρ) + ρ·q)`. This reduces to Landau's M=1 expression.
   - Definitions: clonal ⇔ CCF ≈ 1; subclonal ⇔ CCF ≪ 1.

### Formula check
The repo description (`docs/algorithms/Oncology/Clonal_Subclonal_Classification.md` §2.2) states `f(c) = α·M·c / (2(1−α) + α·q)`, uniform-prior binomial posterior on a 100-point grid c ∈ [0.01, 1], clonal iff `P(CCF>0.95) > 0.5`. Every element matches the two primary sources verbatim. The M-general factor is correctly attributed to Satas Eq. 1 and reduces to Landau (M=1).

### Edge-case semantics check
- Grid lower bound 0.01 (a detected mutation is in ≥ 1 cancer cell) — sourced to Landau.
- Probabilistic (posterior mass), not point-estimate, classification — sourced to Landau; a near-1 point estimate with shallow coverage can still be subclonal (demonstrated by case M5b: VAF 0.5, N=100 → P(CCF>0.95)=0.443 < 0.5 → Subclonal).
- Multiplicity M raises CCF for the same VAF — sourced to Satas.
- Point-estimate helper threshold CCF > 0.95 **strict** (0.95 excluded) — sourced to Landau ("> 0.95").
- Domain validation (purity ∈ (0,1]; a ∈ [0,N]; q ≥ 1; M ∈ [1,q]; CCF ∈ [0,1]) is standard and consistent with the model (division by purity, binomial counts).

### Independent cross-check (numbers)
A standalone Python reimplementation of the Landau grid posterior (written from the sourced formula, not from the repo code) reproduced the Evidence/TestSpec dataset **exactly**:

| Case | a | N | q | M | ρ | CCF mean | P(CCF>0.95) | Status |
|------|---|---|---|---|----|----------|-------------|--------|
| M1 | 300 | 300 | 2 | 1 | 1.0 | 0.999486 | 1.000000 | Clonal |
| M2 | 400 | 1000 | 2 | 1 | 0.8 | 0.972455 | 0.864167 | Clonal |
| M3 | 240 | 1000 | 2 | 1 | 0.8 | 0.601297 | 0.000000 | Subclonal |
| M4 | 200 | 1000 | 2 | 1 | 1.0 | 0.401198 | 0.000000 | Subclonal |
| M5 | 100 | 100 | 2 | 2 | 1.0 | 0.994330 | 0.998016 | Clonal |
| M5b | 50 | 100 | 2 | 1 | 1.0 | 0.924301 | 0.443319 | Subclonal |

I also confirmed the **clonal/subclonal classification is robust** to the grid convention: re-running with 100 midpoints and with a 10000-point grid leaves the Clonal/Subclonal call unchanged for all six cases (only the descriptive CCF mean and P shift slightly). The exact CCF-mean/P values asserted by the tests are tied to the specific "100 regular endpoint" grid — a faithful reading of Landau's "regular grid of 100 c values" in [0.01,1].

### Findings / divergences
None. Description is mathematically and biologically correct and fully sourced. The registry `ploidy` scalar was superseded by per-variant `LocalCopyNumber` q (documented assumption, API-shape only) — consistent with Landau's per-locus q. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `ClassifyClonality` 4663–4689; `IdentifyClonalMutations` 4702–4725; `ClassifyOne` 4730–4772; `BinomialLikelihoodKernel` 4778–4795; `ValidateClonalityVariant` 4798–4823; constants `ClonalCcfThreshold=0.95` (4573), `ClonalProbabilityThreshold=0.5` (4579), `CcfGridPointCount=100` (4586), `CcfGridLowerBound=0.01` (4592), `CcfGridUpperBound=1.0` (4595), `NormalDiploidCopyNumber=2.0` (68).

### Formula realised correctly?
Yes. `denominator = 2(1−ρ) + ρ·q`; `alleleFractionPerUnitCcf = ρ·M / denominator`; grid `c = 0.01 + step·i`, `step = (1−0.01)/99`, 100 points; `f = min(1, perUnit·c)`; binomial kernel `f^a·(1−f)^(N−a)` in log-space (the constant C(N,a) cancels under normalisation — correct); posterior = weight/sum; `ccfMean = Σ c·posterior`; `probabilityClonal = Σ_{c>0.95} posterior`; clonal iff `probabilityClonal > 0.5`. This is exactly the validated Landau model with the Satas M factor.

### Cross-verification table recomputed vs code
The 18 (now 21) unit tests assert the exact sourced values in the table above and all pass against the actual compiled code. My independent Python computation matches both the tests and the code to 1e-6, so the code reproduces the external reference, not just the test echoes.

### Variant/delegate consistency
`ClonalityVariant` has two constructors (4-arg explicit M, 3-arg default M=1); the 3-arg path is exercised (M1–M4) and the 4-arg path (M5, S4). `IdentifyClonalMutations` uses the same `ClonalCcfThreshold` constant as the posterior classifier — consistent.

### Numerical robustness
Log-space likelihood avoids underflow for large N. The all-zero-posterior fallback (`weightSum ≤ 0`) is a defensive guard that is **unreachable with valid inputs**: `perUnit = ρM/(2(1−ρ)+ρq) ≤ 1` (since M ≤ q and ρ ≤ 1), so f < 1 on interior grid points always yields nonzero likelihood and weightSum > 0. No div-by-zero (purity > 0 enforced; denominator > 0).

### Test quality audit (HARD gate)
- **Sourced, not echoed:** the six modelled cases trace to my independent Python recomputation of the Landau model this session, not to the implementation. A deliberately-wrong implementation would fail them (exact `.Within(1e-6)` on CCF and P).
- **No green-washing:** exact-value assertions with tight tolerances; no Greater/AtLeast/range used where an exact value is known; no skips, no widened tolerances.
- **Coverage:** every public method/overload exercised; both classification branches (clonal/subclonal); multiplicity-load-bearing case (M5b); invariants INV-1..INV-6; empty inputs (C1/C2).
- **Gap found and fixed (Stage-B test-coverage defect, 0 code change):** the **lower-bound** validation branches were untested — `AltReads<0`, `Multiplicity<1`, and `IdentifyClonalMutations` `CCF<0` (only the upper bounds and NaN were asserted). Added one assertion each (`V(-1,100,2)`, `V(50,100,2,0)`, `-0.1`), sourced to the documented domain. Fixture 18 → 21.
- **Honest green:** full unfiltered `dotnet test` = **6661 passed, 0 failed**; `dotnet build` 0 errors (4 pre-existing NUnit2007 warnings in unrelated `ApproximateMatcher_EditDistance_Tests.cs`; the clonality source/test files build warning-free).

### Findings / defects
No code or description defect. One test-coverage gap (untested lower-bound validation) — fixed in-session. **Stage B: PASS.**

## Verdict & follow-ups
- **Stage A: PASS** — description matches Landau (2013) and Satas (2021) verbatim, independently re-confirmed from primary sources; cross-check reproduces all dataset values.
- **Stage B: PASS** — code realises the validated formula exactly; tests assert externally-sourced values; coverage gap (lower-bound validation) closed.
- **End-state: ✅ CLEAN** — fully functional; no remaining defect.
- **Test-quality gate: PASS** (after fix).
