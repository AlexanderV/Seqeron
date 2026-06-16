# Validation Report: ONCO-SIG-002 — Mutational Signature Fitting (NNLS Refitting + Cosine Similarity)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CosineSimilarity(a, b)`, `OncologyAnalyzer.FitSignatures(catalog, signatures)`, `OncologyAnalyzer.ReconstructCatalog(signatures, exposures)`, `SignatureFitResult` (record struct)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

| Source | URL | What it confirmed |
|--------|-----|-------------------|
| Non-negative least squares (Wikipedia, cites Lawson & Hanson 1974, Ch. 23) | https://en.wikipedia.org/wiki/Non-negative_least_squares | NNLS objective `argmin_x ‖Ax − y‖₂²  s.t. x ≥ 0`; Lawson-Hanson active-set algorithm with gradient `w = Aᵀ(y − Ax)`, main loop `while R≠∅ and max(w_R) > ε`, inner-loop step `α = min(x_i/(x_i − s_i))` over negative components, removing non-positive variables back to R; convex problem with convex feasible set ⇒ an interior non-negative solution is the unique global optimum (NNLS = unconstrained LS when already feasible). |
| Cosine similarity (Wikipedia) | https://en.wikipedia.org/wiki/Cosine_similarity | `S_C(A,B) = (A·B)/(‖A‖‖B‖) = ΣAᵢBᵢ / (√ΣAᵢ²·√ΣBᵢ²)`; range [0,1] for non-negative vectors (0 = orthogonal, 1 = identical direction); divides by the product of norms ⇒ undefined for a zero-norm vector. |
| Blokzijl et al. 2018, MutationalPatterns, Genome Medicine 10:33 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5922316/ | Verbatim cosine formula `sim_AB = ΣAᵢBᵢ/(√ΣAᵢ²·√ΣBᵢ²)`, "value between 0 and 1; identical when 1, independent when 0". Verbatim NNLS `min_x ‖S·x − d‖²₂, x ≥ 0` (S = signature matrix, x = weights, d = 96-channel count vector). Cosine-similarity threshold = 0.95 marks successful reconstruction. |
| Rosenthal et al. 2016, deconstructSigs, Genome Biology 17:31 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4762164/ | Reconstruction `R = S·W`; minimises SSE between tumour profile T and reconstruction; "any coefficient must be greater than 0, as negative contributions make no biological sense"; weights "normalized between 0 and 1". |

### Formula check

- **Cosine similarity** (doc §2.2, `OncologyAnalyzer.cs:2322`): matches Blokzijl 2018 and Wikipedia verbatim — dot product over the product of Euclidean norms; range [0,1] for non-negative inputs. ✅
- **NNLS objective** (doc §2.2, §4): `minₓ ‖S·x − d‖₂², x ≥ 0` matches Blokzijl 2018 and the NNLS primary statement verbatim. ✅
- **Lawson-Hanson active set** (doc §4.2): gradient `w = Sᵀ(d − Sx)`, max-gradient index selection, normal-equations LS on the passive set, bounded inner step `α = min x_i/(x_i − s_i)` — matches the cited algorithm step-for-step. ✅
- **Reconstruction** `S·x` and **proportion normalisation** (÷ Σ exposures) match deconstructSigs `R = S·W` and "weights normalized between 0 and 1". ✅

### Edge-case semantics check

- Zero catalog d = 0 ⇒ only feasible minimiser of ‖Sx‖², x≥0 is x = 0 (sourced; convexity + non-negativity). ✅
- Zero-norm cosine ⇒ genuinely undefined per the formula (÷0); the doc states 0.0 as a *documented assumption* (Assumption 2), not a sourced value. Acceptable and flagged. 🟡 (assumption, declared)
- Identical ⇒ 1, orthogonal ⇒ 0 (sourced). ✅
- Negative unconstrained coefficient ⇒ clamp to 0 + refit on the remaining active set (sourced from the active-set algorithm). ✅

### Independent cross-check (numbers)

Recomputed every expected value independently with Python (`scipy.optimize.nnls`, NumPy) and by hand from the formulas — NOT from the repo code:

| Case | Independent result | Spec/test expected |
|------|--------------------|--------------------|
| cos([1,2,3],[1,2,3]) | 1.0 | 1.0 (M1) ✅ |
| cos([1,0],[0,1]) | 0.0 | 0.0 (M2) ✅ |
| cos([1,1],[1,0]) | 0.7071067811865475 | 1/√2 (M3) ✅ |
| cos([3,4],[6,8]) | 1.0 | 1.0 (M4) ✅ |
| NNLS S=I, d=[3,5] | [3,5] | [3,5] (M5) ✅ |
| NNLS sig1=[1,0],sig2=[1,1], d=[0,1] | [0, 0.5] | [0,0.5] (M6) ✅ |
| proportions of [3,5] | [0.375, 0.625] | [0.375,0.625] (M9) ✅ |
| INV-5: d=[4,1,7], S=[[1,0],[1,1],[0,1]] | fitSSE 33.33 ≤ ‖d‖²=66 | ≤ (S3) ✅ |
| recon S·x, sig=[1,0]/[1,1], x=[2,3] | [5,3] | [5,3] ✅ |
| imperfect fit sig=[1,1,1], d=[3,0,0] | exposure 1, recon [1,1,1], recon cos = 1/√3 = 0.5773502691896258 | (new test) ✅ |

### Findings / divergences

None material. The only divergences from a pure "sourced" value are the two declared assumptions (zero-norm cosine = 0.0; proportion normalisation = ÷Σ), both documented and both standard. Stage A: **PASS**.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `CosineSimilarity` (2335–2368): null guards, empty/length-mismatch `ArgumentException`, single-pass dot + both norms, zero-norm ⇒ 0.0.
- `ReconstructCatalog` (2383–2409): validates signatures, exposure-count match; computes Σⱼ sig[j][k]·x[j].
- `FitSignatures` (2431–2452): validates, solves NNLS, reconstructs, normalises, computes reconstruction cosine; returns `SignatureFitResult`.
- `SolveNonNegativeLeastSquares` (2489–2590): faithful Lawson-Hanson active set (gradient select, normal-equations LS on passive set via Gaussian elimination with partial pivoting, bounded inner step, KKT termination at ε=1e-12, finite-iteration safety cap).
- Helpers `ComputeResidual`, `ComputeGradient`, `SolveLeastSquaresOnPassiveSet`, `SolveLinearSystem`, `NormalizeExposures`, `ValidateSignatures`.

### Formula realised correctly?

Yes. The cosine code is the formula verbatim. The NNLS solver realises the Lawson-Hanson active set exactly as the cited algorithm describes (gradient `Sᵀ(d−Sx)`, passive/active sets, bounded step, KKT stop). Reconstruction = S·x; normalisation = ÷Σ exposures (all-zero when Σ=0). Reconstruction cosine reuses `CosineSimilarity(d, S·x)`. All independently cross-checked values above match the actual code (the full suite passes against them).

### Cross-verification table recomputed vs code

Every value in the Stage-A table was verified to be produced by the actual code via the passing test suite (M1–M10, S1–S3, C1, guards, plus the three new tests). The new imperfect-fit test confirms the reconstruction-cosine path on a strictly sub-unity, sourced value (1/√3).

### Variant/delegate consistency

`FitSignatures` reuses `ReconstructCatalog` and `CosineSimilarity` (no divergent re-implementation). `ReconstructCatalog` is exercised both directly and via `FitSignatures`. ONCO-SIG-003 (`BootstrapExposures`) reuses this `FitSignatures`/NNLS, consistent.

### Numerical robustness

Single-pass accumulation; squared norms compared to exact 0.0 before dividing (no ÷0); Gaussian elimination with partial pivoting; singular passive-set matrix (collinear signatures) leaves that component at 0 rather than throwing (documented). ε=1e-12 and a finite outer-iteration cap guard termination. No overflow concern on count-scale inputs.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** every numeric expectation traces to a formula/source and was reproduced this session via scipy/NumPy/hand — not read off the implementation. A deliberately-wrong NNLS (e.g. returning the unconstrained [−1,1]) would fail M6/M10; a wrong cosine would fail M1–M4. Not green-washed.
- **No weakening:** exact `Within(1e-9..1e-12)` on canonical values. The two inequality assertions (`GreaterThanOrEqualTo(0.0)` in M10, `LessThanOrEqualTo` in S3) are *invariant* checks (INV-4, INV-5) that accompany exact value assertions — legitimate, not range-widening.
- **Coverage gaps found and filled (this session):** the existing 14 tests covered M1–M10, S1–S3, C1, null/length/dim/no-signature guards, and direct reconstruction. Missing were: (a) the reconstruction-cosine path on an imperfect/under-determined fit (all prior fits reconstructed exactly, cosine = 1); (b) the empty-vector `ArgumentException` branch of `CosineSimilarity`; (c) null-argument guards on `ReconstructCatalog`. Added three tests with exact sourced values (imperfect fit recon cos = 1/√3; empty-vector throw; null `ReconstructCatalog`).
- **Honest green:** FULL unfiltered suite = **6640 passed, 0 failed, 0 skipped** (was 6637 + 3 new); `dotnet build` 0 errors, no new warnings in the changed test file.

**Gate result: PASS.**

### Findings / defects

No code defect. Pre-existing tests were sound; three coverage gaps were closed with exact, externally-sourced expectations. No code change was required.

## Verdict & follow-ups

- **Stage A: PASS.** Formulas (cosine similarity, NNLS objective, Lawson-Hanson active set, reconstruction, proportion normalisation) match Blokzijl 2018, Rosenthal 2016, and the Lawson-Hanson primary algorithm verbatim. Two declared, standard assumptions (zero-norm cosine = 0.0; ÷Σ proportions).
- **Stage B: PASS.** Code faithfully realises the validated description; all values independently cross-checked against scipy.optimize.nnls and hand computation; tests strengthened to cover the reconstruction-quality (<1), empty-vector, and null-reconstruction paths.
- **End-state: ✅ CLEAN.** No defect; algorithm fully functional. Test-quality gate PASS. Full suite green (6640/0).
