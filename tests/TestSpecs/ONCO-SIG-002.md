# Test Specification: ONCO-SIG-002

**Test Unit ID:** ONCO-SIG-002
**Area:** Oncology
**Algorithm:** Mutational Signature Fitting / Refitting (NNLS decomposition + cosine similarity)
**Status:** ☑ Validated (2026-06-16; +3 coverage tests: imperfect-fit recon cosine = 1/√3, empty-vector cosine throw, null ReconstructCatalog)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Blokzijl et al. (2018), MutationalPatterns, Genome Medicine 10:33 | 3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5922316/ | 2026-06-14 |
| 2 | Rosenthal et al. (2016), deconstructSigs, Genome Biology 17:31 | 3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4762164/ | 2026-06-14 |
| 3 | Lawson & Hanson (1974), Solving Least Squares Problems, Ch. 23 (active-set NNLS) | 4 | https://en.wikipedia.org/wiki/Non-negative_least_squares | 2026-06-14 |
| 4 | Pan & Wang (2020), iMutSig | 3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7702159/ | 2026-06-14 |

### 1.2 Key Evidence Points

1. Cosine similarity sim = ΣAᵢBᵢ / (√ΣAᵢ² · √ΣBᵢ²), value in [0,1] — Source 1 (§ Mutational profile similarity); Source 4.
2. Signature refitting solves min_x ‖S·x − d‖₂², x ≥ 0 (NNLS) — Source 1 (§ Finding the contribution of known signatures).
3. Reconstructed profile = S·x (= S·W); residual T − S·W minimised; coefficients > 0; weights normalised between 0 and 1 — Source 2.
4. NNLS solved by the Lawson-Hanson active-set algorithm; when the unconstrained LS solution is already non-negative, NNLS returns it unchanged — Source 3.
5. Reconstruction cosine similarity is the quality check; ≥ 0.95 = successful reconstruction — Source 1.

### 1.3 Documented Corner Cases

- Zero observed catalog (d = 0): NNLS minimiser is x = 0; reconstruction = 0 (Source 1/3).
- Cosine similarity undefined for a zero-norm vector (division by 0) — degenerate case (Source 1/4 formula).
- Identical vectors ⇒ cosine = 1; orthogonal (disjoint support) ⇒ cosine = 0 (Source 1/4).
- Unconstrained LS coefficient negative ⇒ clamped to 0, refit on the remaining active set (Source 3).

### 1.4 Known Failure Modes / Pitfalls

1. Dividing by zero norm in cosine similarity for an empty / zero catalog — must be handled — Source 1/4 formula.
2. Returning negative exposures (biologically meaningless) — NNLS constraint x ≥ 0 forbids — Source 2/3.
3. Dimension mismatch between catalog length and signature row count — invalid input — derived from S·x shape.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CosineSimilarity(a, b)` | OncologyAnalyzer | Canonical | Source 1/4 formula |
| `FitSignatures(catalog, signatures)` | OncologyAnalyzer | Canonical | NNLS min ‖Sx−d‖², x≥0 (Source 1/3) |
| `ReconstructCatalog(signatures, exposures)` | OncologyAnalyzer | Canonical | S·x (Source 2) |
| `SignatureFitResult` (record) | OncologyAnalyzer | Internal | exposures + proportions + reconstruction cosine |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | CosineSimilarity ∈ [0,1] for non-negative inputs | Yes | Source 1/4 (range 0–1) |
| INV-2 | CosineSimilarity(a,a) = 1 for any non-zero a | Yes | Source 1/4 (identical ⇒ 1) |
| INV-3 | CosineSimilarity is scale-invariant: cos(a, k·b) = cos(a, b), k>0 | Yes | Source 1/4 (cosine of angle) |
| INV-4 | All fitted exposures ≥ 0 | Yes | Source 2/3 (x ≥ 0) |
| INV-5 | Fit residual SSE ‖Sx−d‖² ≤ SSE of x=0 (= ‖d‖²) | Yes | Source 1/3 (x=0 is feasible; NNLS is the minimiser) |
| INV-6 | Normalised exposures sum to 1 when Σexposures > 0 (else all 0) | Yes | Source 2 (weights normalised 0–1 as proportions) |
| INV-7 | When the unconstrained LS solution is non-negative, NNLS = that solution | Yes | Source 3 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Cosine identical | a=b=[1,2,3] | 1.0 | Source 1/4: 14/(√14·√14)=1 |
| M2 | Cosine orthogonal | [1,0] vs [0,1] | 0.0 | Source 1/4: 0/(1·1)=0 |
| M3 | Cosine general | [1,1] vs [1,0] | 0.70710678118654752 | Source 1/4: 1/(√2·1) |
| M4 | Cosine scale-invariant | [3,4] vs [6,8] | 1.0 | Source 1/4: 50/(5·10)=1 |
| M5 | NNLS identity recovery | S=[[1,0],[0,1]], d=[3,5] | exposures=[3,5] | Source 1/3 (unconstrained LS, non-neg) |
| M6 | NNLS constraint binds | S=[[1,1],[0,1]], d=[0,1] | exposures=[0,0.5] | Source 3 active-set (clamp + refit) |
| M7 | Reconstruction = S·x | S=[[1,0],[0,1]], exposures=[3,5] | recon=[3,5] | Source 2 (R=S·W) |
| M8 | Reconstruction cosine | identity fit of d=[3,5] | reconstructionCosine=1.0 | Source 1 (quality check) |
| M9 | Normalised exposures | exposures [3,5] → proportions | [0.375,0.625], sum=1 | Source 2 (normalise 0–1) |
| M10 | All exposures non-negative | constraint-binding fit M6 | min exposure = 0 ≥ 0 | Source 2/3 (x≥0) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Zero catalog | S=[[1,0],[0,1]], d=[0,0] | exposures=[0,0], recon=[0,0], proportions=[0,0] | Source 1/3 degenerate |
| S2 | Cosine zero vector | [0,0] vs [1,1] | 0.0 (documented) | Assumption 2 |
| S3 | Residual SSE bound | random-ish d, identity S | SSE ≤ ‖d‖² | INV-5 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Property: scale invariance | cos(a, k·a) = 1 for several k | 1.0 | INV-3 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for signature fitting / cosine similarity in the Oncology test folder. ONCO-SIG-001 covers
  only `ClassifySbsContext` / `Build96ContextCatalog` / `EnumerateSbs96Channels` (different methods).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M10 | ❌ Missing | New unit, no prior tests |
| S1–S3 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_FitSignatures_Tests.cs` — all cases.
- **Remove:** none (no pre-existing file).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_FitSignatures_Tests.cs | Canonical, this unit | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | S1 | ❌ Missing | Implemented | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented | ✅ Done |
| 13 | S3 | ❌ Missing | Implemented | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented | ✅ Done |
| 15 | guards | ❌ Missing | null + dimension-mismatch throws | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M10 | ✅ | Exact hand-derived values from sources |
| S1–S3 | ✅ | Degenerate / invariant cases |
| C1 | ✅ | Property test (scale invariance) |
| Guards (null, dim mismatch) | ✅ | ArgumentNullException / ArgumentException |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Exposure normalisation = exposures / Σexposures (proportions summing to 1; all-zero when Σ=0) | M9, S1 |
| 2 | Cosine similarity returns 0.0 when either vector has zero norm | S2 |

---

## 7. Open Questions / Decisions

1. Decision: implement the Lawson-Hanson active-set NNLS (Source 3) as the deterministic solver for the
   MutationalPatterns NNLS objective (Source 1); deconstructSigs' iterative greedy heuristic (Source 2) is NOT
   reproduced (it is a non-deterministic threshold heuristic), but its reconstruction model R=S·W,
   non-negativity, and proportion normalisation are adopted as they coincide with the NNLS minimiser.
