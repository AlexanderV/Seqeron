# Test Specification: ONCO-SIG-002

**Test Unit ID:** ONCO-SIG-002
**Area:** Oncology
**Algorithm:** Mutational Signature Fitting / Refitting (NNLS) **+ De-novo Signature Extraction via NMF**
**Status:** ☐ In Progress (2026-06-23: enhancement — automatic NMF rank selection `SelectRank`, KL/Poisson objective overload of `ExtractSignatures`, and cosine reference matching `MatchToReferenceSignatures`; NNLS refitting + Frobenius extraction unchanged)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

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

---

## Appendix A — De-novo NMF Signature Extraction (`ExtractSignatures`, added 2026-06-23)

### A.1 Authoritative Sources (extension)

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 5 | Lee & Seung (2001), Algorithms for NMF, NIPS 13 | 1 | https://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization (proof guide https://arxiv.org/html/2501.11341v1) | 2026-06-23 |
| 6 | Non-negative matrix factorization, Wikipedia (citing Lee & Seung) | 4 | https://en.wikipedia.org/wiki/Non-negative_matrix_factorization | 2026-06-23 |
| 7 | Alexandrov et al. (2013), Cell Reports 3(1):246–259 | 1 | https://doi.org/10.1016/j.celrep.2012.12.008 | 2026-06-23 |
| 8 | Alexandrov et al. (2020), Nature 578:94–101; COSMIC SBS96 | 1–2 | https://cancer.sanger.ac.uk/signatures/sbs/sbs96/ | 2026-06-23 |

### A.2 Key Evidence Points (extension)

1. NMF model V ≈ W·H, W ≥ 0 (channels × k = signatures), H ≥ 0 (k × samples = exposures); blind-source-separation framing — Source 7.
2. Lee & Seung Frobenius multiplicative updates (Theorem 1): H ← H ⊙ (WᵀV) ⊘ (WᵀWH); W ← W ⊙ (VHᵀ) ⊘ (WHHᵀ) — Source 5/6.
3. Objective ‖V − WH‖²_F is monotonically non-increasing under these updates (Theorem 1) — Source 5/6.
4. Fixed point at exact factorization: multiplicative factors are all-ones when V = WH — Source 6.
5. NMF is non-convex → local minimum only; initialisation-dependent; fixed seed for reproducibility — Source 5/6.
6. Each signature is L1-normalised to sum to 1 (a probability distribution over the channels) — Source 7/8 (COSMIC).

### A.3 Canonical Methods Under Test (extension)

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ExtractSignatures(countMatrix, rank, maxIterations, tolerance, seed)` | OncologyAnalyzer | Canonical | NMF V≈WH via Lee & Seung Frobenius MU rules (Source 5/7) |
| `SignatureExtractionResult` (record) | OncologyAnalyzer | Internal | signatures (L1-normalised) + exposures + residual + iterations + objective history |

### A.4 Invariants (extension)

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-NMF-1 | All extracted signature weights ≥ 0 and all exposures ≥ 0 | Yes | Source 5/6 (MU preserve nonnegativity) |
| INV-NMF-2 | Each extracted signature is L1-normalised (Σ channels = 1) | Yes | Source 7/8 (COSMIC probability distribution) |
| INV-NMF-3 | Objective ‖V − WH‖²_F is monotonically non-increasing across iterations | Yes | Source 5/6 (Theorem 1) |
| INV-NMF-4 | Exactly factorable V = W₀·H₀ ⇒ residual → 0 at convergence | Yes | Source 6 (factors of ones at V=WH) |
| INV-NMF-5 | Planted signatures recovered up to column permutation & positive scaling (separable W₀,H₀) | Yes | Source 7 (blind source separation); separability uniqueness |

### A.5 Test Cases (extension)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Exact reconstruction | V = W₀·H₀, rank 2, tol 0 | FinalResidual < 1e-6 | Source 6 (V=WH fixed point) |
| M2 | Product recovers V | each (W·H)[c,s] ≈ V[c,s] | within 1e-3 | Source 7 (V≈WH) |
| M3 | Planted recovery | separable W₀,H₀ | each planted sig cosine ~1 with some extracted (perm/scale) | Source 7; separability |
| M4 | Nonnegativity | any input | min W ≥ 0 and min H ≥ 0 | Source 5/6 (INV-NMF-1) |
| M5 | L1 normalisation | any input | each signature sums to 1 | Source 7/8 (INV-NMF-2) |
| M6 | Monotonicity | objective history | non-increasing (within 1e-9 slack) | Source 5/6 (INV-NMF-3) |
| S1 | SBS-96 planted | 96 channels, k=2 | low relative residual; signatures 96-long, sum to 1 | Source 7/8 |
| S2 | Determinism | same seed twice | identical factors & iteration count | Source 5 (seeded init) |
| S3 | Scale absorption | exactly factorable | exposures total > 0 (scale absorbed, not zeroed) | Source 7 (W·H invariant) |
| V1–V10 | Validation | null/empty/zero-samples/ragged/negative/NaN/rank<1/rank>channels/maxIter≤0/tol<0 | throws ArgumentNullException / ArgumentException | documented failure modes |

### A.6 Coverage (extension)

Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ExtractSignatures_Tests.cs` (19 tests). All cases
M1–M6, S1–S3, V1–V10 ✅ Implemented and passing. Work queue remaining = 0.

### A.7 Assumptions (extension)

| # | Assumption | Used In |
|---|-----------|---------|
| 3 | NMF objective = squared Euclidean (Frobenius), Lee & Seung Theorem 1 (KL/Poisson variant not implemented) | M1–M6 |
| 4 | Deterministic seeded non-negative random initialisation (seed default 42) | S2 |

### A.8 Open Questions / Decisions (extension)

1. Decision: implement the **Frobenius** objective (clean monotone-non-increase property; faithful Lee & Seung
   Theorem 1). KL/Poisson variant documented as a not-implemented refinement.
2. Out of scope (residual limitation): automatic rank selection / model-stability selection (many NMF restarts
   choosing k by reproducibility + error, à la SigProfiler). Extraction is at a **caller-specified** k only.
3. Planted-truth recovery (M3) uses **separable** W₀, H₀ (pure-pixel condition) so the nonnegative factorisation
   is unique up to permutation/scaling, making ground truth recoverable; verification runs to tight convergence
   (tolerance 0) since the default relative-change stop halts before full stationarity.

---

## Appendix B — ONCO-SIG-002 enhancement: rank selection + KL/Poisson objective + cosine reference matching (2026-06-23)

Closes all three clauses of LIMITATIONS §2 ("Automatic NMF rank / model-stability selection … de-novo
extraction uses the Frobenius, not the Poisson/KL, objective … extracted signatures are not auto-matched to
COSMIC references").

### B.1 Authoritative Sources (enhancement)

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 9 | Lee & Seung (2001), Theorem 2 (KL/divergence updates) | 1 | https://arxiv.org/html/2501.11341v1 ; https://proceedings.neurips.cc/paper/2000/file/f9d1152547c0bde01830b7e8bd60024c-Paper.pdf | 2026-06-23 |
| 10 | Brunet et al. (2004), PNAS 101(12):4164 (consensus / cophenetic) | 1 | https://doi.org/10.1073/pnas.0308531101 (defns via https://nimfa.biolab.si/nimfa.models.nmf.html , https://search.r-project.org/CRAN/refmans/NMF/html/cophcor.html) | 2026-06-23 |
| 11 | Islam et al. (2022), SigProfilerExtractor, Cell Genomics 2(11):100179 | 1/3 | https://doi.org/10.1016/j.xgen.2022.100179 ; https://github.com/AlexandrovLab/SigProfilerExtractor | 2026-06-23 |
| 12 | Rousseeuw (1987), silhouette width, J. Comput. Appl. Math. 20:53 | 1 | https://doi.org/10.1016/0377-0427(87)90125-7 | 2026-06-23 |

### B.2 Key Evidence Points (enhancement)

1. KL/divergence objective `D(V‖WH)=Σ(V·log(V/WH)−V+WH)`; Theorem-2 updates `H_aμ←H_aμ·(Σ_i W_ia V_iμ/(WH)_iμ)/Σ_i W_ia`, `W_ia←W_ia·(Σ_μ H_aμ V_iμ/(WH)_iμ)/Σ_μ H_aμ`; monotone non-increasing — Source 9.
2. Connectivity matrix C_ij = 1 iff samples i,j share the argmax metagene; consensus = mean connectivity over runs; cophenetic correlation = Pearson(consensus-distance, average-linkage cophenetic distance); rank = first where cophenetic begins to fall — Source 10.
3. Stability = per-signature average silhouette width across replicate runs; select rank with average stability ≥ 0.80 and minimum ≥ 0.20, plus reconstruction error — Source 11; silhouette s=(b−a)/max(a,b) — Source 12.
4. De-novo signatures matched to references by maximising cosine similarity (Hungarian) — Source 11.

### B.3 Canonical Methods Under Test (enhancement)

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ExtractSignatures(countMatrix, rank, NmfObjective, …)` | OncologyAnalyzer | Canonical | KL/Poisson (Theorem 2) and Frobenius (Theorem 1) variants (Source 9) |
| `SelectRank(countMatrix, minRank, maxRank, …)` | OncologyAnalyzer | Canonical | cophenetic + silhouette stability + error rank selection (Source 10/11) |
| `MatchToReferenceSignatures(extracted, references)` | OncologyAnalyzer | Canonical | greedy best-cosine label per extracted signature (Source 11) |
| `RankStability`, `RankSelectionResult`, `SignatureMatch` (records) | OncologyAnalyzer | Internal | diagnostics / outputs |

### B.4 Invariants (enhancement)

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-KL-1 | KL divergence D(V‖WH) is monotonically non-increasing across iterations | Yes | Source 9 (Theorem 2) |
| INV-KL-2 | KL extraction recovers separable planted signatures (cosine ~1) | Yes | Source 7/9 |
| INV-RS-1 | Rank 1 ⇒ cophenetic correlation = 1.0 (all-ones consensus) | Yes | Source 10 |
| INV-RS-2 | Mean reconstruction error is non-increasing in rank on factorable data | Yes | Source 7/11 |
| INV-RS-3 | Selection is deterministic for a fixed base seed (derived seed sequence) | Yes | fixed seed convention |
| INV-MT-1 | Cosine match of a positively-scaled/exact copy of a reference = that reference, cosine 1.0 | Yes | Source 1/11 (cosine scale-invariance) |

### B.5 Test Cases (enhancement)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M-KL1 | KL monotonicity | D(V‖WH) history | non-increasing (within 1e-7 slack) | Source 9 |
| M-KL2 | KL planted recovery | separable k0=2 | both true sigs cosine > 0.999 | Source 7/9 |
| M-KL3 | KL normalisation | any input | nonneg; each signature sums to 1 | Source 8/9 |
| M-KL4 | KL determinism | same seed twice | identical factors | seeded init |
| M-KL5 | Frobenius backward-compat | default overload == explicit Frobenius | equal FinalResidual | this overload routing |
| M-RS1 | Rank selection true k | planted k0=2, ranks 1..3 | SelectedRank = 2; 3 diagnostic rows ascending | Source 10/11 |
| M-RS2 | Rank-1 cophenetic | rank 1 row | cophenetic = 1.0 | Source 10 |
| M-RS3 | Error vs rank | Frobenius, ranks 1..3 | err(2) ≤ err(1) | Source 7/11 |
| M-RS4 | Selection determinism | same seed | identical report | fixed seed |
| S-RS5 | Single candidate | minRank=maxRank=2 | SelectedRank=2; 1 row | range degenerate |
| M-MT1 | Cosine match scaled/unrelated | 5·ref0, ref1 | ref0 cosine 1.0; ref1 cosine 1.0; cross < 0.5 | Source 1/11 |
| M-MT2 | One match per extracted | exact copy of ref1 | ReferenceIndex 1, cosine 1.0 | Source 11 |
| V (enh) | Validation | null/minRank<1/maxRank<minRank/maxRank>channels/runs<1/threshold range; match null/mismatch/empty refs | throws | contract |

### B.6 Coverage (enhancement)

Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_SelectRank_Tests.cs` (21 tests). All cases
M-KL1..5, M-RS1..4, S-RS5, M-MT1..2, and validation guards ✅ Implemented and passing. Work queue remaining = 0.

### B.7 Assumptions (enhancement)

| # | Assumption | Used In |
|---|-----------|---------|
| 5 | Consensus clustering = argmax over H columns + average-linkage cophenetic distance (Brunet 2004 via nimfa/renozao-NMF) | M-RS1, M-RS2 |
| 6 | Per-signature stability = average silhouette width with cosine distance (SigProfiler "silhouette of the cluster"; Rousseeuw 1987) | M-RS1 |
| 7 | Matching = greedy best-cosine per extracted signature (per-signature reduction of SigProfiler Hungarian) | M-MT1, M-MT2 |

### B.8 Open Questions / Decisions (enhancement)

1. KL/Poisson is offered as a caller-selectable objective (`NmfObjective.KullbackLeibler`); the original
   `ExtractSignatures` Frobenius behaviour is preserved byte-for-byte via overload routing (M-KL5).
2. `SelectRank` defaults to the KL objective (the SigProfiler choice) but accepts either; per-rank cophenetic,
   average/minimum stability, and mean reconstruction error are all returned for auditability.
3. The reference set for `MatchToReferenceSignatures` is caller-supplied (COSMIC not embedded), per the existing
   `FitSignatures` convention; "auto-matching" = labelling each extracted signature with its closest reference.
