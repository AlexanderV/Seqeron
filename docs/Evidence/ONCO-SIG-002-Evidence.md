# Evidence Artifact: ONCO-SIG-002

**Test Unit ID:** ONCO-SIG-002
**Algorithm:** Mutational Signature Fitting / Refitting (NNLS decomposition + cosine similarity)
**Date Collected:** 2026-06-14

---

## Online Sources

### MutationalPatterns (Blokzijl et al. 2018, Genome Medicine 10:33)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5922316/
**Accessed:** 2026-06-14 (retrieved via WebFetch of the PMC article page)
**Authority rank:** 3 (reference implementation in an established, peer-reviewed bioinformatics library)

**Key Extracted Points:**

1. **Cosine similarity (verbatim, section "Mutational profile similarity"):**
   `simAB = α = Σᵢ₌₁ⁿ Aᵢ Bᵢ / ( √(Σᵢ₌₁ⁿ Aᵢ²) · √(Σᵢ₌₁ⁿ Bᵢ²) )`, comparing two mutational profiles A and B
   as non-negative vectors with n mutation types; "The result ranges from 0 (independent) to 1 (identical)."
2. **NNLS refitting (verbatim, section "Finding the contribution of known signatures in mutation catalogues"):**
   `min_x ‖ S · x − d ‖₂²,  where x ≥ 0`, where S = signature matrix, x = signature contribution
   (weight) vector, d = the original 96-mutation count vector for a sample. Described as "a constrained
   version of the least-squares problem where the weights are not allowed to become negative," solved with
   an active-set method.
3. **Reconstruction quality:** the cosine similarity between the original profile and the reconstructed
   profile (S·x) is the quality check; profiles with similarity ≥ 0.95 indicate successful reconstruction.

### deconstructSigs (Rosenthal et al. 2016, Genome Biology 17:31)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4762164/
**Accessed:** 2026-06-14 (retrieved via WebFetch of the PMC article page)
**Authority rank:** 3 (reference implementation, peer-reviewed)

**Key Extracted Points:**

1. **Reconstruction model:** the reconstructed profile is `S·W` (signatures matrix times weights); the
   residual `R = T − (S·W)` between tumour profile T and reconstruction is minimised.
2. **Non-negativity (verbatim):** "any coefficient must be greater than 0, as negative contributions make
   no biological sense."
3. **Normalisation (verbatim):** "the weights W are normalized between 0 and 1," functionally producing
   proportions of signature contribution.
4. **Error metric (verbatim):** weights are chosen by "minimizing the sum-squared error (SSE) between the
   mutational profile of the tumor sample T and the mutational signature S."

### Non-negative least squares — Lawson & Hanson active-set algorithm (Wikipedia, citing Lawson & Hanson 1974)

**URL:** https://en.wikipedia.org/wiki/Non-negative_least_squares
**Accessed:** 2026-06-14 (retrieved via WebFetch; primary source: Lawson C.L., Hanson R.J. (1974),
*Solving Least Squares Problems*, Prentice-Hall, Ch. 23)
**Authority rank:** 4 (Wikipedia citing the primary Lawson-Hanson text; the algorithm itself is the cited primary)

**Key Extracted Points:**

1. **Problem statement (verbatim):** `argmin_x ‖ A x − y ‖₂²  subject to  x ≥ 0`.
2. **Active-set algorithm (verbatim steps):**
   - Initialise: P = ∅; R = {1,…,n}; **x** = 0; **w** = Aᵀ(**y** − A**x**).
   - Main loop while R ≠ ∅ and max(**w_R**) > ε: find j∈R at max(**w_R**); move j from R to P.
   - Let A_P = A restricted to P; set s_P = ((A_P)ᵀ A_P)⁻¹ (A_P)ᵀ **y**; s_R = 0.
   - Inner loop while min(s_P) ≤ 0: α = min( x_i / (x_i − s_i) ) over i∈P with s_i ≤ 0;
     **x** = **x** + α(**s** − **x**); move indices with x_j ≤ 0 from P to R; recompute s_P, s_R = 0.
   - Set **x** = **s**; update **w** = Aᵀ(**y** − A**x**).
   - Terminate when R = ∅ or max(**w_R**) ≤ ε.
3. **Equivalence:** when the unconstrained least-squares solution is already non-negative, NNLS returns it
   unchanged (no constraint is active).

### iMutSig (Pan & Wang 2020, F1000Research) — cosine similarity for signature comparison

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7702159/
**Accessed:** 2026-06-14 (retrieved via WebSearch result summary of the PMC article)
**Authority rank:** 3 (peer-reviewed application)

**Key Extracted Points:**

1. **Cosine similarity (verbatim form):** `CS(P,C) = P·C / (‖P‖ · ‖C‖)`, the cosine of the angle between
   two signature vectors; the standard metric for the pairwise similarity of mutational signature catalogues,
   value in [0, 1].

---

## Online Sources — De-novo NMF extraction (ONCO-SIG-002 extension, 2026-06-23)

### Lee D.D. & Seung H.S. (2001) — "Algorithms for Non-negative Matrix Factorization" (NIPS 13)

**URL:** https://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization
**Supplementary proof guide retrieved (full text, HTML):** https://arxiv.org/html/2501.11341v1
**Accessed:** 2026-06-23 (retrieved via WebFetch of the arXiv HTML proof guide)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points (verbatim formulas from the retrieved proof guide):**

1. **Euclidean (Frobenius) objective — Theorem 1:** "The Euclidean distance ‖V−WH‖ is nonincreasing under the
   update rules":
   - H update: `H_aμ ← H_aμ · (Wᵀ V)_aμ / (Wᵀ W H)_aμ`
   - W update: `W_ia ← W_ia · (V Hᵀ)_ia / (W H Hᵀ)_ia`
   "The objective remains constant only when W and H reach a stationary point."
2. **Generalized Kullback-Leibler objective — Theorem 2:** "The divergence D(V∥WH) is nonincreasing under the
   update rules":
   - H update: `H_aμ ← H_aμ · Σ_i [W_ia V_iμ /(WH)_iμ] / Σ_k W_ka`
   - W update: `W_ia ← W_ia · Σ_μ [H_aμ V_iμ /(WH)_iμ] / Σ_ν H_aν`
3. **Multiplicative structure / nonnegativity:** both updates are ratios applied to current values, "ensuring
   non-negativity preservation and monotonic convergence toward local optima."

### Non-negative matrix factorization — Wikipedia (citing the Lee & Seung primary)

**URL:** https://en.wikipedia.org/wiki/Non-negative_matrix_factorization
**Accessed:** 2026-06-23 (retrieved via WebFetch)
**Authority rank:** 4 (encyclopedic, citing the Lee & Seung primary)

**Key Extracted Points:**

1. **Frobenius updates (confirming the primary):** `H ← H .* (WᵀV) ./ (WᵀWH)` and `W ← W .* (VHᵀ) ./ (WHHᵀ)`,
   "applied iteratively until W and H converge."
2. **Objective:** `F(W,H) = ‖V − WH‖_F²`; "the multiplicative factors for W and H … are matrices of ones when
   V = WH" — a fixed point at exact factorization; the cost is non-increasing each iteration.
3. **Non-convexity (verbatim):** "current algorithms are sub-optimal in that they only guarantee finding a
   local minimum, rather than a global minimum of the cost function" because NMF is non-convex.

### Alexandrov L.B. et al. (2013) — "Deciphering Signatures of Mutational Processes Operative in Human Cancer", Cell Reports 3(1):246–259

**URL (abstract):** https://pubmed.ncbi.nlm.nih.gov/23318258/
**URL (publisher):** https://www.cell.com/cell-reports/fulltext/S2211-1247(12)00433-0 (DOI 10.1016/j.celrep.2012.12.008)
**Accessed:** 2026-06-23 (abstract retrieved via WebFetch of the PubMed page)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **NMF / blind-source-separation framing (from abstract):** "By modeling mutational processes as a blind
   source separation problem, we introduce a computational framework" that decomposes catalogs of somatic
   mutations into distinct mutational signatures and their respective contributions (exposures) across samples.
2. **Factor roles:** signatures are the latent sources (W); per-sample contributions are the mixing weights
   (H) — exactly the W (signatures) and H (exposures) factors of NMF for V ≈ WH.

### Alexandrov L.B. et al. (2020) — "The repertoire of mutational signatures in human cancer", Nature 578:94–101 (SigProfiler); COSMIC SBS96

**URL (COSMIC SBS96):** https://cancer.sanger.ac.uk/signatures/sbs/sbs96/
**Accessed:** 2026-06-23 (retrieved via WebFetch + WebSearch summary)
**Authority rank:** 1–2 (peer-reviewed paper / curated database)

**Key Extracted Points:**

1. **96 channels:** single-base substitutions are categorized into 96 SBS mutation types (6 × 4 × 4); matches
   `OncologyAnalyzer.Sbs96ChannelCount = 96` already in the file.
2. **Signature normalization:** each COSMIC SBS signature "is essentially a probability distribution across
   these 96 mutation types"; COSMIC reports signatures as "relative proportions of mutations" — i.e. each
   signature is L1-normalized to sum to 1. SigProfiler signature extraction is NMF-based.

---

## Documented Corner Cases and Failure Modes — De-novo NMF (2026-06-23)

### From Lee & Seung (2001) / Wikipedia

1. **Non-convexity → local minimum:** NMF is non-convex; multiplicative updates converge only to a local
   minimum/stationary point, dependent on initialization. Ground-truth recovery is therefore only guaranteed
   (up to column permutation/scaling) on data actually generated by a nonnegative factorization, and only for
   favourable initialisations. Tests use a fixed seed and planted-truth data.
2. **Fixed point at exact factorization:** when V = WH exactly the multiplicative factors are all-ones, so the
   iterates stop — the solver reconstructs an exactly-factorable V to within tolerance.
3. **Division-by-zero guard:** the denominators (WᵀWH, WHHᵀ) can be zero if a row/column collapses; a small
   epsilon is added to the denominator to avoid 0/0.

### From Alexandrov (2013/2020)

4. **Permutation/scale ambiguity:** NMF factors are identified only up to a permutation of the k components
   and a positive diagonal rescaling between W and H. L1-normalizing the signature columns (W) and absorbing
   the scale into H fixes the scaling ambiguity per the COSMIC convention.

---

## Test Datasets — De-novo NMF (2026-06-23)

### Dataset: Planted-truth small factorization (synthetic, exact)

**Source:** Lee & Seung (2001) fixed-point property (V = W₀H₀ ⇒ exact NMF reconstruction).

| Parameter | Value |
|-----------|-------|
| W₀ (channels × k) | non-negative, k = 2, channels = 4 |
| H₀ (k × samples) | non-negative, samples = 3 |
| V = W₀·H₀ | exactly factorable rank-2 nonnegative matrix |
| Expected reconstruction residual ‖V − WH‖_F | ≈ 0 (within iteration tolerance) |
| Expected recovery | planted signatures up to column permutation & scaling |

### Dataset: SBS-96 planted signatures (synthetic, exact)

**Source:** Alexandrov (2013) 96-channel model; planted W₀ with 96 channels, k = 2.

| Parameter | Value |
|-----------|-------|
| Channels | 96 (SBS-96) |
| k (rank) | 2 |
| Samples | 5 |
| V | W₀·H₀, exactly factorable |
| Expected | low residual; each extracted signature column sums to 1 after L1 normalization |

---

## Assumptions — De-novo NMF (2026-06-23)

1. **ASSUMPTION: Objective = squared Euclidean (Frobenius).** SigProfiler uses a Poisson/KL objective for the
   final extraction, but Lee & Seung give BOTH the Frobenius and KL multiplicative updates with proven monotone
   non-increase (Theorems 1 & 2). This implementation uses the **Frobenius** objective (Theorem 1) — a faithful,
   citable NMF with the clean monotone-non-increase property used for verification. The KL/Poisson variant is
   documented as a not-implemented refinement. Both are nonnegative factorisations of the same V ≈ WH model; the
   choice invents no constant.
2. **ASSUMPTION: Deterministic seeded initialisation.** Lee & Seung do not prescribe an initialisation;
   nonnegative random initialisation is standard. A fixed RNG seed is used so results are reproducible (mirrors
   the repo's existing `DefaultBootstrapSeed = 42` convention).

---

## Recommendations for Test Coverage — De-novo NMF (2026-06-23)

1. **MUST Test:** Exactly-factorable V = W₀·H₀ reconstructs with residual ‖V−WH‖_F ≈ 0 (within tolerance). —
   Evidence: Lee & Seung fixed-point property; Wikipedia "matrices of ones when V = WH".
2. **MUST Test:** Planted signatures recovered up to column permutation and positive scaling (cosine ≈ 1 after
   matching columns). — Evidence: Alexandrov (2013) blind-source-separation; NMF permutation/scale ambiguity.
3. **MUST Test:** Nonnegativity of W and H, and L1 column-normalization invariant (each signature sums to 1).
   — Evidence: Lee & Seung nonnegativity preservation; COSMIC "probability distribution across the 96 types".
4. **MUST Test:** Objective ‖V−WH‖_F² is monotonically non-increasing across iterations. — Evidence: Lee &
   Seung Theorem 1.
5. **SHOULD Test:** Reconstruction W·H approximates V; exposures absorb the scale removed by normalization. —
   Rationale: V ≈ WH model.
6. **SHOULD Test:** Input validation (null catalog, empty, ragged rows, k ≤ 0, k > channels, negative entries,
   maxIterations ≤ 0). — Rationale: documented failure modes / robustness.
7. **COULD Test:** Determinism — same seed ⇒ identical factors on repeated runs. — Rationale: reproducibility.

---

## Documented Corner Cases and Failure Modes

### From MutationalPatterns / NNLS (Lawson-Hanson)

1. **Zero observed catalog (d = 0):** the only minimiser of ‖S·x‖² with x ≥ 0 is x = 0 (all exposures 0);
   the reconstruction is the zero vector.
2. **Reconstruction-similarity gate:** cosine similarity < 0.95 between original and reconstructed profile
   signals unmeasured mutational processes (a quality warning, not an error).

### From the cosine-similarity definition

1. **Zero vector:** cosine similarity is undefined when either vector has zero norm (division by 0); a
   zero-norm input has no direction. Callers must treat a zero observed catalog as a degenerate case.
2. **Identical vectors:** cosine similarity = 1 exactly (parallel vectors).
3. **Orthogonal (disjoint support) vectors:** dot product 0 → cosine similarity 0.

---

## Test Datasets

### Dataset: Cosine similarity worked values (hand-derived from the MutationalPatterns / iMutSig formula)

**Source:** Blokzijl et al. (2018), formula simAB = ΣAB / (√ΣA² · √ΣB²); iMutSig CS(P,C) = P·C/(‖P‖‖C‖).

| Parameter | Value |
|-----------|-------|
| A = [1,2,3], B = [1,2,3] (identical) | sim = 14 / (√14 · √14) = 1.0 |
| A = [1,0], B = [0,1] (orthogonal) | sim = 0 / (1·1) = 0.0 |
| A = [1,1], B = [1,0] | sim = 1 / (√2 · 1) = 0.70710678118654752 |
| A = [3,4], B = [3,4] scaled B=[6,8] | sim = 50 / (5·10) = 1.0 (scale-invariant) |

### Dataset: NNLS decomposition worked values (hand-derived from min ‖Sx−d‖², x≥0)

**Source:** MutationalPatterns NNLS formula; Lawson-Hanson active-set algorithm.

| Parameter | Value |
|-----------|-------|
| S = [[1,0],[0,1]] (identity), d = [3,5] | x = [3,5]; reconstruction [3,5]; residual 0; cos(d,recon)=1 |
| S = [[1,1],[0,1]], d = [0,1] (unconstrained x₁=−1 negative) | NNLS x = [0, 0.5]; reconstruction [0.5,0.5] |
| S = [[1,0],[0,1]], d = [0,0] | x = [0,0]; reconstruction [0,0] |
| Normalised exposures of x=[3,5] (deconstructSigs proportions) | [0.375, 0.625] |

Derivation of the constraint-binding case S=[[1,1],[0,1]], d=[0,1]:
unconstrained normal equations SᵀS x = Sᵀd give [[1,1],[1,2]]x = [0,1] ⇒ x=[−1, 1]; x₁<0 so signature 1 is
clamped to 0; refit signature 2 alone: x₂ = (s₂·d)/(s₂·s₂) = ([1,1]·[0,1]) / ([1,1]·[1,1]) = 1/2 = 0.5.
NNLS solution = [0, 0.5].

---

## Assumptions

1. **ASSUMPTION: Exposure normalisation to proportions.** deconstructSigs "normalizes weights between 0 and
   1." Two normalisations are common: dividing by the sum of weights (proportions that sum to 1) and dividing
   by total mutations. We expose the raw NNLS exposures (counts) AND a proportion form = exposures / Σexposures
   (sums to 1 when Σ > 0, all-zero when Σ = 0). This is a presentation form; it does not change the fitted
   exposures, which are source-defined by the NNLS minimiser.
2. **ASSUMPTION: Cosine similarity of a zero vector.** No source defines cosine similarity for a zero-norm
   vector (division by 0). We return 0.0 for the pair when either norm is 0 (no shared direction / no signal),
   and document it; this only affects the degenerate empty-catalog case.

---

## Recommendations for Test Coverage

1. **MUST Test:** cosine similarity of identical, orthogonal, and the [1,1]/[1,0] (= 1/√2) vectors against
   exact hand-derived values; scale invariance. — Evidence: MutationalPatterns / iMutSig formula.
2. **MUST Test:** NNLS recovers exact exposures for an identity/orthonormal signature matrix; clamps a
   would-be-negative coefficient to 0 and refits (S=[[1,1],[0,1]], d=[0,1] ⇒ [0,0.5]); non-negativity of all
   exposures. — Evidence: MutationalPatterns NNLS formula; Lawson-Hanson algorithm.
3. **MUST Test:** reconstruction = S·x and reconstruction cosine = 1 for an exactly representable catalog;
   normalised exposures sum to 1 and equal the hand-derived proportions. — Evidence: deconstructSigs R=S·W.
4. **SHOULD Test:** zero catalog ⇒ all exposures 0, reconstruction 0; dimension-mismatch and null inputs
   throw. — Rationale: documented degenerate / failure modes.
5. **COULD Test:** invariant that residual SSE of the NNLS fit ≤ the SSE of the all-zero fit. — Rationale:
   NNLS is a minimiser, so it cannot be worse than x=0.

---

## References

1. Blokzijl F, Janssen R, van Boxtel R, Cuppen E (2018). MutationalPatterns: comprehensive genome-wide
   analysis of mutational processes. Genome Medicine 10:33. https://pmc.ncbi.nlm.nih.gov/articles/PMC5922316/
2. Rosenthal R, McGranahan N, Herrero J, Taylor BS, Swanton C (2016). deconstructSigs: delineating
   mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma
   evolution. Genome Biology 17:31. https://pmc.ncbi.nlm.nih.gov/articles/PMC4762164/
3. Lawson CL, Hanson RJ (1974). Solving Least Squares Problems. Prentice-Hall, Ch. 23 (active-set NNLS), via
   https://en.wikipedia.org/wiki/Non-negative_least_squares
4. Pan W, Wang X (2020). iMutSig: a web application to identify the most similar mutational signature using
   shiny. https://pmc.ncbi.nlm.nih.gov/articles/PMC7702159/
5. Lee DD, Seung HS (2001). Algorithms for Non-negative Matrix Factorization. Advances in Neural Information
   Processing Systems 13 (NIPS 2000). https://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization
   (supplementary proof guide: https://arxiv.org/html/2501.11341v1)
6. Non-negative matrix factorization. Wikipedia. https://en.wikipedia.org/wiki/Non-negative_matrix_factorization
7. Alexandrov LB, Nik-Zainal S, Wedge DC, Campbell PJ, Stratton MR (2013). Deciphering Signatures of Mutational
   Processes Operative in Human Cancer. Cell Reports 3(1):246–259. https://doi.org/10.1016/j.celrep.2012.12.008
8. Alexandrov LB et al. (2020). The repertoire of mutational signatures in human cancer. Nature 578:94–101.
   https://doi.org/10.1038/s41586-020-1943-3 ; COSMIC SBS96: https://cancer.sanger.ac.uk/signatures/sbs/sbs96/

---

## Online Sources — ONCO-SIG-002 enhancement: rank selection + KL/Poisson objective + reference matching (2026-06-23)

### Lee & Seung (2001) — KL (divergence/Poisson) multiplicative updates (Theorem 2)

**Retrieved (full text, HTML proof guide):** https://arxiv.org/html/2501.11341v1
**Primary PDF:** https://proceedings.neurips.cc/paper/2000/file/f9d1152547c0bde01830b7e8bd60024c-Paper.pdf
**Accessed:** 2026-06-23 (WebFetch of both; formulas extracted verbatim)
**Authority rank:** 1

**Key Extracted Points (verbatim):**

1. **KL (divergence) objective:** `D(V‖WH) = Σ_{ij} ( V_ij log(V_ij/(WH)_ij) − V_ij + (WH)_ij )`. This is the
   generalized KL divergence; it is the Poisson negative-log-likelihood up to a V-only constant — the objective
   SigProfiler uses for mutational-signature extraction.
2. **KL update for H (Theorem 2):** `H_aμ ← H_aμ · ( Σ_i W_ia V_iμ/(WH)_iμ ) / ( Σ_k W_ka )`.
3. **KL update for W (Theorem 2):** `W_ia ← W_ia · ( Σ_μ H_aμ V_iμ/(WH)_iμ ) / ( Σ_v H_av )`.
4. **Monotonicity:** "Both update sequences guarantee nonincreasing cost functions at each iteration" — so
   D(V‖WH) is monotone non-increasing under the KL updates (verification handle).

### Brunet, Tamayo, Golub & Mesirov (2004) — consensus clustering & cophenetic correlation (PNAS 101(12):4164)

**Primary:** https://www.pnas.org/doi/10.1073/pnas.0308531101 (PMC https://pmc.ncbi.nlm.nih.gov/articles/PMC384712/ — gated)
**Definitions retrieved verbatim from reference implementations citing Brunet (2004):**
https://nimfa.biolab.si/nimfa.models.nmf.html ; https://search.r-project.org/CRAN/refmans/NMF/html/cophcor.html
**Accessed:** 2026-06-23 (WebFetch)
**Authority rank:** 1 (paper) + 3 (impl docs for verbatim text)

**Key Extracted Points (verbatim):**

1. **Connectivity matrix:** "entry C_ij is 1 iff sample i and sample j belong to the same cluster, 0 otherwise.
   Sample assignment is determined by its largest metagene expression value." (cluster = argmax over the k
   exposure rows of that sample's H column.)
2. **Consensus matrix:** "the mean connectivity matrix across multiple runs of the factorization", used to
   "measure the stability of the clusters obtained by NMF."
3. **Cophenetic correlation coefficient:** the "Pearson correlation between the samples' distances induced by
   the consensus matrix (seen as a similarity matrix) and their cophenetic distances from a hierarchical
   clustering" (average linkage). "In a perfect consensus matrix, the coefficient equals 1."
4. **Rank-selection rule:** "Select the first rank, where the magnitude of the cophenetic correlation
   coefficient begins to fall."

### Alexandrov 2013 / SigProfilerExtractor (Islam et al. 2022) — stability (silhouette) + cosine matching

**URL:** https://www.cell.com/cell-reports/fulltext/S2211-1247(12)00433-0 ; https://github.com/AlexandrovLab/SigProfilerExtractor ;
https://www.cell.com/cell-genomics/fulltext/S2666-979X(22)00124-0
**Accessed:** 2026-06-23 (WebSearch summaries + WebFetch of the GitHub page)
**Authority rank:** 1 (paper) + 3 (reference implementation)

**Key Extracted Points (verbatim where quoted):**

1. **Stability via silhouette:** "stability is calculated for each signature by computing the silhouette value
   of the cluster corresponding to that signature"; partition clustering of the replicate factorizations is
   evaluated "in silhouette values."
2. **Selection rule:** "solutions stable if signatures have an average stability above 0.80 with no individual
   signature having stability below 0.20"; defaults `stability=0.8`, `min_stability=0.2`. "Model selection …
   based on each solution's stability and its reconstruction of the original data."
3. **Cosine matching:** "SigProfilerExtractor implements the Hungarian algorithm to pair consensus vectors from
   two matrices; the Hungarian algorithm maximizes the total cosine similarities of all paired vectors between
   two matrices." De-novo signatures are matched/decomposed to COSMIC references by cosine similarity.

### Rousseeuw (1987) — silhouette width definition

**URL:** https://doi.org/10.1016/0377-0427(87)90125-7
**Accessed:** 2026-06-23 (citation for the silhouette formula s(i) = (b−a)/max(a,b))
**Authority rank:** 1

**Key Extracted Points:** silhouette s(i) = (b(i) − a(i)) / max(a(i), b(i)), where a(i) is the mean
intra-cluster distance and b(i) the lowest mean distance to any other cluster; s ∈ [−1, 1].

---

## Documented Corner Cases — ONCO-SIG-002 enhancement (2026-06-23)

1. **Rank 1 cophenetic = 1 trivially:** with k = 1 every sample is in the single cluster, so every connectivity
   matrix is all-ones, the consensus is all-ones, and cophenetic correlation = 1 by construction (Brunet 2004).
   Rank 1 is therefore uninformative for stability-based selection.
2. **KL with zero entries:** `V_ij = 0` contributes `lim x→0 x log x = 0` to D; `(WH)_ij` must be floored to
   avoid log(·/0). `V_ij = 0` still adds `(WH)_ij` to D.
3. **Single candidate rank:** automatic selection over a range [kMin,kMax] degenerates to that rank.
4. **Cosine scale/permutation invariance:** matching a positively-scaled or channel-permuted copy of a
   reference returns cosine 1.0 to that reference (cosine is scale-invariant).

---

## Test Datasets — ONCO-SIG-002 enhancement (2026-06-23)

### Dataset: Planted-truth k0 = 2 (KL recovery + rank selection)

**Source:** Alexandrov (2013) blind-source-separation framing; V = W₀·H₀.

| Parameter | Value |
|-----------|-------|
| Channels | 4 |
| True signatures W₀ (L1) | s1 = [0.70,0.10,0.10,0.10], s2 = [0.10,0.10,0.70,0.10] |
| Samples | 6 (3 dominated by s1, 3 by s2, with mild mixing) |
| True rank k0 | 2 |
| Expected | KL recovers s1,s2 with cosine ≈ 1 (after matching); auto-rank reports per-k stability+error and selects k = 2 |

### Dataset: Cosine reference-matching cases

**Source:** cosine scale-invariance (Blokzijl 2018); SigProfiler cosine matching.

| Parameter | Value |
|-----------|-------|
| references | refA = [0.7,0.1,0.1,0.1], refB = [0.1,0.1,0.7,0.1] |
| query = 5·refA (scaled) | best match = refA, cosine = 1.0 |
| query = refB | best match to refA has low cosine (< 0.5), best overall = refB cosine = 1.0 |

---

## Assumptions — ONCO-SIG-002 enhancement (2026-06-23)

1. **ASSUMPTION: consensus clustering = argmax over H columns + average-linkage cophenetic distance.** Taken
   verbatim from the nimfa / renozao-NMF reference implementations of Brunet (2004) (PNAS PMC full text gated).
   Authority rank 3, consistent across both implementations.
2. **ASSUMPTION: per-signature stability = average silhouette width (cosine distance) of the run-signatures
   clustered to that consensus signature.** SigProfiler states "silhouette value of the cluster corresponding
   to that signature"; the silhouette is the standard Rousseeuw (1987) definition; cosine distance = 1 − cosine
   similarity is the signature-space metric used throughout this class.
3. **ASSUMPTION: matching = greedy best-cosine per extracted signature** (label each with its single closest
   reference + the cosine), the per-signature reduction of SigProfiler's global Hungarian assignment, which is
   exactly the "label each extracted signature with its closest reference and the cosine" contract requested.

---

## Recommendations for Test Coverage — ONCO-SIG-002 enhancement (2026-06-23)

1. **MUST:** KL update monotonically decreases D(V‖WH) each iteration — Lee & Seung (2001) Theorem 2.
2. **MUST:** KL extraction recovers planted s1/s2 (cosine ≈ 1 after matching) on k0 = 2 data — Alexandrov 2013.
3. **MUST:** auto rank selection picks k0 = 2 and reports per-k stability + reconstruction error — Brunet 2004 +
   Alexandrov 2013/SigProfiler.
4. **MUST:** cophenetic = 1.0 at rank 1 (all-ones consensus) — Brunet 2004.
5. **MUST:** cosine matching maps scaled/permuted reference to itself (cosine = 1.0) and unrelated to low cosine
   — cosine scale-invariance; SigProfiler.
6. **SHOULD:** determinism (same fixed seed-sequence ⇒ identical KL result, identical selection report).
7. **SHOULD:** validation (null/empty/ragged/negative; rank/k-range; runs < 1; empty reference) throws.

---

## Change History

- **2026-06-14**: Initial documentation (NNLS refitting).
- **2026-06-23**: Added de-novo NMF signature extraction (Lee & Seung Frobenius MU rules; Alexandrov NMF
  framing; COSMIC L1 signature normalization) for the ONCO-SIG-002 extension.
- **2026-06-23**: Added ONCO-SIG-002 enhancement evidence — Lee & Seung KL/Poisson Theorem 2 updates, Brunet
  (2004) consensus/cophenetic rank selection, Alexandrov 2013 / SigProfiler silhouette stability + cosine
  reference matching (closing all three LIMITATIONS clauses).
