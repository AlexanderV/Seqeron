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

---

## Change History

- **2026-06-14**: Initial documentation.
