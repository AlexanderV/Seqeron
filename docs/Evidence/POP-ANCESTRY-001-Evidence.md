# Evidence Artifact: POP-ANCESTRY-001

**Test Unit ID:** POP-ANCESTRY-001
**Algorithm:** Ancestry Estimation (supervised / projection ADMIXTURE — EM with fixed reference allele frequencies)
**Date Collected:** 2026-06-13

---

## Online Sources

### Alexander, Novembre & Lange (2009) — "Fast model-based estimation of ancestry in unrelated individuals" (ADMIXTURE)

**URL:** https://faculty.eeb.ucla.edu/Novembre/AlexanderEtAl_GR_2009.pdf (open author copy of *Genome Research* 19(9):1655–1664, DOI 10.1101/gr.094052.109)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper)

**Retrieval method:** WebSearch "supervised ADMIXTURE ancestry estimation EM algorithm fixed allele frequencies likelihood model" → fetched the UCLA-hosted PDF, converted to text with `pdftotext -layout`, and read the Methods section verbatim.

**Key Extracted Points:**

1. **Data encoding (Methods, "A statistical model"):** "Let g_ij represent the observed number of copies of allele 1 at marker j of person i. Thus, g_ij equals 2, 1, or 0 accordingly, as i has genotype 1/1, 1/2, or 2/2 at marker j."
2. **Parameters:** "Population k contributes a fraction q_ik of individual i's genome. Allele 1 at SNP j has frequency f_kj in population k." Q = {q_ik} is I×K, F = {f_kj} is K×J.
3. **Log-likelihood (Equation 2):** `L(Q,F) = Σ_i Σ_j { g_ij · ln( Σ_k q_ik f_kj ) + (2 − g_ij) · ln( Σ_k q_ik (1 − f_kj) ) }` (up to an additive constant).
4. **Constraints:** "0 ≤ f_kj ≤ 1, q_ik ≥ 0, and Σ_k q_ik = 1."
5. **FRAPPE EM update for ancestry (Equation 4):** `q_ik^{n+1} = (1/2J) · Σ_j [ g_ij · a^n_ijk + (2 − g_ij) · b^n_ijk ]` where `a^n_ijk = q^n_ik f^n_kj / (Σ_m q^n_im f^n_mj)` and `b^n_ijk = q^n_ik (1 − f^n_kj) / (Σ_m q^n_im (1 − f^n_mj))`.
6. **Convergence criterion (Equation 5):** declare convergence once `L(Q^{n+1},F^{n+1}) − L(Q^n,F^n) < ε`; ADMIXTURE's default is `ε = 10^{-4}` ("we choose ε = 10⁴ [printed; = 10⁻⁴] as the default stopping criterion in ADMIXTURE"); FRAPPE used the looser ε = 1.
7. **Per-iteration complexity:** "The computational complexity for each iteration of this algorithm is O(IJK²)."
8. **Label invariance:** "the log-likelihood (Equation 2) is invariant under permutations of the labels of the ancestral populations" (K! equivalent global maxima).

### Alexander & Lange (2011) — "Enhancements to the ADMIXTURE algorithm for individual ancestry estimation" (supervised mode)

**URL:** https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-246 (*BMC Bioinformatics* 12:246, DOI 10.1186/1471-2105-12-246)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper)

**Retrieval method:** WebSearch "ADMIXTURE supervised learning mode .pop file allele frequencies fixed reference individuals manual" returned this article; the search result snippet quoted the supervised-mode semantics (see point 1). Full text is gated behind a Springer auth redirect (recorded as a retrieval limitation), so the supervised semantics below are corroborated by the official ADMIXTURE manual (next source), which was retrieved in full.

**Key Extracted Points:**

1. **Supervised analysis:** when reference individuals of known ancestry are available they are used as training samples; estimation of the admixed individuals' ancestries is then a supervised-learning problem with "less uncertainty in allele frequencies" and "shorter run times owing to the reduced number of parameters to estimate."

### ADMIXTURE 1.4 Software Manual — Alexander, Shringarpure, Novembre & Lange

**URL:** https://dalexander.github.io/admixture/admixture-manual.pdf
**Accessed:** 2026-06-13
**Authority rank:** 3 (canonical reference-implementation documentation)

**Retrieval method:** WebSearch returned the manual; fetched the PDF and converted with `pdftotext -layout`; read §2.10 (Supervised analysis) and §2.14 (Projection analysis) verbatim.

**Key Extracted Points:**

1. **Supervised mode (§2.10):** "it is not uncommon that some or all of the individuals in our data sample will have known ancestries, allowing us to set some rows in the matrix Q to known constants. … It is assumed that all reference samples have 100% ancestry from some ancestral population."
2. **Projection mode (§2.14):** allele frequencies learned from a reference panel are supplied as **fixed** input to estimate (project) ancestry of new samples: "Use learned allele frequencies as (fixed) input to next step … Run projection ADMIXTURE … `admixture -P study.bed 2`." This is exactly the "estimate Q given fixed F" task implemented here.

---

## Documented Corner Cases and Failure Modes

### From Alexander et al. (2009)

1. **Label permutation non-identifiability:** the likelihood has at least K! equivalent maxima (Eq. 2 is invariant under relabeling of ancestral populations). With *fixed, labelled* reference panels (this unit) the labels are pinned, so this is not a problem here — but it means results are only meaningful relative to the supplied reference labels.
2. **EM monotonicity / slow convergence:** "FRAPPE's EM algorithm converges slowly, as do many EM algorithms." The EM is an ascent algorithm — successive log-likelihoods are non-decreasing (basis of the Eq. 5 stopping rule).

### From ADMIXTURE manual

1. **Reference panels are 100% single-ancestry:** in supervised/projection mode each reference population's allele-frequency vector represents a pure ancestral source (§2.10).

---

## Test Datasets

### Dataset: Two-population symmetric panel, one diagnostic individual (derivation from Eq. 2 / Eq. 4)

**Source:** Derived directly from Alexander et al. (2009) Equations 2 and 4 (allele frequencies fixed). The EM was iterated by hand-equivalent computation following Eq. 4 literally (not from the implementation).

K = 2 populations (A, B); J = 2 SNPs. Fixed allele frequencies f_kj (allele-1 frequency):

| Population \ SNP | SNP1 | SNP2 |
|------------------|------|------|
| A | 0.8 | 0.2 |
| B | 0.2 | 0.8 |

Individual genotype g = [2, 0] (allele-1 counts). Start q⁰ = (0.5, 0.5).

| Iteration n | q_A | q_B | log-likelihood L |
|-------------|-----|-----|------------------|
| 0 | 0.5 | 0.5 | −2.772588722239781 |
| 1 | **0.8** | **0.2** | −1.542649923248 |
| 2 | 0.941176470588 | 0.058823529412 | −1.073055946379 |
| 3 | 0.984615384615 | 0.015384615385 | −0.938996389738 |
| ∞ | 1.0 | 0.0 | (→ −0.8929…, monotone increasing) |

After 1 EM iteration q = (0.8, 0.2) **exactly**; L is strictly increasing each step; q sums to 1 every iteration; the estimate converges to (1.0, 0.0).

### Dataset: Single-SNP closed-form one-iteration checks (derivation from Eq. 4)

**Source:** Eq. 4 with J = 1, K = 2, start q⁰ = (0.5, 0.5), f = [[0.9],[0.1]].

| Genotype g | One-iteration q (A, B) | Derivation |
|------------|------------------------|------------|
| g = 2 (1/1) | (0.9, 0.1) | denom_F = 0.5; q_A = ½·[2·(0.5·0.9/0.5)] = ½·1.8 = 0.9 |
| g = 1 (1/2) | (0.5, 0.5) | symmetric panels + heterozygote ⇒ a,b contributions cancel; stays uniform |
| g = 0 (2/2) | (0.1, 0.9) | by symmetry with g = 2 (1−f panel mirror) |

### Dataset: Uninformative identical panels

**Source:** Eq. 4. With f_A = f_B for every SNP, the numerators and denominators in a and b are proportional, so a uniform q is a fixed point: q⁰ = (0.5,0.5) → (0.5,0.5) for any genotype.

---

## Assumptions

1. **ASSUMPTION: Default iteration / convergence handling.** Alexander et al. specify the convergence rule (Eq. 5, ε = 10⁻⁴) but the public API exposes a fixed `maxIterations` budget. We run the EM (Eq. 4) for up to `maxIterations` and additionally stop early when the log-likelihood gain falls below ε = 10⁻⁴ (Eq. 5). This does not change the fixed-point/output for the tested cases (the EM converges to the same maximum); it only bounds work. Source-backed by Eq. 5; the *budget* parameter is API shape, not a correctness constant.
2. **ASSUMPTION: Missing-genotype encoding is excluded.** Genotypes are the integer allele-1 counts {0,1,2}. Any value outside {0,1,2} is treated as missing and contributes nothing to the Eq. 4 sums (its SNP is skipped for that individual). ADMIXTURE handles missing data but the manual does not give the per-individual EM term in closed form here; skipping is the standard treatment (a missing site provides no likelihood term in Eq. 2). Marked as assumption.

---

## Recommendations for Test Coverage

1. **MUST Test:** one EM iteration from uniform start on the symmetric two-population panel yields q = (0.8, 0.2) exactly. — Evidence: Alexander et al. (2009) Eq. 4 derivation (Test Dataset 1).
2. **MUST Test:** single-SNP g=2 with f=[0.9,0.1] → (0.9,0.1) after one iteration; g=0 → (0.1,0.9). — Evidence: Eq. 4 closed form (Test Dataset 2).
3. **MUST Test:** proportions sum to 1 (Σ_k q_ik = 1) for every returned individual. — Evidence: Alexander et al. (2009) constraint Σ_k q_ik = 1.
4. **MUST Test:** identical panels keep a uniform individual at uniform (uninformative fixed point). — Evidence: Eq. 4 (Test Dataset 3).
5. **MUST Test:** convergence — many iterations drive the diagnostic individual to its source population (q → (1,0)). — Evidence: Eq. 4 monotone ascent (Test Dataset 1).
6. **SHOULD Test:** log-likelihood (Eq. 2) is non-decreasing across iterations (EM ascent property / Eq. 5 basis). — Rationale: invariant from the EM theory.
7. **SHOULD Test:** empty individuals or empty reference panels → empty result; individual whose genotype length ≠ J is skipped. — Rationale: documented validation behavior.
8. **COULD Test:** label/order independence — permuting reference panels permutes the returned proportions consistently (Eq. 2 label invariance). — Rationale: identifiability note.

---

## References

1. Alexander DH, Novembre J, Lange K. 2009. Fast model-based estimation of ancestry in unrelated individuals. Genome Research 19(9):1655–1664. https://doi.org/10.1101/gr.094052.109 (open copy: https://faculty.eeb.ucla.edu/Novembre/AlexanderEtAl_GR_2009.pdf)
2. Alexander DH, Lange K. 2011. Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. BMC Bioinformatics 12:246. https://doi.org/10.1186/1471-2105-12-246
3. Alexander DH, Shringarpure SS, Novembre J, Lange K. ADMIXTURE 1.4 Software Manual (§2.10 Supervised analysis, §2.14 Projection analysis). https://dalexander.github.io/admixture/admixture-manual.pdf

---

## Change History

- **2026-06-13**: Initial documentation.
