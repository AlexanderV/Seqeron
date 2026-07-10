---
type: source
title: "Evidence: POP-ANCESTRY-001 (Supervised / projection ADMIXTURE — ancestry proportions by EM with fixed reference allele frequencies)"
tags: [validation, population-genetics]
doc_path: docs/Evidence/POP-ANCESTRY-001-Evidence.md
sources:
  - docs/Evidence/POP-ANCESTRY-001-Evidence.md
source_commit: 9e7930d3c6f0f119ea8d74d1a72b1581f0850ac4
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: POP-ANCESTRY-001

The validation-evidence artifact for test unit **POP-ANCESTRY-001** — **Ancestry Estimation**
in the *supervised / projection ADMIXTURE* mode: estimate an individual's ancestry proportions
**Q** by EM given **fixed** reference-panel allele frequencies **F**. It is one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the model,
EM update, invariants, worked oracles, and corner cases are synthesized in the dedicated concept
[[ancestry-estimation-admixture]]. This is the **first POP-\* (population-genetics) unit** in the
wiki. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Alexander, Novembre & Lange (2009), "Fast model-based estimation of ancestry in unrelated
    individuals"** (*Genome Research* 19(9):1655–1664, DOI 10.1101/gr.094052.109; rank 1) — the
    ADMIXTURE model. Genotype `g_ij ∈ {0,1,2}` = allele-1 copy count; parameters `Q = {q_ik}`
    (I×K ancestry fractions) and `F = {f_kj}` (K×J allele-1 frequencies); the **log-likelihood
    (Eq. 2)** `L(Q,F) = Σ_i Σ_j { g_ij·ln(Σ_k q_ik f_kj) + (2−g_ij)·ln(Σ_k q_ik (1−f_kj)) }`
    under constraints `0 ≤ f_kj ≤ 1, q_ik ≥ 0, Σ_k q_ik = 1`; the **FRAPPE EM update (Eq. 4)**
    `q_ik^{n+1} = (1/2J)·Σ_j [ g_ij·a^n_ijk + (2−g_ij)·b^n_ijk ]` with
    `a^n_ijk = q^n_ik f^n_kj / (Σ_m q^n_im f^n_mj)` and
    `b^n_ijk = q^n_ik (1−f^n_kj) / (Σ_m q^n_im (1−f^n_mj))`; the **convergence rule (Eq. 5)**
    stop when `ΔL < ε` (ADMIXTURE default `ε = 10⁻⁴`; FRAPPE used ε = 1); per-iteration cost
    **O(IJK²)**; and the **label-permutation invariance** (Eq. 2 unchanged under relabeling the K
    ancestral populations → K! equivalent maxima).
  - **Alexander & Lange (2011), "Enhancements to the ADMIXTURE algorithm for individual ancestry
    estimation"** (*BMC Bioinformatics* 12:246, DOI 10.1186/1471-2105-12-246; rank 1) —
    **supervised** mode: reference individuals of known ancestry become training samples, so
    estimating the admixed individuals' Q is a supervised-learning problem with "less uncertainty
    in allele frequencies" and shorter run times (fewer parameters). Full text was gated behind a
    Springer auth redirect (recorded retrieval limitation); the supervised semantics are
    corroborated by the manual below.
  - **ADMIXTURE 1.4 Software Manual** (Alexander, Shringarpure, Novembre & Lange; rank 3) —
    **§2.10 Supervised analysis**: known-ancestry rows of Q are pinned to constants and "all
    reference samples have 100% ancestry from some ancestral population." **§2.14 Projection
    analysis**: allele frequencies learned from a reference panel are supplied as **fixed** input
    to estimate (project) Q of new samples (`admixture -P study.bed 2`) — exactly the "estimate Q
    given fixed F" task this unit implements.
- **Datasets (documented oracles):**
  - **Two-population symmetric panel, one diagnostic individual** (derived from Eq. 2 / Eq. 4,
    F fixed): K = 2 (A, B), J = 2 SNPs, `f_A = [0.8, 0.2]`, `f_B = [0.2, 0.8]`; genotype
    `g = [2, 0]`; start `q⁰ = (0.5, 0.5)`. After **one** EM iteration `q = (0.8, 0.2)` **exactly**
    (L = −1.5426…); iteration 2 → (0.9412, 0.0588); iteration 3 → (0.9846, 0.0154); converges to
    `(1.0, 0.0)`. L is strictly increasing each step and q sums to 1 every iteration.
  - **Single-SNP closed-form one-iteration checks** (Eq. 4, J = 1, K = 2, `q⁰ = (0.5,0.5)`,
    `f = [[0.9],[0.1]]`): `g = 2 → (0.9, 0.1)`; `g = 1` (heterozygote on a symmetric panel) →
    `(0.5, 0.5)` (a,b contributions cancel, stays uniform); `g = 0 → (0.1, 0.9)`.
  - **Uninformative identical panels** (`f_A = f_B` at every SNP): a uniform q is a fixed point —
    `(0.5,0.5) → (0.5,0.5)` for any genotype (numerators/denominators in a and b are
    proportional).

## Deviations and assumptions

Two documented **assumptions** (neither changes the fixed point / output for the tested cases):

1. **Default iteration / convergence handling.** Eq. 5 specifies `ΔL < ε` (ε = 10⁻⁴), but the
   public API exposes a fixed `maxIterations` budget. The unit runs Eq. 4 for up to
   `maxIterations` and additionally stops early once the log-likelihood gain falls below ε = 10⁻⁴.
   The budget is API shape, not a correctness constant; the EM converges to the same maximum.
2. **Missing-genotype encoding is excluded.** Genotypes are integer allele-1 counts {0,1,2}; any
   value outside {0,1,2} is treated as missing and contributes nothing to the Eq. 4 sums (that SNP
   is skipped for the individual). ADMIXTURE handles missing data but the manual does not give the
   per-individual EM term in closed form; skipping is the standard treatment (a missing site
   provides no Eq. 2 likelihood term).

**Identifiability note (from Eq. 2):** the likelihood has ≥ K! equivalent maxima under
population-label permutation. With the **fixed, labelled** reference panels of this unit the labels
are pinned, so non-identifiability does not bite here — but results are meaningful only relative to
the supplied reference labels. **EM ascent:** successive log-likelihoods are non-decreasing (the
basis of the Eq. 5 stopping rule), though EM converges slowly.

Recommended coverage (MUST): one EM iteration on the symmetric panel → `q = (0.8, 0.2)` exactly;
single-SNP `g=2/f=[0.9,0.1] → (0.9,0.1)` and `g=0 → (0.1,0.9)`; `Σ_k q_ik = 1` for every returned
individual; identical panels keep a uniform individual uniform; many iterations drive the
diagnostic individual to its source population `q → (1,0)`. SHOULD: L (Eq. 2) non-decreasing across
iterations; empty individuals / empty reference panels → empty result; genotype length ≠ J →
individual skipped. COULD: label/order independence (permuting reference panels permutes returned
proportions consistently). No source contradictions — Alexander et al. (2009), Alexander & Lange
(2011), and the ADMIXTURE manual are mutually consistent; the manual's supervised/projection modes
are the applied specialization of the 2009 model.
