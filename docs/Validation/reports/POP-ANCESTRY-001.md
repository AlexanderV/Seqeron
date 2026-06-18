# Validation Report: POP-ANCESTRY-001 — Ancestry Estimation (supervised / projection ADMIXTURE)

- **Validated:** 2026-06-15   **Area:** PopGen
- **Canonical method(s):** `PopulationGeneticsAnalyzer.EstimateAncestry(individuals, referencePops, maxIterations)` (private helpers `EstimateIndividualAncestry`, `AncestryLogLikelihood`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm

- **Alexander, Novembre & Lange (2009), *Genome Research* 19(9):1655–1664 (ADMIXTURE)** — fetched the open author copy
  (https://faculty.eeb.ucla.edu/Novembre/AlexanderEtAl_GR_2009.pdf), saved as PDF, converted with `pdftotext -layout`,
  and read the "A statistical model" / "FRAPPE's EM algorithm" Methods sections **verbatim**. Confirmed line-for-line:
  - **Genotype encoding:** "Let g_ij represent the observed number of copies of allele 1 at marker j of person i.
    Thus, g_ij equals 2, 1, or 0 accordingly, as i has genotype 1/1, 1/2, or 2/2 at marker j."
  - **Eq. 2 (log-likelihood):** `L(Q,F) = Σ_i Σ_j { g_ij ln(Σ_k q_ik f_kj) + (2−g_ij) ln(Σ_k q_ik (1−f_kj)) }` "up to an additive constant."
  - **Eq. 4 (FRAPPE EM ancestry update):** `q_ik^{n+1} = (1/2J) Σ_j [ g_ij a^n_ijk + (2−g_ij) b^n_ijk ]`,
    with `a^n_ijk = q^n_ik f^n_kj / (Σ_m q^n_im f^n_mj)` and `b^n_ijk = q^n_ik(1−f^n_kj)/(Σ_m q^n_im(1−f^n_mj))`.
  - **Eq. 5 (convergence):** declare convergence once `L(Q^{n+1},F^{n+1}) − L(Q^n,F^n) < ε`; "we choose ε = 10⁻⁴ as the
    default stopping criterion in ADMIXTURE" (printed as "104" in the OCR = 10⁻⁴; FRAPPE used ε = 1).
  - **Constraints:** "0 ≤ f_kj ≤ 1, q_ik ≥ 0, and Σ_k q_ik = 1."
  - **Label invariance:** "the log-likelihood (Equation 2) is invariant under permutations of the labels of the ancestral
    populations … at least K! equivalent global maxima."
  - **Complexity:** "O(IJK²)" per iteration (joint Q+F estimation).
- **Alexander & Lange (2011), *BMC Bioinformatics* 12:246** and **ADMIXTURE 1.4 Manual §2.10/§2.14** — cited in the
  Evidence doc for the supervised/projection semantics (F fixed, only Q estimated). Full text gated; the supervised
  reduction (estimate Q given fixed F) is exactly what Eq. 4 computes when F is not updated, which the 2009 paper fully
  specifies. Acceptable corroboration.

### Formula check

Every formula in `Ancestry_Estimation.md` §2.2, the TestSpec §1.2, and the Evidence doc reproduces the source equations
exactly (symbols, 1/2J normalization, a/b definitions, ε = 10⁻⁴, constraints). No divergence.

### Edge-case semantics check

Documented and sourced/derived: empty inputs → empty; genotype length ≠ J → individual skipped (ASM-04 index
alignment); genotype ∉ {0,1,2} → missing, SNP contributes no Eq. 2 term (standard treatment, marked ASSUMPTION 2);
all-missing → uniform prior 1/K; identical panels → uniform fixed point (INV-04). The two ASSUMPTIONs (maxIterations
budget + Eq. 5 early stop; missing-genotype skip) are honestly flagged as API-shape / standard-treatment, not
correctness constants.

### Independent cross-check (numbers)

Re-implemented Eq. 4 + Eq. 2 from scratch in Python this session and reproduced **every** Evidence/TestSpec value:

| Case | Source-derived expected | My independent computation |
|------|------------------------|----------------------------|
| M1 iter1 (f symmetric, g=[2,0]) | q=(0.8, 0.2) | (0.8, 0.19999999999999998) ✓ |
| M1 L₀,L₁,L₂,L₃ | −2.7725887, −1.5426499, −1.0730559, −0.9389963 | identical to 7+ dp ✓ |
| M2 single SNP g=2, f=(0.9,0.1) | (0.9, 0.1) | (0.9, 0.1) ✓ |
| M3 single SNP g=0 | (0.1, 0.9) | (0.0999999…, 0.9) ✓ |
| M7 heterozygote g=1 | (0.5, 0.5) | (0.5, 0.5) ✓ |
| M5 convergence (1000 iters) | (1.0, 0.0) | (1.0, 0.0) ✓ |
| M6 identical panels | (0.5, 0.5) | (0.5, 0.5) ✓ |
| S6 missing g=(-1,2), f=(0.3,0.9)/(0.7,0.1) | (0.9, 0.1) | (0.9, 0.1) ✓ |

### Findings / divergences

None. The biology and mathematics are correct and faithfully match the cited primary source.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs`
- `EstimateAncestry` (L1263–1290): materializes inputs; empty inds/refs → `yield break`; skips individuals whose
  genotype length ≠ panel SNP count; emits one `AncestryProportion` keyed by population id.
- `EstimateIndividualAncestry` (L1296–1361): uniform init 1/K; per iteration computes mix1 = Σ q_m f_mj,
  mix0 = Σ q_m(1−f_mj), accumulates `g·a + (2−g)·b` per population, divides by 2·J (J = informative SNP count);
  out-of-range g skipped; if no informative SNP → break (keeps uniform); early stop when ΔL < 1e-4.
- `AncestryLogLikelihood` (L1367–1396): Eq. 2 over informative SNPs.
- `AncestryLogLikelihoodTolerance = 1e-4` (L1231) = Eq. 5 default ε.

### Formula realised correctly?

Yes — code is a line-by-line transcription of Eq. 4 (the `a`/`b` responsibilities, the 1/2J normalizer) and Eq. 2.
The 1/2J normalization is what guarantees Σ_k q_ik = 1 exactly (INV-01) without renormalization, matching the source.
F is never mutated → genuinely the supervised/projection case.

### Cross-verification table recomputed vs code

Ran the 15 existing tests + 2 added: all pass with the exact source-derived values above. M5 with the impl's Eq. 5
early stop halts at iter 9 with q_A = 0.99999618 — within the test's 1e-3 bound (the exact MLE 1.0 is the limit, so a
tolerance, not an exact value, is correct at finite iterations).

### Variant/delegate consistency

Single public entry point; no `*Fast`/delegate variants. `AncestryProportion` record echoes `IndividualId` (verified by S4).

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** M1–M3, M6, M7, S6 assert exact values derived from Eq. 4 / Eq. 2 with `Within(1e-10)`.
  M1 (0.8 exact) is **discriminating**: the documented wrong per-population-responsibility update yields 0.941 at
  iter 1 (TestSpec §7), so M1 would fail against it. Not a tautology.
- **No green-washing:** no weakened assertions, no widened tolerances, no skips. M5's `Within(1e-3)` and S1's
  `GreaterThanOrEqualTo` are justified (asymptotic limit; monotonicity property — no exact finite-iter value exists).
- **Cover all logic:** all 7 MUST, 6 SHOULD, 1 COULD covered. **Two coverage gaps closed this session** (test-only,
  no code change): added `EstimateAncestry_ArbitraryInput_ProportionsInUnitInterval` (locks INV-02: 0 ≤ q ≤ 1, which
  previously had no standalone test) and `EstimateAncestry_ZeroIterations_ReturnsUniformPrior` (locks the
  `maxIterations = 0` boundary of the documented `maxIterations ≥ 0` contract).
- **Honest green:** FULL unfiltered suite = **6550 passed, 0 failed** (1 pre-existing benchmark skip, unrelated);
  `dotnet build` 0 errors (4 warnings are pre-existing, in `ApproximateMatcher_EditDistance_Tests.cs`, untouched).

**Gate result: PASS.**

### Findings / defects

No algorithm defect. Two minor test-coverage gaps (INV-02 bound, maxIterations=0 boundary) were completely closed by
adding two source/contract-locked tests. No assertion weakened.

## Verdict & follow-ups

- **Stage A: PASS** — description matches Alexander et al. (2009) Eq. 2/4/5 and the supervised/projection semantics
  verbatim; all corner cases sourced or honestly flagged as assumptions.
- **Stage B: PASS** — implementation is a faithful transcription; all expected values independently reproduced;
  tests are exact and discriminating.
- **End-state: CLEAN.** No defects; the two added tests strengthen coverage without changing behaviour.
