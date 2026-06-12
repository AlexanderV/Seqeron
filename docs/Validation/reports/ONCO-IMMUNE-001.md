# Validation Report: ONCO-IMMUNE-001 — Immune Infiltration Estimation

- **Validated:** 2026-06-12   **Area:** Oncology / Tumor Immunology
- **Canonical method(s):** `ImmuneAnalyzer.EstimateInfiltration(...)` (ssGSEA + ESTIMATE tumor purity), `ImmuneAnalyzer.DeconvoluteImmuneCells(...)` (NNLS cell-type deconvolution)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End state:** CLEAN (no defect; honest scope confirmed)

This unit is **not** TMB/MSI/neoantigen. It computes two tumor-immunology metrics:
1. **ESTIMATE-style immune/stromal infiltration** via single-sample GSEA (ssGSEA) integral form + tumor-purity estimate.
2. **Immune cell-type deconvolution** via Non-Negative Least Squares (NNLS), distinct from CIBERSORT's ν-SVR.

---

## Stage A — Description

### Sources opened & what they confirm
- **Yoshihara et al. (2013), Nat Commun 4:2612 (ESTIMATE)** — confirmed via Nature page, hacksig/CRAN refman, and Aging-US methods: the tumor-purity formula is exactly
  **Tumor_purity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)**.
  Both coefficients match the code constants `EstimatePurityCoefficientA/B` to the last digit. External sources also note this calibration was derived from **Affymetrix** data and is not directly transferable to RNA-seq ESTIMATE scores (scope note, see below).
- **Barbie et al. (2009), Nature 462:108–112 + Hänzelmann et al. (2013), GSVA, BMC Bioinformatics 14:7** — confirmed via GSVA Bioconductor vignette and the GenePattern/MSigDB ssGSEA module docs: ssGSEA computes a per-sample enrichment score as the **(weighted) sum of the difference between the empirical CDFs of genes inside vs. outside the gene set** (the *integral/area* of the running-sum statistic, not the max deviation as in classic GSEA), with the tail-weight exponent **τ = 0.25** for `method="ssgsea"`. The hit positions are weighted by rank^τ. The code's `SsGseaTau = 0.25` and the integral-of-running-sum construction match this.
- **Newman et al. (2015), Nat Methods 12:453 (CIBERSORT) + Abbas et al. (2009), PLoS One 4:e6098 + Lawson & Hanson (1995)** — confirmed: the deconvolution model is **m = S·f**, solved as **min ‖m − S·f‖² s.t. f ≥ 0**, then normalized so **Σf = 1**. NNLS (Lawson–Hanson active-set) is the LLSR/NNLS baseline benchmarked in Newman (2015); CIBERSORT itself uses ν-SVR. The 22 LM22 cell phenotypes are correctly enumerated.

### Formula check
| Claim | Source | Code | Match |
|---|---|---|---|
| Purity = cos(a + b·ESTIMATE), a=0.6049872018, b=0.0001467884 | Yoshihara 2013 | `ComputeTumorPurity` | ✅ exact |
| ESTIMATE score = Immune + Stromal | Yoshihara 2013 | `estimateScore = immuneScore + stromalScore` | ✅ |
| ssGSEA τ = 0.25, integral of running sum, rank^τ hits | Barbie 2009 / GSVA | `ComputeSsGseaScore` | ✅ |
| NNLS min‖m−Sf‖² s.t. f≥0, then Σf=1 | Lawson–Hanson / Abbas 2009 | `SolveNnls` + normalization | ✅ |

### Edge-case semantics (all sourced/defined)
- Empty profile → ImmuneScore=0, StromalScore=0, TumorPurity=cos(a). ✅
- No overlapping genes → scores 0 (ssGSEA guards `nHits==0` / `nMiss==0`). ✅
- No overlapping deconvolution genes → all fractions 0, OverlappingGenes=0. ✅
- Purity formula can exceed [0,1] at extreme scores → clamped via `Math.Clamp`. ✅ (matches documented corner case)
- Null → `ArgumentNullException`. ✅

### Independent cross-check (hand computation, ssGSEA integral)
Profile {A=100, B=1, C=0.5}, gene set {A,C}. Ranked desc: A(rank3), B(rank2), C(rank1). N=3, N_H=2, N_miss=1, missStep=1, TW=3^¼+1.
Walk: hit A → RS=3^¼/TW; miss B → RS=3^¼/TW−1; hit C → RS=(3^¼+1)/TW−1=0.
Integral = 2·3^¼/TW − 1 = **(3^¼ − 1)/(3^¼ + 1) ≈ 0.13655**. This is the value asserted by test M14a, and it is distinct from the ≈0.5799 an expression-weighted variant would give. Single-hit-at-top → integral 1.5 (M14b); single-hit-at-bottom → −1.5 (M14c). All three independently recomputed and correct.

### Findings / divergences (Stage A notes)
1. **Honest-scope note (not a defect):** The default signature matrix is a **simplified 5-marker × 22-cell-type** matrix, not the full 547-gene LM22; deconvolution is **NNLS, not ν-SVR (CIBERSORT)**; and the ESTIMATE purity coefficients are **Affymetrix-calibrated** while this implementation's single-sample ssGSEA integral is on a different (much smaller) numeric scale than the cohort-scaled ESTIMATE score. Consequently the absolute TumorPurity number is **not clinically meaningful**. This is *honestly declared*: the XML docs and the spec's Assumption Register (assumptions 1–3) explicitly state these are simplifications and that users should supply their own signature matrices / data for production use. There is **no** "clinical/diagnostic/FDA/validated" advertising in the source (`grep` confirmed none). Tests assert only the *formula identities* and *mathematical invariants* (linearity, Σf=1, [0,1] clamp), never clinical accuracy. → PASS-WITH-NOTES, correctly scoped.

---

## Stage B — Implementation

- **Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs`
  - `EstimateInfiltration` (348–378), `ComputeSsGseaScore` (513–574), `ComputeTumorPurity` (580–584)
  - `DeconvoluteImmuneCells` (409–496), `SolveNnls` (591–686), passive-set LSQ + Gaussian elimination (719–837), Pearson/RMSE helpers (842–884).

### Formula realised correctly?
- **ssGSEA:** ranks descending, computes `totalHitWeight = Σ rank^τ`, then walks accumulating `rank^τ/totalHitWeight` on hits (Phit increments, cumulative→1) and `−1/nMiss` on misses (Pmiss decrements, cumulative→1), summing the running difference into `integral`. This is exactly the integral of (Phit − Pmiss). ✅
- **Purity:** `cos(a + b·estimateScore)` clamped to [0,1]. ✅
- **NNLS:** Lawson–Hanson active-set with gradient `w = Aᵀ(b − Ax)`, passive-set normal-equations solve via Gaussian elimination w/ partial pivoting, feasibility back-off via α-ratio, final non-negativity cleanup, then Σf=1 normalization with a `sum > 0` div guard. ✅

### Cross-verification table recomputed vs code (tests run)
| ID | Expected (sourced) | Code result | Match |
|----|--------------------|-------------|-------|
| M1 | scores 0, purity=cos(a) | ✅ | ✅ |
| M5 | pure CD8 → f_CD8=1.0, corr=1, RMSE=0 | ✅ | ✅ |
| M6 | B+CD8 50:50 → 0.5/0.5 | ✅ | ✅ |
| S1 | 75:25 → 0.75/0.25 | ✅ | ✅ |
| M9/INV-3 | purity=cos(a+b·score) ∈ [0,1] | ✅ | ✅ |
| M10/INV-4 | ESTIMATE = immune+stromal | ✅ | ✅ |
| M14a | (3^¼−1)/(3^¼+1)≈0.13655 | ✅ | ✅ |
| M14b/c | +1.5 / −1.5 integral | ✅ | ✅ |
| INV-1/2 | f≥0, Σf=1 (6 cell types) | ✅ | ✅ |
| C1/C3 | 22 cell types, 5 genes each | ✅ | ✅ |

### Variant/delegate consistency
Two canonical static methods, no `*Fast`/instance variants to reconcile. Constants in the test file mirror the source constants exactly.

### Numerical robustness
- Div guards: ssGSEA returns 0 on `nHits==0 || nMiss==0` and `totalHitWeight==0`; normalization guarded by `sum > 0`; Pearson denominator `< 1e-15 → 0`; Gaussian pivot `< 1e-15` skipped; RMSE `n==0 → 0`. No div-by-zero on stated ranges. Negative (log-transformed) expression handled (C2/C2b). ✅

### Test quality audit
33 tests, all asserting **exact sourced values** with tight tolerances (1e-10 identities, 1e-6 computed). M14a/b/c are genuine discriminating references (rank-vs-expression weighting; integral-vs-max-deviation). No tautological "does not throw"-only MUST tests. Deterministic. Edge cases (empty/no-overlap/null/negative/extreme) covered.

### Findings / defects
None. Build succeeds (0 warnings). `ImmuneAnalyzer_ImmuneInfiltration_Tests`: 33/33 pass. Full `Seqeron.Genomics.Tests`: **4486 passed, 0 failed** (matches baseline).

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — every formula/threshold/coefficient matches authoritative primary sources exactly; the only notes are *honestly-declared* simplifications (5-marker matrix vs LM22; NNLS vs ν-SVR; Affymetrix-calibrated purity on a non-cohort-scaled ssGSEA score → absolute purity not clinically meaningful, but never advertised as such).
- **Stage B: PASS** — code faithfully realises the validated formulas; all cross-checks recomputed and matched; tests are real and deterministic.
- **End state: CLEAN** — no defect found; nothing changed. Project builds green; unit tests 33/33; full suite 4486/4486.
