# Validation Report: ONCO-IMMUNE-001 — Immune Infiltration Estimation

- **Validated:** 2026-06-24   **Area:** Oncology / Tumor Immunology
- **Canonical method(s):** `ImmuneAnalyzer.EstimateInfiltration(...)` (ssGSEA enrichment + ESTIMATE tumor purity), `ImmuneAnalyzer.DeconvoluteImmuneCells(...)` (NNLS cell-type deconvolution)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End state:** CLEAN (no defect; code & tests unchanged since `a3b9e83b`)

This unit computes two tumor-immunology metrics:
1. **ESTIMATE-style immune/stromal infiltration** via single-sample GSEA (ssGSEA, integral/area form, τ=0.25) plus a tumor-purity estimate `cos(a + b·ESTIMATE)`.
2. **Immune cell-type deconvolution** via Non-Negative Least Squares (NNLS, Lawson–Hanson active set), explicitly *not* CIBERSORT's ν-SVR.

This is a re-validation in a fresh context. `git log` on `src/.../ImmuneAnalyzer.cs` shows a single commit (`a3b9e83b feat(oncology): ONCO-IMMUNE-001 ...`); the implementation and the 33-test file are byte-identical to the state validated under the prior report (`cb113ce`).

---

## Stage A — Description

### Sources opened & what they confirm
- **Yoshihara et al. (2013), Nat Commun 4:2612 (ESTIMATE)** — tumor-purity formula `Tumor_purity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`. Both coefficients match the code constants `EstimatePurityCoefficientA/B` to the last digit. ESTIMATE score is defined as immune + stromal. The calibration is Affymetrix-derived (scope note, below).
- **Barbie et al. (2009), Nature 462:108–112 + Hänzelmann et al. (2013), GSVA, BMC Bioinformatics 14:7** — opened the GSVA Bioconductor vignette and the GSVA/ssGSEA method description. Authoritative statement of the ssGSEA procedure: *"the genes are replaced by their **ranks** according to their absolute expression L={r1,…,rN}. The list is ordered from the highest rank N to the lowest 1. An enrichment score ES(G,S) is obtained by a sum (integration) of the difference between a weighted ECDF of the genes in the signature P_wG and the ECDF of the remaining genes P_NG."* The hit weight is **|r_j|^τ with τ=0.25 for `method="ssgsea"`**, where **r_j is the integer rank order statistic** (N…1), not the expression value. This is the decisive Stage-A point — see Finding 2.
- **Newman et al. (2015), Nat Methods 12:453 (CIBERSORT) + Abbas et al. (2009), PLoS One 4:e6098 + Lawson & Hanson (1995)** — deconvolution model `m = S·f`, solved `min ‖m − S·f‖² s.t. f ≥ 0`, then normalized so `Σf = 1`. NNLS (Lawson–Hanson active-set) is the LLSR/NNLS baseline; CIBERSORT itself uses ν-SVR. 22 LM22-style phenotypes correctly enumerated.

### Formula check
| Claim | Source | Code | Match |
|---|---|---|---|
| Purity = cos(a + b·ESTIMATE), a=0.6049872018, b=0.0001467884 | Yoshihara 2013 | `ComputeTumorPurity` | ✅ exact |
| ESTIMATE score = Immune + Stromal | Yoshihara 2013 | `estimateScore = immune + stromal` | ✅ |
| ssGSEA τ=0.25; hit weight = `rank^τ`, rank=N−i (highest expr → N); score = integral of (P_hit − P_miss) | Barbie 2009 / GSVA | `ComputeSsGseaScore` | ✅ |
| NNLS min‖m−Sf‖² s.t. f≥0, then Σf=1 | Lawson–Hanson / Abbas 2009 | `SolveNnls` + normalization | ✅ |

### Edge-case semantics (all sourced/defined)
- Empty profile → ImmuneScore=0, StromalScore=0, TumorPurity=cos(a). ✅
- No overlapping genes → scores 0 (guards `nHits==0 || nMiss==0`, `totalHitWeight==0`). ✅
- No overlapping deconvolution genes → all fractions 0, OverlappingGenes=0. ✅
- Purity can exceed [0,1] at extreme scores → `Math.Clamp(_,0,1)`. ✅
- Null → `ArgumentNullException`. ✅

### Independent cross-check (hand computation, this session)
Profile {A=100, B=1, C=0.5}, gene set {A,C}. Ranked desc: A(i=0,rank=3), B(i=1,rank=2), C(i=2,rank=1). N=3, nHits=2, nMiss=1, missStep=1, TW=3^¼+1^¼=3^¼+1.
Walk: hit A → RS=3^¼/TW, integral=3^¼/TW; miss B → RS=3^¼/TW−1, integral+= ; hit C → RS=(3^¼+1)/TW−1=0, integral+=0.
**Integral = 2·3^¼/TW − 1 = (3^¼ − 1)/(3^¼ + 1) ≈ 0.136548.** Matches test M14a exactly and is distinct from the ≈0.57992 an expression-value-weighted variant would yield.
Single hit at top rank: hit(+1)→RS 1, miss(−0.5)→0.5, miss(−0.5)→0 ⇒ integral **1.5** (M14b). Single hit at bottom: miss(−0.5)→−0.5, miss(−0.5)→−1.0, hit(+1)→0 ⇒ integral **−1.5** (M14c). All three recomputed by hand and confirmed.

### Findings / divergences (Stage A)
1. **Honest-scope note (not a defect):** default signature matrix is a simplified 5-marker × 22-cell-type matrix, *not* the full ~547-gene LM22; deconvolution is NNLS, *not* ν-SVR; ESTIMATE purity coefficients are Affymetrix-calibrated, applied to a single-sample ssGSEA integral on a different numeric scale than the cohort-scaled ESTIMATE score — so the absolute TumorPurity number is not clinically meaningful. This is explicitly declared in the XML docs and the spec's Assumption Register; no clinical/diagnostic advertising is present. Tests assert only formula identities and mathematical invariants, never clinical accuracy. → PASS-WITH-NOTES, correctly scoped.
2. **ssGSEA weighting confirmed correct (the key Stage-A question).** Real GSVA/Barbie ssGSEA weights hits by `|rank|^τ` where the rank is the integer rank-order statistic (N…1), τ=0.25 — *not* by the expression value^τ, and *not* by max-deviation (which is classic GSEA). The code's `rank = n − i` with `Math.Pow(rank, 0.25)` and the integral-of-running-sum construction is exactly this. The earlier TestSpec correction (from `|expr|^τ` to `rank^τ`) was therefore the right call. This is a genuine, source-grounded match, not a coincidence.
3. **Minor divergence (documented, not a defect):** GSVA's ssGSEA has a final normalization step (divide the integral by the range of values across the sample cohort). This implementation is strictly single-sample and omits cohort normalization (there is no cohort), so its enrichment scores are the un-normalized integral. This is consistent with the single-sample, formula-identity scope and is why the absolute purity is not on the ESTIMATE numeric scale (Finding 1). Tests assert the un-normalized integral, which is the correct value for the per-sample computation as implemented.

---

## Stage B — Implementation

- **Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs`
  - `EstimateInfiltration` (348–378), `ComputeSsGseaScore` (513–574), `ComputeTumorPurity` (580–584)
  - `DeconvoluteImmuneCells` (409–496), `SolveNnls` (591–686), passive-set LSQ (719–764) + Gaussian elimination (769–837), Pearson/RMSE helpers (842–884).

### Formula realised correctly?
- **ssGSEA:** ranks descending; `totalHitWeight = Σ (n−i)^τ` over hits; walks accumulating `(n−i)^τ/totalHitWeight` on hits and `−1/nMiss` on misses, summing the running difference into `integral`. Exactly the integral of (P_hit − P_miss) with rank-order weighting. ✅
- **Purity:** `cos(a + b·estimateScore)` clamped to [0,1]. ✅
- **NNLS:** Lawson–Hanson active-set; gradient `w = Aᵀ(b − Ax)`; passive-set normal-equations via Gaussian elimination w/ partial pivoting; feasibility back-off via α-ratio; final non-negativity cleanup; then Σf=1 normalization guarded by `sum > 0`. ✅

### Cross-verification table recomputed vs code (tests executed this session)
| ID | Expected (sourced) | Code result | Match |
|----|--------------------|-------------|-------|
| M1 | scores 0, purity=cos(a) | ✅ | ✅ |
| M5 | pure CD8 → f_CD8=1.0, corr=1, RMSE=0 | ✅ | ✅ |
| M6 | B+CD8 50:50 → 0.5/0.5 | ✅ | ✅ |
| S1 | 75:25 → 0.75/0.25 | ✅ | ✅ |
| M9/INV-3 | purity=cos(a+b·score) ∈ [0,1] | ✅ | ✅ |
| M10/INV-4 | ESTIMATE = immune+stromal | ✅ | ✅ |
| M14a | (3^¼−1)/(3^¼+1) ≈ 0.136548 | ✅ | ✅ |
| M14b/c | +1.5 / −1.5 integral | ✅ | ✅ |
| INV-1/2 | f≥0, Σf=1 (6 cell types) | ✅ | ✅ |
| C1/C3 | 22 cell types, 5 genes each | ✅ | ✅ |

### Variant/delegate consistency
Two canonical static methods; no `*Fast`/instance variants to reconcile. Test-file constants mirror the source constants exactly.

### Numerical robustness
Div guards: ssGSEA → 0 on `nHits==0 || nMiss==0` and `totalHitWeight==0`; normalization guarded by `sum > 0`; Pearson denominator `< 1e-15 → 0`; Gaussian pivot `< 1e-15` skipped; RMSE `n==0 → 0`. Negative (log-transformed) expression handled (C2/C2b). No div-by-zero on stated ranges. ✅

### Test quality audit
33 tests, all asserting exact sourced values (1e-10 for identities, 1e-6 for computed). M14a/b/c are genuine discriminating references (rank-vs-expression weighting; integral-vs-max-deviation) — they would catch the two most likely ssGSEA bugs. No tautological "does not throw"-only MUST tests. Deterministic. Edge cases (empty/no-overlap/null/negative/extreme) covered.

### Findings / defects
None. Build succeeds (0 warnings). `ImmuneAnalyzer_ImmuneInfiltration_Tests`: **33/33 pass**. No code changed this session, so the full suite was not re-run (baseline unchanged).

---

## Update 2026-06-24 — limitation fix: opt-in absolute purity (`EstimateTumorPurity`)

The "relative, not clinically-absolute purity" limitation (Finding 1 / Finding 3) is now addressed by an **opt-in** addition; the default 5-marker/ssGSEA `EstimateInfiltration` path is unchanged.

- **New public method:** `ImmuneAnalyzer.EstimateTumorPurity(double estimateScore)`. Applies the verbatim Yoshihara (2013) closed-form transform `purity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)` to a **caller-supplied, Affymetrix/ABSOLUTE-calibrated ESTIMATE score** (the original ESTIMATE numeric scale), producing an absolute purity rather than the relative single-sample value of the `InfiltrationResult.TumorPurity` field.
- **Domain handling (reference-grounded):** mirrors the ESTIMATE/`tidyestimate` reference implementation `purity = ifelse(purity < 0, NA, purity)` — when the cosine evaluates negative (cos argument past π/2, ESTIMATE score ≳ 6579.6) the result is `double.NaN` (out of the calibrated domain), **not** a clamped 0. Calibration is Affymetrix-only (nonlinear least squares vs ABSOLUTE on TCGA), so it is invalid for RNA-seq-derived scores — documented in the XML doc, algorithm doc, and Evidence.
- **Sources retrieved this session (URLs):** `https://search.r-project.org/CRAN/refmans/hacksig/html/hack_estimate.html` and `https://www.aging-us.com/article/203714/text` (formula + coefficients verbatim, two independent sources); `https://raw.githubusercontent.com/KaiAragaki/tidyestimate/main/R/estimate_score.R` (reference R: `cos(0.6049872018 + 0.0001467884 * estimate)` + `ifelse(purity < 0, NA, purity)` + `is_affymetrix` gate); CIBERSORT download page search (LM22 gated behind academic registration).
- **Tests added (E1–E7):** exact hand-computed cosine values — purity(0)=0.8225093766958238, purity(1000)=0.7304773970805112, purity(3000)=0.5015970942006772, purity(6000)=0.0849761233112934 — `Within(1e-10)`; NaN for out-of-domain score 7000/6600; strict monotone-decreasing across −2000…6000; closed-form identity at 2500. Fixture now 40/40 green; full suite green.
- **Honest residual (per STOP RULE):** the CIBERSORT **LM22** 547-gene × 22-cell-type signature matrix and **ν-SVR** solver (Newman et al., 2015) were **not** implemented — LM22 is gated behind academic registration on the CIBERSORT website (not cleanly retrievable as plaintext this session) and a faithful ν-SVR depends on the trained matrix. Not fabricated; the default deconvolution remains the NNLS/LLSR baseline on a representative 5-marker matrix.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — every formula/coefficient matches authoritative primary sources exactly. The decisive ssGSEA rank-order weighting (`rank^τ`, τ=0.25, integral form) was independently re-confirmed against the GSVA/Barbie description this session. Notes are honestly-declared simplifications (5-marker matrix vs LM22; NNLS vs ν-SVR; un-normalized single-sample integral on a non-cohort scale ⇒ absolute purity not clinically meaningful, never advertised as such).
- **Stage B: PASS** — code faithfully realises the validated formulas; all cross-checks recomputed by hand and via tests, all match.
- **End state: CLEAN** — no defect found; nothing changed. Unit tests 33/33 green.
