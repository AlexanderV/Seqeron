# Validation Report: ONCO-IMMUNE-001 — Tumour Immune-Infiltration Estimation

- **Validated:** 2026-06-25 (fresh re-validation; supersedes 2026-06-24)   **Area:** Oncology / Tumour Immunology
- **Canonical surface (this unit):**
  - `ImmuneAnalyzer.EstimateInfiltration(...)` — ssGSEA immune/stromal enrichment scoring + relative tumour-purity (Yoshihara/Barbie).
  - `ImmuneAnalyzer.EstimateTumorPurity(double estimateScore)` — absolute Yoshihara (2013) ESTIMATE-purity closed form `cos(a + b·ESTIMATEScore)` with NaN out-of-domain handling.
  - `ImmuneAnalyzer.DeconvoluteImmuneCells(...)` — NNLS/LLSR cell-type deconvolution (Lawson–Hanson / Abbas 2009) on the bundled representative signature matrix.
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End state:** ✅ CLEAN (no defect; no code/test changed this session — re-validation only)

**Scope note.** The CIBERSORT ν-SVR deconvolution (`DeconvoluteImmuneCellsNuSvr`), the LM22-format loader (`LoadSignatureMatrix`), and the bundled ABIS matrix (`LoadBundledAbisSignatureMatrix`) live in the **same source file and test fixture** but constitute the separate, already-CLEAN unit **IMMUNE-NUSVR-001**. They are referenced here but **not re-litigated**; LM22 being caller-supplied (Stanford no-redistribution licence) is the documented, acceptable boundary of that unit. This report validates only ONCO-IMMUNE-001's own surface above.

---

## Stage A — Description

### Sources opened THIS session & what they confirm
- **Yoshihara et al. (2013), Nat Commun 4:2612 (ESTIMATE)** — via two independent reference reimplementations retrieved this session:
  - `hacksig::hack_estimate` CRAN refman — *"Purity = cos(0.6049872018 + 0.0001467884 * ESTIMATE)"* and *"The ESTIMATE score is defined as the combination (i.e. sum) of immune and stroma scores."*
  - `tidyestimate` reference R (`estimate_score.R`) — verbatim `cos(0.6049872018 + 0.0001467884 * estimate)` followed by `ifelse(purity < 0, NA, purity)`.
  Both coefficients match the code constants `EstimatePurityCoefficientA = 0.6049872018` / `EstimatePurityCoefficientB = 0.0001467884` to the last digit. ESTIMATE score = immune + stromal. Calibration is Affymetrix/ABSOLUTE-derived (scope note below).
- **Barbie et al. (2009), Nature 462:108–112 + Hänzelmann et al. (2013) GSVA + GenePattern ssGSEAProjection docs** — ssGSEA definition confirmed: genes are replaced by their **ranks** L={r₁…r_N} ordered N→1; the enrichment score is the **sum (integration)** of the difference between the weighted ECDF of in-set genes `P_wG` and the ECDF of the remaining genes `P_NG`; the hit weight is `|rank|^τ` with **τ = 0.25** for `method="ssgsea"`.
  - **Decisive detail** confirmed against the GSVA source (`R/ssgsea.R`): `R <- apply(X, 2, function(x) as.integer(rank(x)))` then `Ra <- abs(R)^alpha` (alpha=0.25). The weight is therefore on the **integer rank-order statistic** (1…N), **not** the expression value and **not** the max-deviation of classic GSEA. The code's `rank = n − i` (highest expression → rank N) with `Math.Pow(rank, 0.25)` and the integral-of-running-sum is exactly this.

### Formula check
| Claim | Source (this session) | Code | Match |
|---|---|---|---|
| Purity = cos(a + b·ESTIMATE), a=0.6049872018, b=0.0001467884 | Yoshihara 2013 / hacksig / tidyestimate | `EstimateTumorPurity` + `ComputeTumorPurity` | ✅ exact |
| ESTIMATE score = Immune + Stromal | hacksig refman | `estimateScore = immune + stromal` | ✅ |
| Out-of-domain (cos < 0) → NA | tidyestimate `ifelse(purity<0,NA,purity)` | `EstimateTumorPurity`: `purity < 0 ? NaN : purity` | ✅ exact |
| ssGSEA τ=0.25; hit weight = `rank^τ` (integer rank N…1); score = Σ(P_hit − P_miss) | Barbie 2009 / GSVA `ssgsea.R` | `ComputeSsGseaScore` | ✅ |
| NNLS min‖m−Sf‖² s.t. f≥0, then Σf=1 | Lawson–Hanson / Abbas 2009 | `SolveNnls` + normalization | ✅ |

### Edge-case semantics (sourced/defined)
- **`EstimateInfiltration`** empty profile → Immune=Stromal=ESTIMATE=0, TumorPurity = cos(a) (clamped path). ✅
- No overlapping genes → scores 0 (guards `nHits==0 || nMiss==0`, `totalHitWeight==0`). ✅
- **`EstimateTumorPurity`** mirrors the R reference: returns `NaN` when cos<0 (arg past π/2, score ≳ 6579.6), **not** a clamped 0 — a clamped 0 would falsely read as a valid all-stroma estimate. The default-path `InfiltrationResult.TumorPurity` instead clamps to [0,1]; this is an intentional difference (the un-normalised single-sample ssGSEA integral never reaches the NaN domain, so the relative field stays in range). ✅
- `DeconvoluteImmuneCells`: no overlapping genes → all fractions 0, OverlappingGenes=0; pure type → exact identity. ✅
- Null → `ArgumentNullException` (both methods). ✅

### Independent cross-check (hand computation THIS session — exact numbers)
Yoshihara purity `cos(0.6049872018 + 0.0001467884·s)`, NaN if negative:

| ESTIMATE score s | cos argument | purity | test |
|---|---|---|---|
| 0 | 0.604987201800 | **0.8225093766958238** | E1 ✅ |
| 1000 | 0.751775601800 | **0.7304773970805112** | E2 ✅ |
| 2500 | 0.971958201800 | 0.5636831556684971 | E7 (identity) ✅ |
| 3000 | 1.045352401800 | **0.5015970942006772** | E3 ✅ |
| 6000 | 1.485717601800 | **0.0849761233112934** | E5 (below cutoff) ✅ |
| 6600 | 1.573790641800 | cos = −0.00299 → **NaN** | E5 (above cutoff) ✅ |
| 7000 | 1.632506001800 | cos = −0.06167 → **NaN** | E4 ✅ |
| −2000 | 0.311410401800 | 0.9519023675717762 | E6 (monotone) ✅ |

Negative-cosine cutoff (arg = π/2): **s = (π/2 − a)/b ≈ 6579.601** — consistent with E5's "~6579.6" boundary.

ssGSEA controlled small input — profile {A=100, B=1, C=0.5}, set {A,C}, τ=0.25, ranked desc A(rank3)/B(rank2)/C(rank1), nHits=2, nMiss=1, missStep=1, TW=3^¼+1:
- hit A → RS=3^¼/TW, integral+=RS; miss B → RS−=1; hit C → RS+=1/TW=0; integral = **(3^¼−1)/(3^¼+1) = 0.13646973766…** (M14a — matches the test's runtime closed-form expression; *note: the prior report's prose figure "≈0.136548" was a rounding slip — the true value is 0.136470*).
- Single hit at top rank (3 genes): RS 1 → 0.5 → 0, integral **+1.5** (M14b — distinct from the 1.0 a max-deviation GSEA would give).
- Single hit at bottom rank: RS −0.5 → −1.0 → 0, integral **−1.5** (M14c).

### Findings / divergences (Stage A)
1. **Honest-scope note (not a defect):** the `DeconvoluteImmuneCells` default signature is a representative 5-marker × 22-cell-type matrix, *not* the full ~547-gene LM22; this default deconvolution is NNLS, *not* ν-SVR (ν-SVR + LM22-loader = the separate IMMUNE-NUSVR-001). The relative `InfiltrationResult.TumorPurity` applies the Affymetrix-calibrated cosine to an un-normalised single-sample ssGSEA integral, so that field's absolute number is not clinically meaningful — explicitly declared in the XML docs; tests assert only formula identities + invariants, never clinical accuracy. `EstimateTumorPurity` is the opt-in path that applies the cosine to a caller-supplied, ESTIMATE-scale Affymetrix score (the clinically-calibrated route). → PASS-WITH-NOTES, correctly scoped.
2. **ssGSEA weighting re-confirmed correct (decisive Stage-A point):** `rank^0.25` (integer rank order N…1), integral form. Re-verified against the GSVA `ssgsea.R` source this session — a genuine source-grounded match, not a code echo.
3. **`DefaultImmuneSignatureGenes` / `DefaultStromalSignatureGenes` are now the full 141-gene ESTIMATE signatures** (extracted from ESTIMATE R `SI_geneset.gmt`), correcting the earlier "5-marker default for infiltration" wording — the 5-marker matrix is only the *deconvolution* default. Not a defect; the gene-list content is a documented data table, and tests exercise the algorithm via controlled custom sets (same code path) plus the default-signature ordering/identity tests.

---

## Stage B — Implementation

- **Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs`
  - `EstimateInfiltration` (402–432), `EstimateTumorPurity` (473–481), `ComputeTumorPurity` (1256–1260), `ComputeSsGseaScore` (911–972)
  - `DeconvoluteImmuneCells` (512–599), `SolveNnls` (1267–1362) + passive-set LSQ (1395–1440) + Gaussian elimination (1445–1513), Pearson/RMSE helpers (1518–1560).

### Formula realised correctly?
- **ssGSEA:** ranks descending; `totalHitWeight = Σ (n−i)^0.25` over hits; walks accumulating `(n−i)^0.25/totalHitWeight` on hits and `−1/nMiss` on misses, summing the running difference into `integral`. Exactly Σ(P_hit − P_miss) with integer-rank-order weighting. ✅
- **`EstimateTumorPurity`:** `cos(a + b·estimateScore)`, then `purity < 0 ? NaN : purity` — verbatim Yoshihara + tidyestimate NA handling. ✅
- **`InfiltrationResult.TumorPurity`** (`ComputeTumorPurity`): same cosine, `Math.Clamp(_,0,1)` (relative-path design). ✅
- **NNLS:** Lawson–Hanson active-set; gradient `w = Aᵀ(b − Ax)`; passive-set normal equations via Gaussian elimination w/ partial pivoting; feasibility back-off via α-ratio; non-negativity cleanup; Σf=1 normalization guarded by `sum > 0`. ✅

### Cross-verification table recomputed vs code (fixture executed this session)
| ID | Expected (sourced/hand-computed) | Code | Match |
|----|----|----|----|
| E1 | purity(0) = 0.8225093766958238 | ✅ | ✅ |
| E2 | purity(1000) = 0.7304773970805112 | ✅ | ✅ |
| E3 | purity(3000) = 0.5015970942006772 | ✅ | ✅ |
| E4 | purity(7000) = NaN | ✅ | ✅ |
| E5 | purity(6000)=0.0849761233112934; purity(6600)=NaN | ✅ | ✅ |
| E6 | strictly decreasing over −2000…6000 | ✅ | ✅ |
| E7 | purity(2500) = cos(a+b·2500) identity | ✅ | ✅ |
| M14a | (3^¼−1)/(3^¼+1) ≈ 0.13646974 | ✅ | ✅ |
| M14b/c | +1.5 / −1.5 integral | ✅ | ✅ |
| M1 | empty → scores 0, purity=cos(a) | ✅ | ✅ |
| M10/INV-4 | ESTIMATE = immune + stromal | ✅ | ✅ |
| M9/INV-3 | purity = clamp(cos(a+b·score)) ∈ [0,1] | ✅ | ✅ |
| M5/M6/S1 | pure/50:50/75:25 NNLS exact fractions | ✅ | ✅ |
| INV-1/2 | f≥0, Σf=1 | ✅ | ✅ |

### Variant/delegate consistency
`EstimateInfiltration`, `EstimateTumorPurity`, `DeconvoluteImmuneCells` are independent static methods (no `*Fast`/instance variants). `EstimateTumorPurity` and the `InfiltrationResult.TumorPurity` field share the same cosine; E7 confirms the closed-form identity. Test-file constants mirror the source constants exactly.

### Numerical robustness
ssGSEA → 0 on `nHits==0 || nMiss==0` and `totalHitWeight==0`; purity NaN guard (no clamp masking); NNLS normalization guarded by `sum > 0`; Pearson denominator `< 1e-15 → 0`; Gaussian pivot `< 1e-15` skipped; RMSE `n==0 → 0`. Negative (log-transformed) expression handled (C2/C2b). ✅

### Test quality audit
The ONCO-IMMUNE-001 tests (M1, M3, M5, M6, M9–M14, S1–S4, C1–C3, INV-1…5, E1–E7) assert exact sourced/hand-computed values (1e-10 for identities, 1e-6 for computed), are deterministic, and cover every Stage-A path: purity at several scores, out-of-domain → NaN, the NaN-vs-clamp distinction, immune-vs-stromal scoring, ssGSEA rank-vs-expression and integral-vs-max-deviation discriminators (M14a/b/c), empty/no-overlap/null/negative/extreme. No tautological "does not throw"-only MUST tests on the canonical surface; E1–E7 trace to the Yoshihara formula and the tidyestimate NA rule, not to code echoes. Coverage of the canonical surface is complete. (The fixture's NSVR-*/ABIS-* tests belong to IMMUNE-NUSVR-001.)

Fixture run this session: **65/65 pass** (ONCO-IMMUNE-001 + the co-located IMMUNE-NUSVR-001 tests), 0 failed. No code/test changed → full unfiltered suite not re-run (per protocol, doc-only sessions keep the prior green baseline).

### Findings / defects
None on the ONCO-IMMUNE-001 canonical surface.

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — purity coefficients, NA out-of-domain rule, ESTIMATE=immune+stromal, and the ssGSEA integer-rank `rank^0.25` integral all match authoritative reference implementations retrieved this session, verbatim. Notes are honest, declared scope simplifications (representative deconvolution matrix vs LM22; relative single-sample purity vs absolute `EstimateTumorPurity`).
- **Stage B: PASS** — code faithfully realises the validated formulas; every cross-check recomputed by hand and via the fixture matches.
- **End state: ✅ CLEAN** — no defect; nothing changed. ONCO-IMMUNE-001 tests green within the 65/65 fixture. LM22/exact-CIBERSORT remain the separate IMMUNE-NUSVR-001 (caller-supplied LM22 = accepted boundary).

---

### Historical appendix (prior reports, retained as evidence)
The pre-2026-06-25 report (2026-06-24 and earlier) covered the same unit before the campaign split out IMMUNE-NUSVR-001 and added `EstimateTumorPurity`. Key historical facts retained: the original ESTIMATE-purity/ssGSEA/NNLS implementation was validated CLEAN under `cb113ce`/`a3b9e83b`; the opt-in `EstimateTumorPurity` and the ν-SVR + ABIS additions were limitation fixes (now IMMUNE-NUSVR-001). The 2026-06-24 prose figure for M14a ("≈0.136548") was a rounding slip; the source-exact value is **(3^¼−1)/(3^¼+1) = 0.13646973766…**, which is what the test asserts at runtime — no code defect.

## Runtime enforcement (LimitationPolicy)

This unit's guarded branch — ABIS/caller-matrix immune deconvolution and the ESTIMATE→ABSOLUTE purity transform (no CIBERSORT-LM22 parity) — has **minimum access mode `Moderate`** (`Seqeron.Genomics.Core.LimitationCatalog`). Under the default `LimitationPolicy.DefaultMode = Moderate` it is **allowed** (this guarded branch throws only under `Strict`); see [LIMITATIONS.md](../LIMITATIONS.md) › Runtime enforcement. Additive policy layer; the validated contract and `✅ CLEAN` verdict are unchanged.
