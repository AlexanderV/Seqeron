# Test Specification: ONCO-IMMUNE-001

**Test Unit ID:** ONCO-IMMUNE-001
**Area:** Oncology
**Algorithm:** Immune Infiltration Estimation (4 methods)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Newman et al. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods. | 1 (Peer-reviewed) | https://doi.org/10.1038/nmeth.3337 | 2026-03-06 |
| 2 | Yoshihara et al. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature Communications. | 1 (Peer-reviewed) | https://doi.org/10.1038/ncomms3612 | 2026-03-06 |
| 3 | Barbie et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. | 1 (Peer-reviewed) | https://doi.org/10.1038/nature08460 | 2026-03-06 |
| 4 | Subramanian et al. (2005). Gene set enrichment analysis. PNAS. | 1 (Peer-reviewed) | https://doi.org/10.1073/pnas.0506580102 | 2026-03-06 |
| 5 | Abbas et al. (2009). Deconvolution of blood microarray data identifies cellular activation patterns in systemic lupus erythematosus. PLoS One. | 1 (Peer-reviewed) | https://doi.org/10.1371/journal.pone.0006098 | 2026-03-06 |
| 6 | Lawson & Hanson (1995). Solving Least Squares Problems. SIAM Classics in Applied Mathematics. | 2 (Textbook) | ISBN 978-0-89871-356-5 | 2026-03-06 |
| 7 | Wikipedia — CIBERSORT | 4 (Wikipedia) | https://en.wikipedia.org/wiki/CIBERSORT | 2026-03-06 |
| 8 | Becht et al. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology. | 1 (Peer-reviewed) | https://doi.org/10.1186/s13059-016-1070-5 | 2026-03-06 |
| 9 | Hänzelmann et al. (2013). GSVA: gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics. | 1 (Peer-reviewed) | https://doi.org/10.1186/1471-2105-14-7 | 2026-03-06 |
| 10 | Schölkopf et al. (2000). New support vector algorithms. Neural Computation 12(5):1207-1245 (ν-SVR dual). | 1 (Peer-reviewed) | https://doi.org/10.1162/089976600300015565 ; eqs 60–62: https://alex.smola.org/papers/2003/SmoSch03b.pdf | 2026-06-25 |
| 11 | Chen et al. (2018). Profiling Tumor Infiltrating Immune Cells with CIBERSORT. Methods Mol Biol 1711:243-259. | 1 (Peer-reviewed protocol) | https://pmc.ncbi.nlm.nih.gov/articles/PMC5895181/ | 2026-06-25 |
| 12 | CIBERSORT licence (Stanford). No-redistribution, non-commercial terms; LM22 registration gate. | 1 (governing licence) | https://gist.github.com/dhimmel/58dcd9b512e669f20a65ddf73997b733 ; https://cibersort.stanford.edu | 2026-06-25 |
| 13 | scikit-learn 1.6.1 `NuSVR` (libsvm) — cross-implementation ν-SVR reference. | 3 (reference implementation) | https://scikit-learn.org/stable/modules/generated/sklearn.svm.NuSVR.html | 2026-06-25 |

### 1.2 Key Evidence Points

1. **NNLS deconvolution**: The NNLS/LLSR approach solves min ||m - S·f||² s.t. f ≥ 0 to estimate proportions of 22 immune cell types from bulk RNA expression — Abbas et al. (2009), Lawson & Hanson (1995). Newman et al. (2015) benchmarked this as one of several deconvolution methods alongside their ν-SVR-based CIBERSORT.
2. Cell fractions must satisfy: f ≥ 0 and Σf = 1 after normalization — Newman et al. (2015)
3. Deconvolution model: m = S × f, solved as min ||m - S·f||² s.t. f ≥ 0 — Lawson & Hanson (1995), Abbas et al. (2009)
4. ESTIMATE computes immune/stromal scores using ssGSEA enrichment, then derives tumor purity via cos(a + b × estimateScore) — Yoshihara et al. (2013)
5. ssGSEA enrichment score: integral (sum) of the weighted running sum across all ranked positions, with hit weighting by rank^τ (τ=0.25, rank = N−i for gene at descending-sorted position i) — Barbie et al. (2009), Hänzelmann et al. (2013). The GSVA ssGSEA function weights hits by rank (integer position), not by expression value.
6. Tumor purity coefficients: a = 0.6049872018, b = 0.0001467884 — Yoshihara et al. (2013)
7. **ν-SVR deconvolution (CIBERSORT)**: linear-kernel ν-SVR of `m` on the columns of `B`; dual: maximise `−½Σ(α_i−α_i*)(α_j−α_j*)⟨x_i,x_j⟩ + Σy_i(α_i−α_i*)` s.t. `Σ(α_i−α_i*)=0`, `Σ(α_i+α_i*) ≤ Cνℓ`, `α_i,α_i*∈[0,C]`; `w = Σ(α_i−α_i*)x_i` — Schölkopf et al. (2000), eqs 60–62.
8. **ν sweep + selection**: CIBERSORT sweeps ν ∈ {0.25, 0.5, 0.75} and keeps the ν with the lowest RMSE between `m` and `B·f`; then zero-clips negative weights and normalises to sum 1 — Chen et al. (2018), Newman et al. (2015).
9. **LM22 dimensions + licence**: LM22 = 547 genes × 22 cell types; Stanford licence forbids redistribution ("RECIPIENT shall not distribute the Program …"), so LM22 is caller-supplied — Newman et al. (2015); CIBERSORT licence (source 12).

### 1.3 Documented Corner Cases

1. **Empty expression profile**: No genes → no enrichment/deconvolution possible — Mathematical definition
2. **No overlapping genes**: If signature genes don't overlap with input, scores should be zero — Mathematical definition
3. **Collinear cell types**: Closely related cell types (e.g., resting/activated NK) are harder to distinguish via NNLS — Newman et al. (2015), Abbas et al. (2009)
4. **ESTIMATE score extremes**: Purity formula may exceed [0,1] at extreme scores → clamping required — Yoshihara et al. (2013)

### 1.4 Known Failure Modes / Pitfalls

1. **Low signal-to-noise**: High noise degrades deconvolution accuracy — Newman et al. (2015)
2. **Too few overlapping genes**: If most signature genes are missing, deconvolution becomes underdetermined — Newman et al. (2015)
3. **Non-hematopoietic content**: NNLS deconvolution only estimates hematopoietic fractions; residual represents non-immune cells — Newman et al. (2015), Abbas et al. (2009)

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `EstimateInfiltration` | `ImmuneAnalyzer` | **Canonical** | ssGSEA-based immune/stromal scoring (integral form, τ=0.25) + tumor purity per Yoshihara et al. (2013), Barbie et al. (2009), Hänzelmann et al. (2013) |
| `EstimateTumorPurity` | `ImmuneAnalyzer` | **Canonical** | Opt-in absolute-purity transform `cos(0.6049872018 + 0.0001467884 × score)` with NaN-on-negative domain handling per Yoshihara et al. (2013) + ESTIMATE/tidyestimate reference implementation |
| `DeconvoluteImmuneCells` | `ImmuneAnalyzer` | **Canonical** | NNLS-based immune cell type deconvolution per Lawson & Hanson (1995), Abbas et al. (2009) |
| `DeconvoluteImmuneCellsNuSvr` | `ImmuneAnalyzer` | **Canonical** | CIBERSORT-style linear-kernel ν-SVR deconvolution (ν sweep, lowest-RMSE, zero-clip, sum-to-1) per Schölkopf et al. (2000), Newman et al. (2015), Chen et al. (2018) |
| `LoadSignatureMatrix` | `ImmuneAnalyzer` | **Canonical** | LM22-format TSV loader (caller-supplied LM22 — not bundled, Stanford licence) per Newman et al. (2015) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | All cell fractions ≥ 0 | Yes | Lawson & Hanson (1995): NNLS non-negativity constraint |
| INV-2 | Sum of cell fractions = 1.0 (when any overlap exists) | Yes | Abbas et al. (2009): normalization step |
| INV-3 | Tumor purity ∈ [0, 1] | Yes | Yoshihara et al. (2013): clamped cosine formula |
| INV-4 | ESTIMATE score = Immune score + Stromal score | Yes | Yoshihara et al. (2013): definition |
| INV-5 | OverlappingGenes ≥ 0 and ≤ total signature genes | Yes | Mathematical property |
| INV-6 | `EstimateTumorPurity` is monotone-decreasing in the ESTIMATE score over the valid domain | Yes | Yoshihara et al. (2013): cos is decreasing on [0, π] |
| INV-7 | ν-SVR cell fractions ≥ 0 | Yes | Newman et al. (2015): negative weights clipped to 0 |
| INV-8 | ν-SVR cell fractions sum to 1 (when post-clip mass > 0) | Yes | Newman et al. (2015): normalisation step |
| INV-9 | `BestNu` ∈ {0.25, 0.5, 0.75} | Yes | Chen et al. (2018): CIBERSORT ν sweep |
| INV-10 | ν-SVR deconvolution is deterministic | Yes | No randomness in the solver (fixed coordinate-ascent order) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Empty expression → zero scores | `EstimateInfiltration({})` | ImmuneScore=0, StromalScore=0, TumorPurity=cos(a) | Mathematical definition |
| M2 | Empty expression → zero deconvolution | `DeconvoluteImmuneCells({})` | All fractions=0, Correlation=0, RMSE=0 | Mathematical definition |
| M3 | No overlapping genes → zero enrichment | Expression has no immune/stromal genes | ImmuneScore=0, StromalScore=0 | Mathematical definition |
| M4 | No overlapping genes → zero deconvolution | Expression genes not in signature | All fractions=0, OverlappingGenes=0 | Mathematical definition |
| M5 | Single cell type → 100% fraction | Expression = pure CD8 T cell signature | T_cells_CD8 = 1.0, others = 0 (exact, unique genes CD8A/CD8B) | Abbas et al. (2009): linear identity |
| M6 | Two cell types equal mix → 50:50 | Expression = average of B naive + CD8 T | B_naive = 0.5, CD8 = 0.5 (exact, unique genes force solution) | Abbas et al. (2009): linearity |
| M9 | Tumor purity in [0, 1] + formula | Any valid infiltration | 0 ≤ TumorPurity ≤ 1 AND purity = cos(a + b × estimateScore) | Yoshihara et al. (2013) |
| M10 | ESTIMATE score = immune + stromal | Any valid infiltration | EstimateScore = ImmuneScore + StromalScore | Yoshihara et al. (2013) |
| M11 | High immune > low immune ordering | High-immune vs low-immune profiles | HighProfile.ImmuneScore > LowProfile.ImmuneScore | ESTIMATE concept (Yoshihara et al., 2013) |
| M12 | Null expression throws ArgumentNullException | `EstimateInfiltration(null)` | ArgumentNullException | Robustness |
| M13 | Null expression throws for deconvolution | `DeconvoluteImmuneCells(null)` | ArgumentNullException | Robustness |
| M14 | ssGSEA exact value against hand-computed reference | Custom genes, rank-based integral | (a) score = (3^(1/4)−1)/(3^(1/4)+1); (b) top hit = 1.5; (c) bottom hit = −1.5 | Barbie et al. (2009), Hänzelmann et al. (2013): rank-based ssGSEA |
| E1 | `EstimateTumorPurity(0)` exact | Yoshihara transform at score 0 | cos(0.6049872018) = 0.8225093766958238 | Yoshihara et al. (2013); hacksig/tidyestimate |
| E2 | `EstimateTumorPurity(1000)` exact | mid-range score | cos(0.7517756018) = 0.7304773970805112 | Yoshihara et al. (2013) |
| E3 | `EstimateTumorPurity(3000)` exact | high score | cos(1.0453524018) = 0.5015970942006772 | Yoshihara et al. (2013) |
| E4 | Out-of-domain → NaN | `EstimateTumorPurity(7000)` (cos arg > π/2 → negative) | `double.NaN` | tidyestimate `estimate_score()`: `ifelse(purity<0, NA, purity)` |
| E5 | Negative-cosine cutoff boundary | score 6000 (defined) vs 6600 (NaN); cutoff at (π/2−a)/b ≈ 6579.6 | 6000 → 0.0849761233112934; 6600 → NaN | tidyestimate reference impl |
| E6 | Monotone decreasing (INV-6) | increasing scores −2000…6000 | purity strictly decreases | Yoshihara et al. (2013): cos decreasing on [0, π] |
| E7 | Closed-form identity | score 2500 | equals `cos(a + b × 2500)` | Yoshihara et al. (2013) |
| NSVR-M1 | ν-SVR planted-truth recovery | bulk = B·f, f={CD8:0.60, B_naive:0.30, Monocytes:0.10} on default matrix | recovers each planted fraction within 0.025; absent types ≈0 | Newman et al. (2015): linear mixture; Dataset 4 |
| NSVR-M2 | Match scikit-learn/libsvm `NuSVR` reference | disjoint 3×3 matrix; selected ν=0.75 | fractions = [TypeA 0.508497, TypeB 0.179491, TypeC 0.312012] within 2e-3; BestNu=0.75 | sklearn 1.6.1 NuSVR (Dataset 5); Schölkopf et al. (2000) |
| NSVR-M3 | ν-SVR fractions ≥ 0 and Σ = 1 | mixture of 3 cell types | all fractions ≥ 0; Σ = 1 (within 1e-9) | Newman et al. (2015): zero-clip + normalise (INV-7/8) |
| NSVR-M4 | Selected ν ∈ {0.25,0.5,0.75} | any valid mixture | `BestNu` ∈ CibersortNuValues | Chen et al. (2018) (INV-9) |
| NSVR-M5 | Reconstruction correlation near 1 | exact planted mixture | correlation > 0.95 | Newman et al. (2015): m=B·f is linear |
| NSVR-S1 | Determinism | same input twice | identical fractions and BestNu | No randomness (INV-10) |
| NSVR-S2 | No overlapping genes → zeros | genes absent from signature | all fractions=0, OverlappingGenes=0, BestNu=0 | Mathematical definition |
| NSVR-S3 | Null profile → throws | `DeconvoluteImmuneCellsNuSvr(null)` | ArgumentNullException | Robustness |
| NSVR-S4 | LM22-format loader parses TSV | header + gene rows | cell types from header; values parsed (incl. zeros) | Newman et al. (2015): LM22 TSV format |
| NSVR-S5 | Loaded matrix drives deconvolution | load disjoint matrix, deconvolve planted mix | recovered ordering A>B>C; Σ=1 | end-to-end loader→ν-SVR |
| NSVR-C1 | Loader rejects empty input | `LoadSignatureMatrix([])` | FormatException | Format validation |
| NSVR-C2 | Loader rejects header w/o cell types | header = "Gene symbol" only | FormatException | Format validation |
| NSVR-C3 | Loader rejects ragged row | row with wrong column count | FormatException | Format validation |
| NSVR-C4 | Loader rejects non-numeric value | value = "NOT_A_NUMBER" | FormatException | Format validation |
| NSVR-C5 | Loader rejects null lines | `LoadSignatureMatrix(null)` | ArgumentNullException | Robustness |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Unequal mixture → exact proportional fractions | 75% CD8 + 25% B naive | CD8 = 0.75, B_naive = 0.25 (exact, unique genes) | Validates linearity |
| S2 | Extra genes ignored in deconvolution | Profile has 100 extra genes + signature genes | Same result as signature-only profile | Robustness |
| S3 | Perfect reconstruction for pure profile | Pure cell type deconvolution | correlation = 1.0, RMSE = 0 | Quality metric |
| S4 | Custom gene sets accepted | User-provided immune/stromal gene lists | Uses custom sets, not defaults | API flexibility |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | 22 cell types all resolved | Full LM22-like signature | All 22 fractions computed | Scalability |
| C2 | Negative expression values handled | Log-transformed data with negatives | No exception, valid scores | Robustness |
| C3 | Default signature matrix has 22 types | Check DefaultSignatureMatrix.Count | Count = 22 | Completeness |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests found for `ImmuneAnalyzer` anywhere in the codebase.
- No file named `*Immune*Tests*` exists.
- The `Seqeron.Genomics.Oncology` project did not exist prior to this Test Unit processing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M6 (MUST) | ✅ Covered | M5: exact CD8=1.0 (unique genes proof); M6: exact 0.5:0.5 (unique genes proof) |
| M7 | 🔁 Merged | Subsumed by INV-1/2 parametric (Monocytes case) |
| M8 | 🔁 Merged | Subsumed by INV-1/2 parametric (NK_cells_activated case) |
| M9–M13 (MUST) | ✅ Covered | M9: strengthened with purity formula verification |
| M14 (MUST) | ✅ Covered | NEW: ssGSEA rank-based exact reference (3 sub-tests) |
| S1 | ✅ Covered | Strengthened: exact 0.75:0.25 (unique genes proof) |
| S2–S4 (SHOULD) | ✅ Covered | S3: strengthened to correlation=1.0, RMSE=0 |
| C1–C3 (COULD) | ✅ Covered | C2 has 2 sub-tests (EstimateInfiltration + DeconvoluteImmuneCells) |
| INV-1, INV-2 | ✅ Covered | 6 parameterized cases (includes former M7/M8 cell types) |
| INV-3, INV-4 | ✅ Covered | 4 parameterized cases; strengthened with purity formula check |
| INV-5 | ✅ Covered | Range check on overlapping genes |
| INV-6 | 🔁 Merged | Subsumed by M1 + M2 (empty profile tests) |
| NSVR-M1–M5 (MUST) | ❌ Missing | NEW this session: ν-SVR planted-truth, sklearn reference, INV-7/8/9, reconstruction |
| NSVR-S1–S5 (SHOULD) | ❌ Missing | NEW: determinism, no-overlap, null, loader parse, loader→deconvolution |
| NSVR-C1–C5 (loader validation) | ❌ Missing | NEW: empty / no-cell-types / ragged / non-numeric / null |
| INV-7, INV-8, INV-9, INV-10 | ❌ Missing | NEW: covered by NSVR-M3 (≥0, Σ=1), NSVR-M4 (ν set), NSVR-S1 (determinism) |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ImmuneAnalyzer_ImmuneInfiltration_Tests.cs` — all tests for ONCO-IMMUNE-001
- **Remove:** Nothing (no existing tests to remove)

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ImmuneAnalyzer_ImmuneInfiltration_Tests.cs` | Canonical | 55 |

### 5.5 Phase 7–8 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Created | ✅ Done |
| 2 | M2 | ❌ Missing | Created | ✅ Done |
| 3 | M3 | ❌ Missing | Created | ✅ Done |
| 4 | M4 | ❌ Missing | Created | ✅ Done |
| 5 | M5 | ⚠ Weak | Strengthened: exact CD8=1.0 (unique genes CD8A/CD8B proof) | ✅ Done |
| 6 | M6 | ⚠ Weak | Strengthened: exact B=0.5, CD8=0.5 (unique genes proof) | ✅ Done |
| 7 | M7 | 🔁 Duplicate | Removed: subsumed by INV-1/2 parametric (Monocytes case) | 🗑 Removed |
| 8 | M8 | 🔁 Duplicate | Removed: subsumed by INV-1/2 parametric (NK_activated case) | 🗑 Removed |
| 9 | M9 | ⚠ Weak | Strengthened: added purity = cos(a + b × estimateScore) formula check | ✅ Done |
| 10 | M10 | ✅ Covered | No change (exact identity check) | ✅ Done |
| 11 | M11 | ✅ Covered | No change (ordering test) | ✅ Done |
| 12 | M12 | ✅ Covered | No change (exception test) | ✅ Done |
| 13 | M13 | ✅ Covered | No change (exception test) | ✅ Done |
| 14 | M14 | ❌ Missing | NEW: 3 sub-tests verifying ssGSEA rank-based integral values | ✅ Done |
| 15 | S1 | ⚠ Weak | Strengthened: exact CD8=0.75, B=0.25 (unique genes proof) | ✅ Done |
| 16 | S2 | ✅ Covered | No change (exact comparison) | ✅ Done |
| 17 | S3 | ⚠ Weak | Strengthened: correlation=1.0, RMSE=0 (unique genes proof) | ✅ Done |
| 18 | S4 | ✅ Covered | No change (API contract test) | ✅ Done |
| 19 | C1 | ✅ Covered | No change | ✅ Done |
| 20 | C2 | ✅ Covered | No change (C2 + C2b) | ✅ Done |
| 21 | C3 | ✅ Covered | No change | ✅ Done |
| 22 | INV-1/2 | ✅ Covered | No change (6 parameterized, now also covers former M7/M8) | ✅ Done |
| 23 | INV-3/4 | ⚠ Weak | Strengthened: added purity formula verification | ✅ Done |
| 24 | INV-5 | ✅ Covered | No change | ✅ Done |
| 25 | (legacy INV-6 dup) | 🔁 Duplicate | Removed: subsumed by M1 + M2 | 🗑 Removed |
| 26 | E1 | ❌ Missing | NEW: exact `EstimateTumorPurity(0)` = 0.8225093766958238 | ✅ Done |
| 27 | E2 | ❌ Missing | NEW: exact `EstimateTumorPurity(1000)` = 0.7304773970805112 | ✅ Done |
| 28 | E3 | ❌ Missing | NEW: exact `EstimateTumorPurity(3000)` = 0.5015970942006772 | ✅ Done |
| 29 | E4 | ❌ Missing | NEW: out-of-domain (negative cosine) → NaN | ✅ Done |
| 30 | E5 | ❌ Missing | NEW: negative-cosine cutoff boundary (6000 defined, 6600 NaN) | ✅ Done |
| 31 | E6 | ❌ Missing | NEW: monotone-decreasing purity (INV-6) | ✅ Done |
| 32 | E7 | ❌ Missing | NEW: closed-form cosine identity at score 2500 | ✅ Done |
| 33 | NSVR-M1 | ❌ Missing | NEW: planted-truth recovery on 5-marker matrix (within 0.025) | ✅ Done |
| 34 | NSVR-M2 | ❌ Missing | NEW: scikit-learn/libsvm `NuSVR` cross-check, BestNu=0.75 | ✅ Done |
| 35 | NSVR-M3 | ❌ Missing | NEW: fractions ≥ 0 and Σ=1 (INV-7/8) | ✅ Done |
| 36 | NSVR-M4 | ❌ Missing | NEW: BestNu ∈ {0.25,0.5,0.75} (INV-9) | ✅ Done |
| 37 | NSVR-M5 | ❌ Missing | NEW: reconstruction correlation > 0.95 | ✅ Done |
| 38 | NSVR-S1 | ❌ Missing | NEW: determinism (INV-10) | ✅ Done |
| 39 | NSVR-S2 | ❌ Missing | NEW: no-overlap → zeros, BestNu=0 | ✅ Done |
| 40 | NSVR-S3 | ❌ Missing | NEW: null profile → ArgumentNullException | ✅ Done |
| 41 | NSVR-S4 | ❌ Missing | NEW: LM22-format TSV parses correctly | ✅ Done |
| 42 | NSVR-S5 | ❌ Missing | NEW: loaded matrix drives deconvolution end-to-end | ✅ Done |
| 43 | NSVR-C1 | ❌ Missing | NEW: empty input → FormatException | ✅ Done |
| 44 | NSVR-C2 | ❌ Missing | NEW: header without cell types → FormatException | ✅ Done |
| 45 | NSVR-C3 | ❌ Missing | NEW: ragged row → FormatException | ✅ Done |
| 46 | NSVR-C4 | ❌ Missing | NEW: non-numeric value → FormatException | ✅ Done |
| 47 | NSVR-C5 | ❌ Missing | NEW: null lines → ArgumentNullException | ✅ Done |

**Implementation fix (prior session):** ssGSEA `ComputeSsGseaScore` changed from expression-value weighting (`|expr|^τ`) to rank-based weighting (`rank^τ`, rank = N−i) per Barbie et al. (2009) / GSVA package (Hänzelmann et al. 2013). Previous weighting produced scores on wrong scale for ESTIMATE purity coefficients.

**Implementation addition (prior session):** added opt-in public `EstimateTumorPurity(double)` applying the Yoshihara (2013) `cos(a + b·score)` transform with NaN-on-negative-cosine domain handling (mirrors ESTIMATE/tidyestimate `ifelse(purity<0, NA, purity)`). Default 5-marker/ssGSEA `EstimateInfiltration` path unchanged.

**Implementation addition (2026-06-25):** added opt-in public `DeconvoluteImmuneCellsNuSvr(...)` — CIBERSORT-style linear-kernel ν-SVR deconvolution (ν ∈ {0.25,0.5,0.75} sweep, lowest-RMSE selection, z-score standardisation, zero-clip + sum-to-1) per Schölkopf et al. (2000) / Newman et al. (2015), plus `LoadSignatureMatrix(...)` (LM22-format TSV loader). The ν-SVR dual is solved by SMO-style coordinate ascent; verified by planted-truth recovery and a scikit-learn/libsvm `NuSVR` cross-check (agreement < 2e-3). The LM22 matrix itself is **not** bundled (Stanford no-redistribution licence) — caller-supplied via the loader. Default `EstimateInfiltration`, `EstimateTumorPurity`, and `DeconvoluteImmuneCells` (NNLS) paths unchanged.

**Total items:** 47
**✅ Done:** 44 | **🗑 Removed:** 3 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M6 (MUST) | ✅ Covered | 6 tests; M5/M6 strengthened with exact values (unique genes proof) |
| M7, M8 | 🗑 Removed | Subsumed by INV-1/2 parametric (same cell types, same invariants) |
| M9–M13 (MUST) | ✅ Covered | 5 tests; M9 strengthened with purity formula verification |
| M14 (MUST) | ✅ Covered | NEW: 3 sub-tests — rank-based ssGSEA exact references |
| S1–S4 (SHOULD) | ✅ Covered | 4 tests; S1 strengthened (exact 0.75:0.25), S3 strengthened (corr=1.0, RMSE=0) |
| C1–C3 (COULD) | ✅ Covered | 4 tests (C2 has 2 sub-tests), all passing |
| INV-1, INV-2 | ✅ Covered | 6 parameterized cases (also covers former M7/M8 scope) |
| INV-3, INV-4 | ✅ Covered | 4 parameterized cases; strengthened with formula verification |
| INV-5 | ✅ Covered | 1 test |
| INV-6 (monotone purity) | ✅ Covered | E6 |
| E1–E7 (`EstimateTumorPurity`) | ✅ Covered | NEW: 7 tests — exact cosine values, NaN domain handling, boundary, monotonicity, closed-form identity |
| NSVR-M1–M5 (MUST) | ✅ Covered | 5 tests; M2 = scikit-learn/libsvm reference (decisive), M1 = planted-truth recovery |
| NSVR-S1–S5 (SHOULD) | ✅ Covered | 5 tests; determinism, no-overlap, null, loader parse, loader→deconvolution |
| NSVR-C1–C5 (loader validation) | ✅ Covered | 5 tests; empty / no-cell-types / ragged / non-numeric / null |
| INV-7, INV-8, INV-9, INV-10 | ✅ Covered | via NSVR-M3 (≥0, Σ=1), NSVR-M4 (ν set), NSVR-S1 (determinism) |

**Total test methods:** 55 (prior 40 + 15 for ν-SVR `DeconvoluteImmuneCellsNuSvr` + `LoadSignatureMatrix`)

---

## 6. Assumption Register

**Total assumptions:** 0

_No assumptions. All algorithms and data structures are precisely documented with external source references._
_Default immune and stromal gene sets are the complete 141+141 ESTIMATE signatures from Yoshihara et al. (2013),_
_extracted from the official ESTIMATE R package v1.0.11 (inst/extdata/SI\_geneset.gmt)._
_Deconvolution offers two engines: the NNLS/LLSR baseline per Lawson & Hanson (1995) / Abbas et al. (2009), and the opt-in CIBERSORT ν-SVR per Schölkopf et al. (2000) / Newman et al. (2015)._
_Signature matrices and gene sets are configurable via API parameters._

---

## 7. Open Questions / Decisions

1. **LM22 not bundled (DECISION, resolved).** The CIBERSORT LM22 signature matrix is distributed by Stanford under a non-commercial licence that forbids redistribution ("RECIPIENT shall not distribute the Program …") and is gated behind registration. Decision: do NOT embed LM22; implement the ν-SVR algorithm + an LM22-format loader (`LoadSignatureMatrix`), bundle only the pre-existing representative 5-marker matrix for tests/default, and require the caller to supply `LM22.txt` under their own CIBERSORT licence.
2. **No bit-exact CIBERSORT-tool parity (honest residual).** The ν-SVR is verified by planted-truth recovery and a scikit-learn/libsvm `NuSVR` cross-check, not by reproducing the official CIBERSORT tool's exact published per-sample fractions — those additionally require LM22 + the tool's full quantile-normalisation/permutation pipeline, which is out of scope.
