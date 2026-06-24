# Evidence Artifact: ONCO-IMMUNE-001

**Test Unit ID:** ONCO-IMMUNE-001
**Algorithm:** Immune Infiltration Estimation
**Date Collected:** 2026-03-06

---

## Online Sources

### 1. Wikipedia — CIBERSORT

**URL:** https://en.wikipedia.org/wiki/CIBERSORT
**Accessed:** 2026-03-06
**Authority rank:** 4 (Wikipedia with primary source citations)

**Key Extracted Points:**

1. **Definition:** CIBERSORT and CIBERSORTx are bioinformatics tools used to deconvolute cell type proportions and gene expression profiles from bulk RNA sequencing datasets.
2. **Method:** Uses ν-support vector regression (ν-SVR) with a signature matrix (LM22) containing 547 genes distinguishing 22 immune cell types.
3. **Core model:** The deconvolution problem is formulated as: m = S × f, where m is the mixture expression vector, S is the signature matrix, and f is the fraction vector.
4. **Constraint:** Cell type fractions must be non-negative and sum to ≤ 1.0.
5. **Primary citation:** Newman AM et al. (2015). Nature Methods 12(5):453-457. DOI: 10.1038/nmeth.3337

### 2. Newman et al. (2015) — CIBERSORT Original Paper

**URL:** https://doi.org/10.1038/nmeth.3337
**Accessed:** 2026-03-06
**Authority rank:** 1 (Peer-reviewed paper, Nature Methods)

**Key Extracted Points:**

1. **Signature matrix LM22:** Contains 547 genes that distinguish 22 human hematopoietic cell phenotypes: naive B cells, memory B cells, plasma cells, CD8+ T cells, naive CD4+ T cells, resting memory CD4+ T cells, activated memory CD4+ T cells, follicular helper T cells, regulatory T cells (Tregs), gamma delta T cells, resting NK cells, activated NK cells, monocytes, M0 macrophages, M1 macrophages, M2 macrophages, resting dendritic cells, activated dendritic cells, resting mast cells, activated mast cells, eosinophils, and neutrophils.
2. **Deconvolution algorithm:** Linear model: m = S × f; solved using ν-SVR or alternatively non-negative least squares (NNLS).
3. **NNLS formulation:** min ||m - S·f||² subject to f ≥ 0, then normalize f so sum(f) = 1.
4. **Validation:** Tested against flow cytometry ground truth; CIBERSORT outperformed other methods (NNLS, QP, LLSR) in benchmarks.
5. **Statistical significance:** Uses Monte Carlo permutation to generate a p-value for each deconvolution result.

### 3. Yoshihara et al. (2013) — ESTIMATE Algorithm

**URL:** https://doi.org/10.1038/ncomms3612
**Accessed:** 2026-03-06
**Authority rank:** 1 (Peer-reviewed paper, Nature Communications)

**Key Extracted Points:**

1. **ESTIMATE:** Estimation of STromal and Immune cells in MAlignant Tumor tissues using Expression data.
2. **Method:** Uses single-sample GSEA (ssGSEA) to compute enrichment scores for immune and stromal gene signatures.
3. **Immune signature:** 141 genes representing immune cell infiltration markers.
4. **Stromal signature:** 141 genes representing stromal cell markers.
5. **Outputs:** Immune score, stromal score, ESTIMATE score (sum of immune + stromal), tumor purity estimate.
6. **Tumor purity formula:** Purity = cos(0.6049872018 + 0.0001467884 × ESTIMATE_score).
7. **ssGSEA enrichment score:** For each sample, rank genes by expression, compute running sum statistic similar to GSEA but at single-sample level.

**Re-retrieval (2026-06-24) — verbatim purity transform + domain (for `EstimateTumorPurity`):**

- **Search:** `Yoshihara 2013 ESTIMATE tumor purity cos(0.6049872018 + 0.0001467884 ESTIMATEScore) formula`. Snippet (Nature Communications / Verhaak Lab): *"The tumor purity was estimated by the formula: Tumor_purity = cos (0.6049872018 + 0.0001467884\*ESTIMATEScore)."*
- **Fetched** `https://search.r-project.org/CRAN/refmans/hacksig/html/hack_estimate.html` (R `hacksig::hack_estimate` doc, a re-implementation of the ESTIMATE transform). Verbatim: **"Purity = cos(0.6049872018 + 0.0001467884 * ESTIMATE)"**. No additional domain restriction stated there.
- **Fetched** `https://www.aging-us.com/article/203714/text` (Aging journal). Verbatim: **"Tumor_purity = cos (0.6049872018 + 0.0001467884\*ESTIMATEScore)"** — coefficients confirmed identical, third independent source.
- **Search + fetch** of the ESTIMATE/`tidyestimate` reference implementation. The biostars / open-source-biology discussion states: *"This formula was derived using only Affymetrix data, and therefore cannot be used to convert RNAseq-based ESTIMATE score to tumor purity."* The PMC full text confirms: *"established a regression curve for ESTIMATE score and tumor purity based on ABSOLUTE in the TCGA data set by applying the nonlinear least squares method."*
- **Fetched** `https://raw.githubusercontent.com/KaiAragaki/tidyestimate/main/R/estimate_score.R` (reference implementation source). Verbatim R:
  ```r
  if (is_affymetrix) {
      scores <- scores |>
          dplyr::mutate(purity = cos(0.6049872018 + 0.0001467884 * .data$estimate),
                        purity = ifelse(.data$purity < 0, NA, .data$purity))
  }
  ```
  → **Decisive domain rule:** the reference implementation computes `cos(a + b·score)` and sets the purity to **NA when it is negative** (out of the calibrated domain), and only computes it at all when `is_affymetrix = TRUE`. The `tidyestimate` vignette states: *"As the data used to train the model to convert ESTIMATE scores to purity scores were produced by Affymetrix, it is unwise to infer a tumor purity score using the same method for RNAseq data."*
8. **Domain/calibration (for the opt-in absolute-purity method):** the cosine model is fit by nonlinear least squares against ABSOLUTE purity on TCGA Affymetrix data; valid only for Affymetrix-derived scores; negative cosine values are out of domain → NA (mirrored as `double.NaN`). cos is monotone-decreasing on [0, π], so purity strictly decreases with ESTIMATE score over the calibrated range.

### 4. Wikipedia — Tumor Microenvironment

**URL:** https://en.wikipedia.org/wiki/Tumor_microenvironment
**Accessed:** 2026-03-06
**Authority rank:** 4 (Wikipedia with primary source citations)

**Key Extracted Points:**

1. **Immune cell types in TME:** CD8+ cytotoxic T cells, CD4+ T helper cells, regulatory T cells (Tregs), B cells, NK cells, tumor-associated macrophages (M1 pro-inflammatory, M2 anti-inflammatory), myeloid-derived suppressor cells (MDSCs), neutrophils, dendritic cells.
2. **Functional significance:** High CD8+ T cell infiltration generally associated with good prognosis; high Tregs and M2 macrophages associated with immune evasion and poor prognosis.
3. **Primary cited sources:** Anderson & Simon (2020), Joyce & Fearon (2015), Lei et al. (2020).

### 5. Wikipedia — Gene Set Enrichment Analysis

**URL:** https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis
**Accessed:** 2026-03-06
**Authority rank:** 4 (Wikipedia with primary source citations)

**Key Extracted Points:**

1. **GSEA method:** Ranks genes by expression, computes enrichment score as Kolmogorov-Smirnov-like statistic.
2. **ssGSEA variant:** Single-sample extension for scoring individual samples.
3. **Enrichment score formula:** ES = max deviation of running sum = Phit(S,i) - Pmiss(S,i).
4. **Primary citation:** Subramanian A et al. (2005). PNAS 102(43):15545-15550. DOI: 10.1073/pnas.0506580102

### 6. Becht et al. (2016) — MCP-counter

**URL:** https://doi.org/10.1186/s13059-016-1070-5
**Accessed:** 2026-03-06
**Authority rank:** 1 (Peer-reviewed paper, Genome Biology)

**Key Extracted Points:**

1. **Method:** Uses transcriptomic markers specific to immune and stromal cell populations.
2. **Cell types scored:** CD3+ T cells, CD8+ T cells, cytotoxic lymphocytes, NK cells, B lineage, monocytic lineage, myeloid dendritic cells, neutrophils, endothelial cells, fibroblasts.
3. **Scoring:** Computes geometric mean of marker gene expression values for each cell type.
4. **Interpretation:** Produces abundance scores (not proportions); scores are relative within a cell type across samples.

---

### 7. Schölkopf, Smola, Williamson & Bartlett (2000) — ν-SVR formulation

**URL:** https://pubmed.ncbi.nlm.nih.gov/10905814/ (citation + abstract); ν-SVR equations from the Smola & Schölkopf (2004) tutorial PDF https://alex.smola.org/papers/2003/SmoSch03b.pdf (eqs (60)–(62), retrieved and read page-by-page this session)
**Accessed:** 2026-06-25
**Authority rank:** 1 (Peer-reviewed: Neural Computation 12(5):1207-1245, doi:10.1162/089976600300015565)

**How retrieved:** WebSearch ("Schölkopf New Support Vector Algorithms 2000 nu-SVR dual quadratic program"); fetched the PubMed page for the citation; downloaded the Smola/Schölkopf tutorial PDF and read pages 1, 3–8, 12–14 directly for the verbatim equations.

**Key Extracted Points:**

1. **ν-SVR primal (eq 60–61):** minimise `R_ν[f] := R_emp[f] + (λ/2)‖w‖² + νε`, equivalently `min ½‖w‖² + C·(Σ_i(c̃(ξ_i)+c̃(ξ_i*)) + ℓ·ν·ε)` subject to `y_i − ⟨w,x_i⟩ − b ≤ ε + ξ_i`, `⟨w,x_i⟩ + b − y_i ≤ ε + ξ_i*`, `ξ_i, ξ_i* ≥ 0`. The accuracy parameter ε becomes a variable to be minimised, traded off by ν.
2. **ν-SVR dual (eq 62):** maximise `−½ Σ_{i,j}(α_i−α_i*)(α_j−α_j*)·k(x_i,x_j) + Σ_i y_i(α_i−α_i*)` subject to `Σ_i(α_i−α_i*) = 0`, `Σ_i(α_i+α_i*) ≤ C·ν·ℓ`, and `α_i, α_i* ∈ [0, C]`. (Identical to the ε-SVR dual except ε drops out of the objective and the extra ν-budget constraint is added.)
3. **Primal recovery:** `w = Σ_i(α_i−α_i*)·x_i` and `f(x) = Σ_i(α_i−α_i*)·k(x_i,x) + b` (tutorial, after eq 18). For the linear kernel `k(x_i,x_j)=⟨x_i,x_j⟩`.
4. **Theorem 9 (Schölkopf et al. 2000), tutorial p.12:** (1) ν is an upper bound on the fraction of errors; (2) ν is a lower bound on the fraction of support vectors; (3) asymptotically, ν equals both the fraction of SVs and the fraction of errors.

---

### 8. Newman et al. (2015) + CIBERSORT protocol (Chen et al., 2018) — ν-SVR deconvolution method

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5895181/ (Chen B et al., "Profiling Tumor Infiltrating Immune Cells with CIBERSORT", Methods Mol Biol 2018; 1711:243-259)
**Accessed:** 2026-06-25
**Authority rank:** 1 (Peer-reviewed protocol of the Newman 2015 method)

**How retrieved:** WebFetch of the PMC article with a prompt extracting the ν sweep, selection metric, post-processing, and LM22 dimensions.

**Key Extracted Points:**

1. **ν sweep:** "CIBERSORT uses a set of ν values (0.25, 0.5, 0.75) and chooses the value producing the best performance."
2. **Selection metric:** the ν chosen is the one giving "the lowest root mean square between **m** and the deconvolution result **f** × **B**" (RMSE between observed mixture and reconstruction).
3. **Output:** "by default, deconvolution results are expressed as relative fractions normalized to 1" — i.e. the regression coefficients are zero-clipped (negatives → 0) and the remaining coefficients normalised to sum 1.
4. **LM22 dimensions:** "a signature matrix file consisting of 547 genes that accurately distinguish 22 mature human hematopoietic populations" → 547 genes × 22 cell types.
5. **Standardisation:** the mixture and signature expression are z-score standardised before regression (Newman et al., 2015, Online Methods).

---

### 9. CIBERSORT / LM22 LICENCE (decisive for data handling)

**URL:** https://cibersort.stanford.edu (registration gate); licence text verbatim from the archived CIBERSORT licence https://gist.github.com/dhimmel/58dcd9b512e669f20a65ddf73997b733 ; corroborated by immunedeconv `set_cibersort_mat` docs (LM22 must be downloaded by the user from the CIBERSORT website).
**Accessed:** 2026-06-25
**Authority rank:** 1 (the licence governing the artifact)

**How retrieved:** WebSearch ("CIBERSORT LM22 license terms"); WebFetch of the licence gist extracting the restrictive clauses verbatim.

**Key Extracted Points (verbatim licence clauses):**

1. **No redistribution:** *"RECIPIENT shall not distribute the Program or transfer it to any other person or organization without prior written permission from STANFORD."*
2. **Non-commercial only:** *"RECIPIENT shall not use the Program on behalf of any organization that is not a non-profit organization. RECIPIENT shall not use the Program for commercial advantage, or in the course of for-profit activities."*
3. **No modification:** *"RECIPIENT shall NOT make modifications to the Program."*
4. **Access gate:** LM22 (`LM22.txt`) is only obtainable after registration at https://cibersort.stanford.edu; it is not publicly downloadable as plaintext.

**LICENCE DECISION:** LM22 is **NOT redistributable** (contrast Pfam/CC0). Per the mission-critical data-handling rule, LM22 is therefore **not embedded** in this library. Instead: the ν-SVR algorithm and an **LM22-format loader** (`LoadSignatureMatrix`) are implemented; the caller supplies `LM22.txt` under their own CIBERSORT licence. Only a small synthetic/representative signature matrix is bundled (the pre-existing 5-marker `DefaultSignatureMatrix`), used for tests and as a non-LM22 default.

---

## Documented Corner Cases and Failure Modes

### From Newman et al. (2015) — CIBERSORT

1. **Low signal-to-noise:** When expression profiles have high noise, deconvolution accuracy degrades significantly.
2. **Missing signature genes:** If a substantial fraction of signature genes are missing from the expression profile, deconvolution results become unreliable.
3. **Collinear cell types:** Closely related cell types (e.g., resting vs. activated NK cells) are harder to distinguish.
4. **Non-hematopoietic content:** CIBERSORT estimates only hematopoietic fractions; non-immune cells are represented in the residual.

### From Yoshihara et al. (2013) — ESTIMATE

1. **Low purity tumors:** ESTIMATE score may exceed the valid range for tumor purity calculation at extreme values.
2. **Non-solid tumors:** ESTIMATE was designed for solid tumors; hematological malignancies are not appropriate inputs.

### General

1. **Empty expression profile:** No genes → no scores can be computed.
2. **No overlapping genes:** If signature genes don't overlap with the input expression profile, all scores should be zero/undefined.
3. **Negative expression values:** Log-transformed expression data may contain negative values; these are valid inputs.
4. **Single gene overlap:** Minimal overlap should produce minimal but non-zero scores.

---

## Test Datasets

### Dataset 1: Synthetic Pure Immune Cell Type

**Source:** Derived from CIBERSORT LM22 concept — pure cell type should deconvolve to 100% that type.

| Parameter | Value |
|-----------|-------|
| Cell type | CD8 T cells |
| Gene count | 5 marker genes |
| Expected fraction | 1.0 for CD8 T cells, 0.0 for others |

### Dataset 2: Equal Mixture of Two Cell Types

**Source:** Mathematical derivation — average of two signature columns.

| Parameter | Value |
|-----------|-------|
| Cell types | B cells + T cells |
| Mixture ratio | 50:50 |
| Expected fractions | 0.5 each |

### Dataset 3: ESTIMATE-like Immune Scoring

**Source:** Derived from Yoshihara et al. (2013) ssGSEA concept.

| Parameter | Value |
|-----------|-------|
| Gene set | Immune signature genes |
| High immune sample | High expression of immune markers |
| Low immune sample | Low expression of immune markers |
| Expected | High sample immune score > Low sample immune score |

---

### Dataset 4: ν-SVR planted-truth recovery (synthetic mixture)

**Source:** Constructed this session per the CIBERSORT linear-mixture model `m = B·f` (Newman et al., 2015).

| Parameter | Value |
|-----------|-------|
| Signature B | bundled 5-marker × 22-cell-type `DefaultSignatureMatrix` |
| Planted fractions f | T_cells_CD8 = 0.60, B_cells_naive = 0.30, Monocytes = 0.10 |
| Bulk m | `m = B·f` (85 overlapping genes) |
| Expected | `DeconvoluteImmuneCellsNuSvr(m)` recovers f within ≈0.005 per type; all others ≈0; Σ = 1 |
| Observed (this session) | CD8 0.5971, B_naive 0.2989, Monocytes 0.1040; corr 0.99997 |

### Dataset 5: scikit-learn `NuSVR` cross-implementation reference

**Source:** scikit-learn 1.6.1 `NuSVR(kernel='linear', nu, C=1.0)` (libsvm backend), run this session on the standardised problem below.

| Parameter | Value |
|-----------|-------|
| Signature (3 cell types × 3 disjoint markers each) | TypeA {a1:10,a2:8,a3:6}, TypeB {b1:9,b2:7,b3:5}, TypeC {c1:11,c2:4,c3:8} |
| Planted f | TypeA 0.5, TypeB 0.2, TypeC 0.3 |
| Standardisation | per-column z-score of B and z-score of m (population sd) |
| Selected ν (lowest RMSE) | 0.75 |
| sklearn normalised fractions | TypeA 0.508497, TypeB 0.179491, TypeC 0.312012 |
| This implementation | TypeA 0.50846, TypeB 0.17956, TypeC 0.31198 (agreement < 2×10⁻³) |

---

## Assumptions

1. **ASSUMPTION: Simplified signature matrix** — The implementation uses a simplified subset of immune cell signature genes rather than the full LM22 matrix (547 genes × 22 cell types). Justification: A library implementation provides the algorithmic framework; users would supply their own signature matrices for production use.

2. **RESOLVED (2026-06-25): ν-SVR now implemented.** The CIBERSORT ν-support-vector-regression deconvolution (Newman et al., 2015; Schölkopf et al., 2000) is implemented as the opt-in `DeconvoluteImmuneCellsNuSvr` (sweeps ν ∈ {0.25, 0.5, 0.75}, selects lowest-RMSE, zero-clips and normalises to sum 1). `DeconvoluteImmuneCells` retains the NNLS/LLSR baseline (Abbas et al., 2009) unchanged. The ν-SVR solver was verified by (a) planted-truth recovery and (b) a cross-implementation match against scikit-learn 1.6.1 `NuSVR` (libsvm). **Residual:** the LM22 matrix itself is caller-supplied (Stanford licence forbids redistribution — see Online Source 9); exact reproduction of CIBERSORT's published per-sample fractions also requires LM22 + the tool's full quantile-normalisation pipeline, which is out of scope.

3. **ASSUMPTION: Simplified ssGSEA** — EstimateInfiltration uses a simplified ssGSEA scoring approach (mean expression of signature genes with rank-based weighting). Justification: The core ssGSEA algorithm is well-defined but the full implementation requires the GSEA enrichment score computation. A simplified enrichment scoring retains the essential concept.

---

## Recommendations for Test Coverage

1. **MUST Test:** Empty/null expression profile → returns zero/empty results — Evidence: Standard robustness
2. **MUST Test:** Single cell type signature → deconvolves to 100% that type — Evidence: CIBERSORT mathematical identity (Newman et al., 2015)
3. **MUST Test:** Equal mixture of two cell types → ~50:50 fractions — Evidence: CIBERSORT linear model property (Newman et al., 2015)
4. **MUST Test:** Cell fractions sum to 1.0 (invariant) — Evidence: CIBERSORT normalization (Newman et al., 2015)
5. **MUST Test:** All fractions ≥ 0 (invariant) — Evidence: NNLS non-negativity constraint (Newman et al., 2015)
6. **MUST Test:** Immune score ordering: high-immune > low-immune sample — Evidence: ESTIMATE concept (Yoshihara et al., 2013)
7. **MUST Test:** No overlapping genes → zero scores — Evidence: Mathematical definition
8. **SHOULD Test:** Unequal mixture → proportional fractions — Rationale: Validates linearity of deconvolution
9. **SHOULD Test:** Expression profile with genes not in signature → ignored — Rationale: Robustness to extra genes
10. **COULD Test:** Large number of cell types → all fractions computed — Rationale: Scalability
11. **MUST Test (ν-SVR):** Planted-truth `m = B·f` → `DeconvoluteImmuneCellsNuSvr` recovers f — Evidence: Newman et al. (2015) linear mixture (Dataset 4)
12. **MUST Test (ν-SVR):** Match scikit-learn/libsvm `NuSVR` on a small standardised problem — Evidence: cross-implementation reference (Dataset 5); Schölkopf et al. (2000)
13. **MUST Test (ν-SVR):** Selected ν ∈ {0.25, 0.5, 0.75} by lowest RMSE — Evidence: CIBERSORT protocol (Chen et al., 2018)
14. **MUST Test (ν-SVR):** Fractions ≥ 0 and Σ = 1 (zero-clip + normalise) — Evidence: Newman et al. (2015)
15. **MUST Test (loader):** LM22-format TSV parses to cellType→(gene→value); rejects empty/ragged/non-numeric — Evidence: LM22 file format (Newman et al., 2015)

---

## References

1. Newman AM, Liu CL, Green MR, et al. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods 12(5):453-457. https://doi.org/10.1038/nmeth.3337
2. Yoshihara K, Shahmoradgoli M, Martínez E, et al. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature Communications 4:2612. https://doi.org/10.1038/ncomms3612
3. Becht E, Giraldo NA, Lacroix L, et al. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology 17(1):218. https://doi.org/10.1186/s13059-016-1070-5
4. Subramanian A, Tamayo P, Mootha VK, et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. PNAS 102(43):15545-15550. https://doi.org/10.1073/pnas.0506580102
5. Newman AM, Steen CB, Liu CL, et al. (2019). Determining cell type abundance and expression from bulk tissues with digital cytometry. Nature Biotechnology 37(7):773-782. https://doi.org/10.1038/s41587-019-0114-2
6. Schölkopf B, Smola AJ, Williamson RC, Bartlett PL (2000). New support vector algorithms. Neural Computation 12(5):1207-1245. https://doi.org/10.1162/089976600300015565 (PMID 10905814). ν-SVR equations from Smola AJ, Schölkopf B (2004), A tutorial on support vector regression, Statistics and Computing 14:199-222, eqs (60)–(62): https://alex.smola.org/papers/2003/SmoSch03b.pdf
7. Chen B, Khodadoust MS, Liu CL, Newman AM, Alizadeh AA (2018). Profiling Tumor Infiltrating Immune Cells with CIBERSORT. Methods Mol Biol 1711:243-259. https://pmc.ncbi.nlm.nih.gov/articles/PMC5895181/
8. CIBERSORT licence (Stanford University). Verbatim clauses: https://gist.github.com/dhimmel/58dcd9b512e669f20a65ddf73997b733 ; registration/download gate: https://cibersort.stanford.edu
9. Pedregosa F, et al. (2011). Scikit-learn: Machine Learning in Python (NuSVR, libsvm backend). JMLR 12:2825-2830. Used as the cross-implementation ν-SVR reference (v1.6.1).

---

## Change History

- **2026-03-06**: Initial documentation.
- **2026-06-25**: Added Online Sources 7–9 (ν-SVR formulation per Schölkopf 2000 / Smola tutorial eqs 60–62; CIBERSORT ν sweep + RMSE selection + zero-clip/sum-to-1; LM22 licence decision — caller-supplied, not redistributable). Resolved Assumption 2 (ν-SVR implemented as opt-in `DeconvoluteImmuneCellsNuSvr` + `LoadSignatureMatrix`). Added planted-truth + scikit-learn `NuSVR` cross-check datasets.
