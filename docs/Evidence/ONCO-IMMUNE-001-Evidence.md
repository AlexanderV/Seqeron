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

## Assumptions

1. **ASSUMPTION: Simplified signature matrix** — The implementation uses a simplified subset of immune cell signature genes rather than the full LM22 matrix (547 genes × 22 cell types). Justification: A library implementation provides the algorithmic framework; users would supply their own signature matrices for production use.

2. **ASSUMPTION: NNLS instead of ν-SVR** — DeconvoluteImmuneCells uses NNLS (non-negative least squares) rather than ν-SVR for the deconvolution. Justification: NNLS is a well-established mathematical method referenced in the CIBERSORT paper as a baseline; ν-SVR requires additional ML frameworks. Both are valid deconvolution approaches.

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

---

## References

1. Newman AM, Liu CL, Green MR, et al. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods 12(5):453-457. https://doi.org/10.1038/nmeth.3337
2. Yoshihara K, Shahmoradgoli M, Martínez E, et al. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature Communications 4:2612. https://doi.org/10.1038/ncomms3612
3. Becht E, Giraldo NA, Lacroix L, et al. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology 17(1):218. https://doi.org/10.1186/s13059-016-1070-5
4. Subramanian A, Tamayo P, Mootha VK, et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. PNAS 102(43):15545-15550. https://doi.org/10.1073/pnas.0506580102
5. Newman AM, Steen CB, Liu CL, et al. (2019). Determining cell type abundance and expression from bulk tissues with digital cytometry. Nature Biotechnology 37(7):773-782. https://doi.org/10.1038/s41587-019-0114-2

---

## Change History

- **2026-03-06**: Initial documentation.
