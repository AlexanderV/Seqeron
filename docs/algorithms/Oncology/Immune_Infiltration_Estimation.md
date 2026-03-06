# Immune Infiltration Estimation

## Documented Theory

<!-- All content in this section comes from authoritative sources only. -->

### Purpose

Immune infiltration estimation quantifies the degree and composition of immune cell presence within tumor tissue from gene expression data. Two complementary approaches exist:

1. **ESTIMATE** (Estimation of STromal and Immune cells in MAlignant Tumor tissues using Expression data): Computes aggregate immune and stromal enrichment scores using single-sample Gene Set Enrichment Analysis (ssGSEA) on predefined gene signatures, then derives tumor purity (Yoshihara et al., 2013).

2. **CIBERSORT** (Cell-type Identification By Estimating Relative Subsets Of RNA Transcripts): Deconvolutes the proportions of individual immune cell types from bulk expression using a reference signature matrix (LM22: 547 genes × 22 cell types) and ν-support vector regression or non-negative least squares (Newman et al., 2015).

### Core Mechanism

#### ESTIMATE-style Infiltration Scoring

For each gene set $G$ (immune or stromal) and expression profile with $N$ genes:

1. Rank genes by expression value (descending).
2. Compute ssGSEA enrichment score (Barbie et al., 2009):

$$ES = \max_{1 \le i \le N} \left| \sum_{j=1}^{i} \left[ \mathbb{1}(g_j \in G) \cdot \frac{|r_j|}{\sum_{g_k \in G} |r_k|} - \mathbb{1}(g_j \notin G) \cdot \frac{1}{N - |G|} \right] \right|$$

3. Compute ESTIMATE score: $S_{ESTIMATE} = S_{immune} + S_{stromal}$
4. Compute tumor purity (Yoshihara et al., 2013):

$$P = \cos(0.6049872018 + 0.0001467884 \times S_{ESTIMATE})$$

clamped to $[0, 1]$.

#### CIBERSORT-style Deconvolution

The deconvolution model (Newman et al., 2015):

$$\mathbf{m} = \mathbf{S} \cdot \mathbf{f}$$

where $\mathbf{m}$ is the mixture expression vector ($n_{genes} \times 1$), $\mathbf{S}$ is the signature matrix ($n_{genes} \times n_{celltypes}$), and $\mathbf{f}$ is the fraction vector ($n_{celltypes} \times 1$).

NNLS solution (Lawson & Hanson, 1995):

$$\min_{\mathbf{f}} \|\mathbf{m} - \mathbf{S} \cdot \mathbf{f}\|^2 \quad \text{subject to} \quad \mathbf{f} \ge 0$$

then normalize: $\hat{f}_j = f_j / \sum_k f_k$.

### Properties

- **ESTIMATE:** Scores are continuous real numbers; tumor purity is bounded in $[0, 1]$.
- **NNLS Deconvolution:** Fractions are non-negative by construction; normalized to sum to 1.0. Solution exists and is unique when $\mathbf{S}^T\mathbf{S}$ is positive definite.
- **Linear model:** Assumes gene expression is a linear combination of cell-type-specific expression profiles (additive model).

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| ssGSEA Time | $O(N \log N)$ for sorting + $O(N)$ for running sum | Barbie et al. (2009) |
| ssGSEA Space | $O(N)$ | — |
| NNLS Time | $O(p^3 \cdot k)$ where $p$ = passive set size, $k$ = iterations | Lawson & Hanson (1995) |
| NNLS Space | $O(n \cdot m)$ for the signature matrix | — |

---

## Implementation Notes

**Implementation location:** [ImmuneAnalyzer.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs)

- `EstimateInfiltration`: Computes immune and stromal enrichment scores via simplified ssGSEA, then derives ESTIMATE score and tumor purity. Accepts custom gene sets or uses defaults.
- `DeconvoluteImmuneCells`: Deconvolutes immune cell type fractions from bulk expression using NNLS with a configurable signature matrix. Returns fractions, Pearson correlation, and RMSE.

---

## Deviations and Assumptions

| # | Item | Type | Status | Resolution |
|---|------|------|--------|------------|
| 1 | **Simplified signature sets** — Uses 26/18 marker genes instead of full ESTIMATE 141-gene sets | Assumption | ⚠ ASSUMPTION | Library provides algorithmic framework; users supply production gene sets via optional parameters |
| 2 | **NNLS instead of ν-SVR** — Uses NNLS for deconvolution instead of ν-SVR | Assumption | ⚠ ASSUMPTION | NNLS is a baseline method referenced in CIBERSORT paper (Newman et al., 2015); functionally equivalent for most cases |
| 3 | **Simplified LM22** — Uses 5 genes per cell type instead of full LM22 (547 genes × 22 types) | Assumption | ⚠ ASSUMPTION | Full LM22 matrix is proprietary; simplified version demonstrates the algorithm with extensible API |

---

## Sources

- Newman AM, Liu CL, Green MR, et al. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods 12(5):453-457. https://doi.org/10.1038/nmeth.3337
- Yoshihara K, Shahmoradgoli M, Martínez E, et al. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature Communications 4:2612. https://doi.org/10.1038/ncomms3612
- Barbie DA, Tamayo P, Boehm JS, et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature 462:108-112. https://doi.org/10.1038/nature08460
- Subramanian A, Tamayo P, Mootha VK, et al. (2005). Gene set enrichment analysis. PNAS 102(43):15545-15550. https://doi.org/10.1073/pnas.0506580102
- Lawson CL, Hanson RJ (1995). Solving Least Squares Problems. SIAM Classics in Applied Mathematics.
- Newman AM, Steen CB, Liu CL, et al. (2019). Determining cell type abundance and expression from bulk tissues with digital cytometry. Nature Biotechnology 37(7):773-782. https://doi.org/10.1038/s41587-019-0114-2
