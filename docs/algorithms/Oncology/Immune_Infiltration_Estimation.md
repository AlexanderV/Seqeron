# Immune Infiltration Estimation

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-IMMUNE-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Immune infiltration estimation uses bulk gene-expression measurements to quantify either aggregate immune and stromal content or the relative abundance of specific immune cell types within a tumor sample. This document covers two complementary methods implemented in the repository: ESTIMATE-style enrichment scoring and NNLS-based cell-type deconvolution. The first produces aggregate immune, stromal, and tumor-purity scores; the second estimates cell-fraction proportions from a signature matrix. The repository exposes both methods as configurable entry points with default signatures, but the deconvolution side remains a simplified implementation because production-grade behavior still depends on the supplied signature matrix.

## 2. Scientific / Formal Basis

> A = ESTIMATE-style infiltration scoring, B = NNLS-based immune cell deconvolution

### 2.A ESTIMATE-style infiltration scoring

#### Domain Context

The ESTIMATE family of methods infers stromal and immune admixture from tumor gene-expression profiles by testing whether predefined stromal and immune gene signatures are enriched near the top of a ranked expression profile. Yoshihara et al. (2013) combine those enrichment scores into an ESTIMATE score and derive tumor purity from that aggregate score.

#### Core Model

The implemented scoring follows a single-sample GSEA style rank walk: genes are ranked in descending expression order, signature hits receive positive weight, non-signature genes receive negative weight, and the enrichment score is accumulated over the full ranked list. The repository source documents the GSVA-style integral form of ssGSEA with weighting exponent $\tau = 0.25$ (Barbie et al., 2009; Hänzelmann et al., 2013).

The aggregate ESTIMATE score is:

$$
S_{ESTIMATE} = S_{immune} + S_{stromal}
$$

Tumor purity is then computed as (Yoshihara et al., 2013):

$$
P = \cos(0.6049872018 + 0.0001467884 \times S_{ESTIMATE})
$$

with the implementation clamping the final purity value to `[0, 1]`.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-ESTIMATE-01 | The immune and stromal signature genes are appropriate for the tissue and assay being analyzed | Enrichment scores can misrepresent infiltration or stromal content |
| ASM-ESTIMATE-02 | Relative rank ordering of expression values is informative for immune and stromal enrichment | Rank-based enrichment becomes unstable or uninformative when the profile is not comparable across genes |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-ESTIMATE-01 | `EstimateScore = ImmuneScore + StromalScore` | The source computes the aggregate score as the sum of the two ssGSEA scores |
| INV-ESTIMATE-02 | `0 <= TumorPurity <= 1` | The purity helper clamps the cosine-derived value to `[0, 1]` |
| INV-ESTIMATE-03 | Empty overlap for a signature gene set yields a zero contribution from that set | The ssGSEA helper returns `0` when the effective hit set is empty |

### 2.B NNLS-based immune cell deconvolution

#### Domain Context

Cell-type deconvolution treats bulk expression as a mixture of cell-type-specific reference profiles. The approach used here follows the standard linear-mixture formulation associated with digital cytometry and related immune deconvolution methods, while solving the fitting problem with non-negative least squares (Lawson and Hanson, 1995; Abbas et al., 2009; Newman et al., 2015).

#### Core Model

The linear mixture model is:

$$
\mathbf{m} = \mathbf{S} \cdot \mathbf{f}
$$

where $\mathbf{m}$ is the observed mixture expression vector, $\mathbf{S}$ is the signature matrix, and $\mathbf{f}$ is the vector of cell-type fractions.

The repository solves the non-negative least-squares problem:

$$
\min_{\mathbf{f}} \|\mathbf{m} - \mathbf{S}\mathbf{f}\|^2 \quad \text{subject to} \quad \mathbf{f} \ge 0
$$

and then normalizes the solution so that the fraction vector sums to `1` when the raw NNLS solution has positive total mass.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-NNLS-01 | Bulk expression is approximately a linear combination of cell-type-specific expression profiles | Estimated fractions can be biased or uninterpretable |
| ASM-NNLS-02 | The signature matrix is representative of the cell types present in the sample | Missing or mis-specified cell types distort the fitted proportions |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-NNLS-01 | Reported cell fractions are non-negative | NNLS constrains the solution vector to `f >= 0` |
| INV-NNLS-02 | Reported cell fractions sum to `1` when the raw NNLS solution has positive total mass | The implementation explicitly normalizes by `sum(f)` |
| INV-NNLS-03 | `Correlation` and `Rmse` summarize the fit between observed and reconstructed expression | The code reconstructs the mixture profile from the fitted fractions and then computes Pearson correlation and RMSE |

#### Comparison with Related Methods

| Aspect | ESTIMATE-style infiltration scoring | NNLS-based immune cell deconvolution |
|--------|-------------------------------------|--------------------------------------|
| Output granularity | Aggregate immune/stromal scores plus tumor purity | Cell-type fraction vector |
| Input dependency | Gene signatures | Signature matrix |
| Optimization target | Rank-based enrichment | Least-squares fit under non-negativity |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[EstimateInfiltration] expressionProfile` | `IReadOnlyDictionary<string, double>` | required | Gene-expression profile keyed by gene symbol | `null` throws `ArgumentNullException` |
| `[EstimateInfiltration] immuneGenes` | `IReadOnlyList<string>?` | `null` | Optional immune signature gene list | Defaults to the built-in 141-gene ESTIMATE immune signature |
| `[EstimateInfiltration] stromalGenes` | `IReadOnlyList<string>?` | `null` | Optional stromal signature gene list | Defaults to the built-in 141-gene ESTIMATE stromal signature |
| `[DeconvoluteImmuneCells] expressionProfile` | `IReadOnlyDictionary<string, double>` | required | Bulk expression profile for deconvolution | `null` throws `ArgumentNullException` |
| `[DeconvoluteImmuneCells] signatureMatrix` | `IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>>?` | `null` | Optional cell-type signature matrix | Defaults to the built-in 22-cell-type matrix with representative marker genes |
| `[DeconvoluteImmuneCells] maxIterations` | `int` | `1000` | Maximum Lawson-Hanson NNLS iterations | Used by the active-set solver |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[EstimateInfiltration] ImmuneScore` | `double` | ssGSEA-style enrichment score for the immune signature |
| `[EstimateInfiltration] StromalScore` | `double` | ssGSEA-style enrichment score for the stromal signature |
| `[EstimateInfiltration] EstimateScore` | `double` | Sum of immune and stromal scores |
| `[EstimateInfiltration] TumorPurity` | `double` | Cosine-derived purity estimate clamped to `[0, 1]` |
| `[EstimateInfiltration] OverlappingImmuneGenes` | `int` | Number of immune signature genes present in the expression profile |
| `[EstimateInfiltration] OverlappingStromalGenes` | `int` | Number of stromal signature genes present in the expression profile |
| `[DeconvoluteImmuneCells] CellFractions` | `IReadOnlyDictionary<string, double>` | Estimated fractions by cell type |
| `[DeconvoluteImmuneCells] Correlation` | `double` | Pearson correlation between observed and reconstructed expression |
| `[DeconvoluteImmuneCells] Rmse` | `double` | Root mean square reconstruction error |
| `[DeconvoluteImmuneCells] OverlappingGenes` | `int` | Number of signature genes present in the expression profile |

### 3.3 Preconditions and Validation

Both public entry points reject a `null` expression profile with `ArgumentNullException`. `EstimateInfiltration` substitutes built-in immune and stromal signatures when custom lists are not provided. On an empty expression profile, it returns zero immune, stromal, and ESTIMATE scores, but tumor purity is still computed from the zero score by the published cosine formula, yielding approximately `0.8225` after clamping. `DeconvoluteImmuneCells` substitutes a built-in signature matrix when none is provided, and it returns zero fractions, zero correlation, zero RMSE, and zero overlap when there are no overlapping signature genes or no cell types in the matrix.

## 4. Algorithm

### 4.A ESTIMATE-style infiltration scoring

#### High-Level Steps

1. Use the built-in ESTIMATE immune and stromal gene signatures unless custom lists are supplied.
2. Intersect each signature with the expression-profile keys.
3. Rank genes by descending expression value.
4. Compute ssGSEA-style enrichment scores for the immune and stromal overlaps.
5. Add the two enrichment scores to obtain `EstimateScore`.
6. Compute tumor purity with the published cosine formula and clamp it to `[0, 1]`.

#### Decision Rules / Reference Tables

The source uses the GSVA-style integral form of ssGSEA rather than the classic maximum-deviation GSEA statistic. Hit weights are proportional to `rank^tau` with `tau = 0.25`, and miss steps subtract a constant `1 / nMiss` from the running sum.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `EstimateInfiltration` | `O(N log N)` | `O(N)` | Dominated by ranking `N` genes for each scoring pass |

### 4.B NNLS-based immune cell deconvolution

#### High-Level Steps

1. Use the built-in signature matrix unless a custom matrix is supplied.
2. Collect the genes shared between the expression profile and the signature matrix.
3. Assemble the mixture vector and signature matrix for those overlapping genes.
4. Solve the non-negative least-squares problem with the Lawson-Hanson active-set method.
5. Normalize the resulting fractions to sum to `1` when their total is positive.
6. Reconstruct the expression profile and report Pearson correlation and RMSE.

#### Decision Rules / Reference Tables

The built-in signature matrix contains 22 immune cell phenotypes with representative marker genes. The solver operates only on overlapping genes, so cell types remain in the result even when overlap is insufficient for a meaningful fit, but their estimated fractions become zero in the no-overlap branch.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DeconvoluteImmuneCells` | `O(p^3 * k)` | `O(n * m)` | `p` = passive-set size, `k` = NNLS iterations, `n` = overlapping genes, `m` = cell types |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ImmuneAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs)

- `ImmuneAnalyzer.EstimateInfiltration(IReadOnlyDictionary<string, double>, IReadOnlyList<string>?, IReadOnlyList<string>?)`: Computes immune and stromal enrichment plus ESTIMATE score and tumor purity.
- `ImmuneAnalyzer.DeconvoluteImmuneCells(IReadOnlyDictionary<string, double>, IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>>?, int)`: Fits cell-type fractions with NNLS and reports fit diagnostics.

### 5.2 Current Behavior

The repository currently ships the full 141-gene ESTIMATE immune and stromal signatures as defaults. For deconvolution, it ships a 22-cell-type default signature matrix with representative marker genes rather than the full LM22 matrix. The NNLS solver normalizes fitted fractions to sum to `1` and reports both correlation and RMSE against the reconstructed expression profile. Both entry points are configurable, so callers can supply alternative signatures or a more complete reference matrix.

### 5.3 Conformance to Theory / Spec

#### 5.3.A ESTIMATE-style infiltration scoring

**Implemented (verbatim from the cited theory/spec):**

- Rank-based ssGSEA-style enrichment scoring for immune and stromal signatures.
- ESTIMATE score as the sum of immune and stromal enrichment scores.
- Tumor-purity calculation with the Yoshihara et al. cosine formula.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

#### 5.3.B NNLS-based immune cell deconvolution

**Implemented (verbatim from the cited theory/spec):**

- Linear mixture modeling of bulk expression against a signature matrix.
- Non-negative least-squares fitting and normalization of the fitted fraction vector.
- Reconstruction-quality reporting via Pearson correlation and RMSE.

**Intentionally simplified:**

- The built-in default signature matrix uses representative markers rather than the full LM22 reference; **consequence:** out-of-the-box deconvolution is illustrative unless a fuller matrix is supplied.

**Not implemented:**

- ν-SVR-based CIBERSORT fitting; **users should rely on:** external CIBERSORT-family tools or a custom workflow rather than this class.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Representative-marker default matrix instead of full LM22 | Assumption | Default fractions are less detailed than a full LM22-based workflow | accepted | Callers can supply a custom signature matrix |
| 2 | NNLS instead of ν-SVR | Deviation | Deconvolution behavior differs from CIBERSORT's published regression engine | accepted | The source explicitly documents the distinction |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty expression profile for `EstimateInfiltration` | Zero immune/stromal/ESTIMATE scores and tumor purity `≈ 0.8225` | The method evaluates the published purity formula at score `0` |
| No overlapping ESTIMATE signature genes | Corresponding enrichment score is `0` | The ssGSEA helper returns zero when the hit set is empty |
| No overlapping genes for deconvolution | All returned cell fractions are `0`, `Correlation = 0`, `Rmse = 0` | The method exits through a no-overlap branch |

### 6.2 Limitations

The deconvolution path is only as strong as its signature matrix. The built-in ESTIMATE signatures are complete, but the built-in deconvolution matrix is intentionally compact and not a substitute for a production reference such as LM22. The current implementation also uses NNLS rather than ν-SVR, so it should not be described as a direct reimplementation of CIBERSORT.

## 8. References

1. Yoshihara, K., et al. 2013. Inferring tumour purity and stromal and immune cell admixture from expression data. Nature Communications 4:2612. https://doi.org/10.1038/ncomms3612.
2. Barbie, D. A., et al. 2009. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature 462:108-112. https://doi.org/10.1038/nature08460.
3. Hänzelmann, S., R. Castelo, and J. Guinney. 2013. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics 14:7.
4. Lawson, C. L., and R. J. Hanson. 1995. Solving Least Squares Problems. SIAM Classics in Applied Mathematics.
5. Abbas, A. R., et al. 2009. Deconvolution of blood microarray data identifies cellular activation patterns in systemic lupus erythematosus. PLoS One 4:e6098.
6. Newman, A. M., et al. 2015. Robust enumeration of cell subsets from tissue expression profiles. Nature Methods 12(5):453-457. https://doi.org/10.1038/nmeth.3337.
7. Newman, A. M., et al. 2019. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nature Biotechnology 37(7):773-782. https://doi.org/10.1038/s41587-019-0114-2.
8. Subramanian, A., et al. 2005. Gene set enrichment analysis. Proceedings of the National Academy of Sciences 102(43):15545-15550. https://doi.org/10.1073/pnas.0506580102.

