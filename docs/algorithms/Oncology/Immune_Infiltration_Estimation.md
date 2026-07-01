# Immune Infiltration Estimation

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-IMMUNE-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-25 |

## 1. Overview

Immune infiltration estimation uses bulk gene-expression measurements to quantify either aggregate immune and stromal content or the relative abundance of specific immune cell types within a tumor sample. This document covers three complementary methods implemented in the repository: ESTIMATE-style enrichment scoring, NNLS-based cell-type deconvolution, and CIBERSORT-style ν-support-vector-regression (ν-SVR) deconvolution. The first produces aggregate immune, stromal, and tumor-purity scores; the second and third estimate cell-fraction proportions from a signature matrix (NNLS as the LLSR baseline, ν-SVR as the CIBERSORT regression engine). The repository exposes all methods as configurable entry points with default signatures, but the deconvolution side remains a simplified implementation because production-grade behavior still depends on the supplied signature matrix — in particular, the CIBERSORT LM22 matrix is caller-supplied for licence reasons (Section 6.2).

## 2. Scientific / Formal Basis

> A = ESTIMATE-style infiltration scoring, B = NNLS-based immune cell deconvolution, C = CIBERSORT ν-SVR deconvolution

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

### 2.C CIBERSORT ν-SVR deconvolution

#### Domain Context

CIBERSORT (Newman et al., 2015) estimates immune cell-type fractions from bulk expression using the same linear-mixture model as NNLS, but fits it with **ν-support-vector regression** (ν-SVR; Schölkopf et al., 2000) rather than least squares. ν-SVR is an ε-insensitive, robust regression: it ignores residuals inside an automatically-sized ε-tube and bounds the influence of outlier genes, which Newman et al. report makes deconvolution robust to noise, unknown mixture content, and signature mis-specification. CIBERSORT is normally run with the LM22 signature matrix (547 genes × 22 cell types).

#### Core Model

For mixture vector $\mathbf{m}$ and signature matrix $\mathbf{B}$ (genes × cell types), a linear-kernel ν-SVR of $\mathbf{m}$ on the columns of $\mathbf{B}$ recovers a weight vector $\mathbf{f} = \mathbf{w}$. The ν-SVR primal (Schölkopf et al., 2000; Smola & Schölkopf tutorial eqs 60–61) is:

$$
\min_{\mathbf{w}, b, \varepsilon, \xi}\; \tfrac{1}{2}\lVert\mathbf{w}\rVert^2 + C\Big(\tfrac{1}{\ell}\sum_i(\xi_i+\xi_i^{*}) + \nu\,\varepsilon\Big)
\quad\text{s.t.}\quad
\begin{cases} y_i-\langle\mathbf{w},x_i\rangle-b \le \varepsilon+\xi_i\\ \langle\mathbf{w},x_i\rangle+b-y_i \le \varepsilon+\xi_i^{*}\\ \xi_i,\xi_i^{*},\varepsilon \ge 0\end{cases}
$$

whose linear-kernel dual (eq 62) is:

$$
\max_{\alpha,\alpha^{*}}\; -\tfrac{1}{2}\sum_{i,j}(\alpha_i-\alpha_i^{*})(\alpha_j-\alpha_j^{*})\langle x_i,x_j\rangle + \sum_i y_i(\alpha_i-\alpha_i^{*})
\quad\text{s.t.}\quad
\sum_i(\alpha_i-\alpha_i^{*})=0,\;\; \sum_i(\alpha_i+\alpha_i^{*})\le C\nu\ell,\;\; \alpha_i,\alpha_i^{*}\in[0,C]
$$

with primal recovery $\mathbf{w}=\sum_i(\alpha_i-\alpha_i^{*})x_i$. CIBERSORT sweeps $\nu \in \{0.25, 0.5, 0.75\}$, selects the $\nu$ that minimises the RMSE between $\mathbf{m}$ and $\mathbf{B}\cdot\mathbf{f}$, then clips negative weights to 0 and normalises the remainder to sum to 1 (Chen et al., 2018; Newman et al., 2015). Mixture and signature are z-score standardised before regression.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-NUSVR-01 | Bulk expression is approximately a linear combination of cell-type reference profiles | Estimated fractions can be biased or uninterpretable |
| ASM-NUSVR-02 | The signature matrix (e.g. LM22) is representative of the cell types present | Missing or mis-specified cell types distort the fitted proportions |
| ASM-NUSVR-03 | ν ∈ (0, 1]; one of {0.25, 0.5, 0.75} yields a usable ε-tube for the data | An empty support set or degenerate tube yields an all-zero weight vector |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-NUSVR-01 | Reported cell fractions are non-negative | Negative ν-SVR weights are clipped to 0 |
| INV-NUSVR-02 | Reported cell fractions sum to `1` when total post-clip mass is positive | The implementation normalizes by `sum(f)` |
| INV-NUSVR-03 | The selected `BestNu` is one of {0.25, 0.5, 0.75} | The sweep only considers `CibersortNuValues` and keeps the lowest-RMSE fit |
| INV-NUSVR-04 | The dual respects `Σ(α_i−α_i*) = 0`, `Σ(α_i+α_i*) ≤ Cνℓ`, `α_i,α_i* ∈ [0,C]` | The SMO step keeps Σβ=0 exactly and clips each step against the box and the ν budget |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[EstimateInfiltration] expressionProfile` | `IReadOnlyDictionary<string, double>` | required | Gene-expression profile keyed by gene symbol | `null` throws `ArgumentNullException` |
| `[EstimateInfiltration] immuneGenes` | `IReadOnlyList<string>?` | `null` | Optional immune signature gene list | Defaults to the built-in 141-gene ESTIMATE immune signature |
| `[EstimateInfiltration] stromalGenes` | `IReadOnlyList<string>?` | `null` | Optional stromal signature gene list | Defaults to the built-in 141-gene ESTIMATE stromal signature |
| `[EstimateTumorPurity] estimateScore` | `double` | required | Affymetrix/ABSOLUTE-calibrated ESTIMATE score (immune + stromal) | Negative-cosine (out-of-domain) scores yield `NaN` |
| `[DeconvoluteImmuneCells] expressionProfile` | `IReadOnlyDictionary<string, double>` | required | Bulk expression profile for deconvolution | `null` throws `ArgumentNullException` |
| `[DeconvoluteImmuneCells] signatureMatrix` | `IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>>?` | `null` | Optional cell-type signature matrix | Defaults to the built-in 22-cell-type matrix with representative marker genes |
| `[DeconvoluteImmuneCells] maxIterations` | `int` | `1000` | Maximum Lawson-Hanson NNLS iterations | Used by the active-set solver |
| `[DeconvoluteImmuneCellsNuSvr] expressionProfile` | `IReadOnlyDictionary<string, double>` | required | Bulk mixture vector for ν-SVR deconvolution | `null` throws `ArgumentNullException`; a non-finite (NaN/±Infinity) value throws `ArgumentException` |
| `[DeconvoluteImmuneCellsNuSvr] signatureMatrix` | `IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>>?` | `null` | Optional signature matrix | Defaults to the built-in representative matrix (NOT LM22); supply LM22 via `LoadSignatureMatrix` |
| `[DeconvoluteImmuneCellsNuSvr] nuValues` | `IReadOnlyList<double>?` | `null` | ν values to sweep | Defaults to `CibersortNuValues` = {0.25, 0.5, 0.75} |
| `[LoadSignatureMatrix] tsvLines` | `IEnumerable<string>` | required | Lines of an LM22-format TSV (header + one row per gene) | `null` throws `ArgumentNullException`; malformed → `FormatException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[EstimateInfiltration] ImmuneScore` | `double` | ssGSEA-style enrichment score for the immune signature |
| `[EstimateInfiltration] StromalScore` | `double` | ssGSEA-style enrichment score for the stromal signature |
| `[EstimateInfiltration] EstimateScore` | `double` | Sum of immune and stromal scores |
| `[EstimateInfiltration] TumorPurity` | `double` | Cosine-derived purity estimate clamped to `[0, 1]` |
| `[EstimateInfiltration] OverlappingImmuneGenes` | `int` | Number of immune signature genes present in the expression profile |
| `[EstimateInfiltration] OverlappingStromalGenes` | `int` | Number of stromal signature genes present in the expression profile |
| `[EstimateTumorPurity] return` | `double` | Absolute tumor purity in `(0, 1]`, or `NaN` when the score is out of the calibrated (non-negative-cosine) domain |
| `[DeconvoluteImmuneCells] CellFractions` | `IReadOnlyDictionary<string, double>` | Estimated fractions by cell type |
| `[DeconvoluteImmuneCells] Correlation` | `double` | Pearson correlation between observed and reconstructed expression |
| `[DeconvoluteImmuneCells] Rmse` | `double` | Root mean square reconstruction error |
| `[DeconvoluteImmuneCells] OverlappingGenes` | `int` | Number of signature genes present in the expression profile |
| `[DeconvoluteImmuneCellsNuSvr] CellFractions` | `IReadOnlyDictionary<string, double>` | ν-SVR cell-type fractions (zero-clipped, summing to 1) |
| `[DeconvoluteImmuneCellsNuSvr] BestNu` | `double` | The ν (∈ {0.25, 0.5, 0.75}) selected as lowest-RMSE (0 when no fit was performed) |
| `[DeconvoluteImmuneCellsNuSvr] Correlation` / `Rmse` / `OverlappingGenes` | `double` / `double` / `int` | Fit diagnostics on the original scale and overlap count |
| `[LoadSignatureMatrix] return` | `IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>>` | Signature matrix: cell type → (gene → value), ready for either deconvolution method |

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

### 4.C CIBERSORT ν-SVR deconvolution

#### High-Level Steps

1. Use the built-in representative signature matrix unless a custom matrix (e.g. caller-supplied LM22) is provided.
2. Collect the genes shared between the mixture and the signature matrix; build mixture vector `m` and matrix `B`.
3. Z-score standardise `m` and each column of `B`.
4. For each ν ∈ {0.25, 0.5, 0.75}, solve the linear-kernel ν-SVR dual and recover the weight vector `f = w`.
5. Reconstruct `B·f` and keep the ν with the lowest RMSE between `m` and `B·f`.
6. Clip negative weights to 0 and normalise the remainder to sum to 1.
7. Report fractions, the selected ν, and fit diagnostics (correlation, RMSE) on the original scale.

#### Decision Rules / Reference Tables

The ν-SVR dual is solved by an SMO-style pairwise coordinate ascent on the dual variables `β_i = α_i − α_i*`: each step moves a working pair `(β_p, β_q)` by `(+δ, −δ)`, which keeps the equality constraint `Σβ_i = 0` exactly; `δ` is the unconstrained optimum `(g_p − g_q)/(K_pp + K_qq − 2K_pq)` clipped against the box `|β_i| ≤ C` and the ν budget `Σ|β_i| ≤ Cνℓ`. The regularisation constant is `C = 1` (`NuSvrCost`, libsvm default). LM22 is the canonical 547-gene × 22-cell-type signature but is caller-supplied (Section 6.2).

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DeconvoluteImmuneCellsNuSvr` | `O(\|ν\| · (n² + n·t·m))` | `O(n² + n·m)` | `n` = overlapping genes, `m` = cell types, `t` = SMO iterations (≤ 200·n), `\|ν\|` = number of ν swept; dominated by the `n×n` kernel and the SMO loop |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ImmuneAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs)

- `ImmuneAnalyzer.EstimateInfiltration(IReadOnlyDictionary<string, double>, IReadOnlyList<string>?, IReadOnlyList<string>?)`: Computes immune and stromal enrichment plus ESTIMATE score and a *relative* tumor purity.
- `ImmuneAnalyzer.EstimateTumorPurity(double)`: Opt-in closed-form transform from a cohort-/Affymetrix-scaled ESTIMATE score to an *absolute* tumor purity via the Yoshihara et al. (2013) cosine model; returns `NaN` for out-of-domain (negative-cosine) scores.
- `ImmuneAnalyzer.DeconvoluteImmuneCells(IReadOnlyDictionary<string, double>, IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>>?, int)`: Fits cell-type fractions with NNLS and reports fit diagnostics.
- `ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(IReadOnlyDictionary<string, double>, IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>>?, IReadOnlyList<double>?)`: Opt-in CIBERSORT-style ν-SVR deconvolution — sweeps ν ∈ {0.25, 0.5, 0.75}, selects lowest-RMSE, zero-clips and normalises to sum 1.
- `ImmuneAnalyzer.LoadSignatureMatrix(IEnumerable<string>)`: Parses an LM22-format TSV (header + gene rows) into the `cellType → (gene → value)` matrix used by both deconvolution methods.
- `ImmuneAnalyzer.LoadBundledAbisSignatureMatrix()`: Loads the bundled **ABIS-Seq** signature matrix (Monaco et al., 2019, CC BY 4.0; 1296 genes × 17 immune cell types) — a real, published, permissively-licensed reference matrix that makes ν-SVR deconvolution work out-of-the-box (Section 6.2).

### 5.2 Current Behavior

The repository currently ships the full 141-gene ESTIMATE immune and stromal signatures as defaults. For deconvolution, the *default* signature matrix (used when `signatureMatrix` is `null`) is a compact 22-cell-type matrix with representative marker genes rather than the full LM22. In addition, a real published reference matrix now ships: the **ABIS-Seq** matrix of Monaco et al. (2019) (1296 genes × 17 immune cell types, CC BY 4.0) is bundled as an embedded resource and loaded via `LoadBundledAbisSignatureMatrix()`, so callers can run `DeconvoluteImmuneCellsNuSvr(profile, ImmuneAnalyzer.LoadBundledAbisSignatureMatrix())` out-of-the-box. The NNLS solver normalizes fitted fractions to sum to `1` and reports both correlation and RMSE against the reconstructed expression profile. The ν-SVR path adds a faithful linear-kernel ν-SVR solver (Schölkopf et al., 2000) verified against scikit-learn `NuSVR` (libsvm) and by planted-truth recovery (on both the synthetic matrix and the bundled ABIS matrix); it accepts the bundled ABIS matrix, the compact default, or a caller-supplied LM22 via `LoadSignatureMatrix`. All entry points are configurable. The *defaults* of `EstimateInfiltration`, `DeconvoluteImmuneCells`, and `DeconvoluteImmuneCellsNuSvr` are unchanged — bundling ABIS is purely additive.

### 5.3 Conformance to Theory / Spec

#### 5.3.A ESTIMATE-style infiltration scoring

**Implemented (verbatim from the cited theory/spec):**

- Rank-based ssGSEA-style enrichment scoring for immune and stromal signatures.
- ESTIMATE score as the sum of immune and stromal enrichment scores.
- Tumor-purity calculation with the Yoshihara et al. cosine formula `purity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`.
- **Opt-in absolute purity** via `EstimateTumorPurity(double)`: applies the same Yoshihara cosine to a caller-supplied, Affymetrix/ABSOLUTE-calibrated ESTIMATE score, and — mirroring the ESTIMATE/`tidyestimate` reference implementation (`purity = ifelse(purity < 0, NA, purity)`) — returns `NaN` when the cosine is negative (out of the calibrated domain) instead of clamping.

**Intentionally simplified:**

- The default `EstimateInfiltration.TumorPurity` field applies the Yoshihara cosine to the library's single-sample **un-normalised** ssGSEA integral, which is on a different numeric scale than the cohort-/rank-normalised ESTIMATE score the coefficients were calibrated for; **consequence:** that default field is a **relative** indicator. For an **absolute** Affymetrix-calibrated purity, callers supply a true ESTIMATE score to the opt-in `EstimateTumorPurity`.

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

- (none) — ν-SVR-based CIBERSORT fitting is now provided by `DeconvoluteImmuneCellsNuSvr` (Section 5.3.C).

#### 5.3.C CIBERSORT ν-SVR deconvolution

**Implemented (verbatim from the cited theory/spec):**

- The linear-kernel ν-SVR dual of Schölkopf et al. (2000) (eqs 60–62): objective `−½ Σ(α_i−α_i*)(α_j−α_j*)⟨x_i,x_j⟩ + Σ y_i(α_i−α_i*)` under `Σ(α_i−α_i*)=0`, `Σ(α_i+α_i*) ≤ Cνℓ`, `α_i,α_i* ∈ [0,C]`, with `w = Σ(α_i−α_i*)x_i`.
- CIBERSORT's ν sweep {0.25, 0.5, 0.75}, lowest-RMSE selection, z-score standardisation, zero-clip of negative weights, and sum-to-1 normalisation (Newman et al., 2015; Chen et al., 2018).
- An LM22-format TSV loader (`LoadSignatureMatrix`) for the caller-supplied 547-gene × 22-cell-type matrix.
- A **bundled, permissively-licensed reference matrix**: the ABIS-Seq matrix (Monaco et al., 2019, CC BY 4.0; 1296 genes × 17 cell types) ships as an embedded resource and is loaded via `LoadBundledAbisSignatureMatrix()`, so the ν-SVR engine deconvolutes out-of-the-box. Planted-truth recovery on this matrix is verified (well-separated lineages within tolerance; a pure single population recovers exactly).

**Intentionally simplified:**

- The matrix used when `signatureMatrix` is `null` remains the representative 5-marker × 22-cell-type matrix, not LM22; **consequence:** the *default* ν-SVR fractions are illustrative — for a real workflow pass `LoadBundledAbisSignatureMatrix()` (bundled ABIS) or a caller-supplied LM22.

**Not implemented:**

- The CIBERSORT **LM22** matrix specifically (caller-supplied for licence reasons — Stanford no-redistribution, Section 6.2) and exact reproduction of the CIBERSORT tool's published per-sample fractions; **users should rely on:** the bundled ABIS matrix for out-of-the-box deconvolution, or the official CIBERSORT/CIBERSORTx tool when bit-exact parity with its full pipeline (LM22 + quantile normalisation + permutation p-value) is required.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Representative-marker *default* matrix instead of full LM22 | Assumption | The `null`-default fractions are less detailed than a full reference-matrix workflow | accepted | A real reference matrix now ships: `LoadBundledAbisSignatureMatrix()` (Monaco 2019, CC BY 4.0); or supply a custom matrix / LM22 via `LoadSignatureMatrix` |
| 2 | **LM22 specifically** not bundled (Stanford no-redistribution licence) | Deviation | CIBERSORT's canonical LM22 signature is caller-supplied, not embedded | accepted | Stanford licence forbids redistribution (Section 6.2); the permissive ABIS matrix is bundled instead, plus an LM22 loader |
| 3 | No bit-exact CIBERSORT-tool parity | Deviation | Per-sample fractions may differ from the official tool's output | accepted | Verified instead by planted-truth recovery + scikit-learn `NuSVR` cross-check |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty expression profile for `EstimateInfiltration` | Zero immune/stromal/ESTIMATE scores and tumor purity `≈ 0.8225` | The method evaluates the published purity formula at score `0` |
| No overlapping ESTIMATE signature genes | Corresponding enrichment score is `0` | The ssGSEA helper returns zero when the hit set is empty |
| No overlapping genes for deconvolution | All returned cell fractions are `0`, `Correlation = 0`, `Rmse = 0` | The method exits through a no-overlap branch |
| No overlapping genes for ν-SVR | All fractions `0`, `BestNu = 0`, `Correlation = 0`, `Rmse = 0` | Same no-overlap branch as NNLS |
| Non-finite (NaN/±Infinity) expression value for ν-SVR | `ArgumentException` | The linear-mixture model is defined only for finite expression; non-finite input is rejected up front so no NaN/Infinity leaks into the contracted-finite `CellFractions`/`Correlation`/`Rmse` |
| Near-constant / zero-variance mixture or reconstruction for ν-SVR | `Correlation = 0` (undefined), other fields finite | The Pearson helper clamps each variance term to its `≥ 0` floor (catastrophic cancellation can drive it slightly negative) and returns `0` for a non-positive or non-finite variance, so the correlation never becomes `NaN` |
| LM22-format TSV malformed (empty / ragged / non-numeric) | `FormatException` | The loader validates the header and every row |

### 6.2 Limitations

The deconvolution path is only as strong as its signature matrix. The built-in ESTIMATE signatures are complete; the `null`-default deconvolution matrix is intentionally compact, but a real published reference matrix now ships out-of-the-box.

**Bundled ABIS-Seq matrix (Monaco 2019, CC BY 4.0).** The ABIS-Seq immune signature matrix of Monaco et al. (2019) — 1296 genes × 17 immune cell types — is bundled as an embedded resource ([Resources/ABIS_sigmatrixRNAseq.tsv](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/Resources/ABIS_sigmatrixRNAseq.tsv)) and exposed via `LoadBundledAbisSignatureMatrix()`. The matrix is Table S5 (sheet "ABIS-Seq") of the open-access Cell Reports article, which is published under the Creative Commons Attribution 4.0 (CC BY 4.0) licence — *"© 2019 The Authors. This is an open access article under the CC BY license"* — so it may be redistributed with attribution. Attribution: Monaco G, Lee B, Xu W, et al., *Cell Reports* 26(6):1627–1640.e7, 2019, doi:10.1016/j.celrep.2019.01.041 (the provenance/licence is also recorded in the resource file header). This makes `DeconvoluteImmuneCellsNuSvr` work out-of-the-box.

**LM22 licence (CIBERSORT).** The canonical CIBERSORT signature matrix LM22 (547 genes × 22 cell types) is distributed by Stanford under a non-commercial licence that explicitly forbids redistribution — *"RECIPIENT shall not distribute the Program or transfer it to any other person or organization without prior written permission from STANFORD"* — and is gated behind registration at https://cibersort.stanford.edu. It is therefore **not bundled** in this library (unlike the CC BY 4.0 ABIS matrix above). Callers obtain `LM22.txt` under their own CIBERSORT licence and load it with `LoadSignatureMatrix`. The ν-SVR algorithm itself is fully implemented and verified (planted-truth recovery on the bundled ABIS matrix and a scikit-learn/libsvm `NuSVR` cross-check). Bit-exact parity with the official CIBERSORT tool's published per-sample fractions is **not** claimed: that additionally requires LM22 and the tool's full quantile-normalisation/permutation pipeline.

The opt-in `EstimateTumorPurity` transform is calibrated, per Yoshihara et al. (2013), against ABSOLUTE purity on TCGA **Affymetrix** data by nonlinear least squares; it is valid for Affymetrix-derived ESTIMATE scores and should not be applied to RNA-seq-derived scores. It is the caller's responsibility to pass a true ESTIMATE-scale score (e.g. from the ESTIMATE R package); applying it to this library's single-sample un-normalised ssGSEA integral would reproduce only the relative `EstimateInfiltration.TumorPurity` value.

## 8. References

1. Yoshihara, K., et al. 2013. Inferring tumour purity and stromal and immune cell admixture from expression data. Nature Communications 4:2612. https://doi.org/10.1038/ncomms3612.
2. Barbie, D. A., et al. 2009. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature 462:108-112. https://doi.org/10.1038/nature08460.
3. Hänzelmann, S., R. Castelo, and J. Guinney. 2013. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics 14:7.
4. Lawson, C. L., and R. J. Hanson. 1995. Solving Least Squares Problems. SIAM Classics in Applied Mathematics.
5. Abbas, A. R., et al. 2009. Deconvolution of blood microarray data identifies cellular activation patterns in systemic lupus erythematosus. PLoS One 4:e6098.
6. Newman, A. M., et al. 2015. Robust enumeration of cell subsets from tissue expression profiles. Nature Methods 12(5):453-457. https://doi.org/10.1038/nmeth.3337.
7. Newman, A. M., et al. 2019. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nature Biotechnology 37(7):773-782. https://doi.org/10.1038/s41587-019-0114-2.
8. Subramanian, A., et al. 2005. Gene set enrichment analysis. Proceedings of the National Academy of Sciences 102(43):15545-15550. https://doi.org/10.1073/pnas.0506580102.
9. Schölkopf, B., A. J. Smola, R. C. Williamson, and P. L. Bartlett. 2000. New support vector algorithms. Neural Computation 12(5):1207-1245. https://doi.org/10.1162/089976600300015565. (ν-SVR dual: Smola & Schölkopf 2004 tutorial eqs 60–62, https://alex.smola.org/papers/2003/SmoSch03b.pdf.)
10. Chen, B., et al. 2018. Profiling Tumor Infiltrating Immune Cells with CIBERSORT. Methods in Molecular Biology 1711:243-259. https://pmc.ncbi.nlm.nih.gov/articles/PMC5895181/.
11. CIBERSORT licence (Stanford University). Non-commercial, no-redistribution terms; registration at https://cibersort.stanford.edu. Verbatim clauses: https://gist.github.com/dhimmel/58dcd9b512e669f20a65ddf73997b733.
12. Monaco, G., B. Lee, W. Xu, et al. 2019. RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types. Cell Reports 26(6):1627-1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041. (Open access, CC BY 4.0; PMC6367568. Bundled ABIS-Seq matrix = Table S5, sheet "ABIS-Seq", `mmc6.xlsx`.)

