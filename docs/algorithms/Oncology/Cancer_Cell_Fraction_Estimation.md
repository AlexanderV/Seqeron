# Cancer Cell Fraction Estimation and CCF Clustering

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / tumor subclonal reconstruction |
| Test Unit ID | ONCO-CCF-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Cancer cell fraction (CCF) is the fraction of cancer cells in a sequenced sample that carry a given somatic
mutation. It is not observed directly; this unit computes the standard point estimate of CCF from a mutation's
variant allele fraction (VAF), the tumor purity, the local tumor copy number, and the mutation multiplicity, then
clusters a set of CCF values into clones/subclones and identifies the clonal cluster. `EstimateCcf` is an exact
closed-form (deterministic) estimate; `ClusterCcfValues` is the globally optimal one-dimensional k-means partition
(Ckmeans.1d.dp dynamic programming, Wang & Song 2011 [6]).
CCF estimation underpins clonal/subclonal classification (ONCO-CLONAL-001) and tumor phylogeny (ONCO-PHYLO-001).

## 2. Scientific / Formal Basis

> A = CCF point estimation, B = 1D CCF clustering

### 2.A CCF point estimation

#### Domain Context

A heterozygous SNV in a pure, diploid region has VAF ≈ ½·CCF. Tumor purity (admixed normal DNA), copy-number
aberrations, and the mutation multiplicity (number of mutated copies per cancer cell) all distort the VAF→CCF
relationship, so CCF must be inferred from all four quantities [1][3].

#### Core Model

Per McGranahan et al. (2016) the observed mutation copy number is
`n_mut = VAF·(1/ρ)·[ρ·CN_t + CN_n·(1−ρ)]`, and CCF = n_mut / m [3]. With the normal locus diploid (CN_n = 2)
this gives the standard estimate [1][2][3]:

```
CCF = VAF · (ρ·N_T + 2(1−ρ)) / (ρ · m)
```

where VAF is the variant allele fraction, ρ the tumor purity, N_T the local tumor copy number, and m the integer
mutation multiplicity. Equivalently Zheng et al. (2022) define VAF = m·CCF·ρ / (ρ·N_T + 2(1−ρ)), which inverts to
the same expression [2]. Multiplicity itself can be estimated as `m = VAF·(ρ·N_T + 2(1−ρ))/ρ` rounded to the
nearest non-zero integer for clonal copy-number regions [1]; here m is supplied by the caller.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-CCF-01 | Local copy-number state (N_T) and multiplicity (m) are correct for the locus | Biased CCF at mis-segmented/amplified/LOH loci |
| ASM-CCF-02 | Normal admixture is diploid (CN_n = 2) | Mis-scaled denominator if the matched normal is aneuploid |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-CCF-01 | Reported CCF ∈ [0, 1] | Raw value capped at 1; a mutation in all cancer cells has CCF = 1 [3]; registry invariant |
| INV-CCF-02 | CCF strictly increases with VAF (other inputs fixed) | Formula is linear in VAF with positive slope [1] |
| INV-CCF-03 | CCF = 0 ⇔ VAF = 0 | Numerator is VAF·(positive constant) |

### 2.B 1D CCF clustering

#### Domain Context

Mutations sharing a CCF belong to the same clonal population; grouping CCF values recovers the clones/subclones.
The cluster with the highest CCF is the clonal cluster, the rest are subclonal lineages [1].

#### Core Model

The k-means objective [5] is the partition of the values into k clusters minimizing the within-cluster sum of
squares `Σ_j Σ_{x∈S_j} (x − μ_j)²`. Lloyd's alternating assignment/update iteration [5] only reaches a *local*
optimum that depends on seeding — also in one dimension (B24 review 2026-09, F17: with quantile seeding it was
suboptimal on 2367/4999 random 1-D inputs). In one dimension an optimal partition consists of contiguous blocks of
the sorted values, so the global optimum is found exactly by dynamic programming: with D[q][i] the minimum WCSS of
x₁..xᵢ in q clusters, `D[q][i] = min_{q≤j≤i} D[q−1][j−1] + ssq(x_j..x_i)` (Wang & Song 2011, Ckmeans.1d.dp [6]),
followed by backtracking. The result needs no seeding and is fully reproducible.

When k is not known, Ckmeans.1d.dp chooses it in [k_min, k_max] by the Bayesian information criterion
(`select_levels`, R default `estimate.k = "BIC"` [6]): for each K the optimal partition defines a Gaussian mixture
with weights λ_k = n_k/N, block means μ_k and unbiased block variances σ²_k (σ² = 0 → d²/36, singleton → d², d the
gap to the nearest neighbouring value outside the block); `BIC(K) = 2·Σ_i ln Σ_k λ_k·φ(x_i; μ_k, σ²_k) − (3K − 1)·ln N`,
and the smallest K with the maximal BIC is selected.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-CL-01 | The number of clones k is supplied, or chosen in [k_min, k_max] by the Ckmeans.1d.dp Gaussian-mixture BIC | Wrong k / non-Gaussian clusters over-/under-split populations |
| ASM-CL-02 | Clones are separable by 1D CCF distance | Overlapping CCF distributions merge distinct clones |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-CL-01 | Every value assigned to exactly one cluster in [0, k) | k-means produces a partition [5] |
| INV-CL-02 | Each centroid equals the mean of its members | Update step [5] |
| INV-CL-03 | Output is deterministic for a given (values, k) | Exact DP on sorted values; no seeding/RNG |
| INV-CL-05 | WCSS is the global minimum over all partitions into min(k, distinct) clusters | Ckmeans.1d.dp optimality [6] |
| INV-CL-06 | Every returned cluster is non-empty; #clusters = min(k, number of distinct values) | Ckmeans.1d.dp `Kmax = min(k, nUnique)` [6] |
| INV-CL-04 | Clonal cluster index = argmax centroid (= k−1, centroids ascending) | "cluster with the highest CP … deemed clonal" [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| vaf | double | required | Variant allele fraction | [0, 1] |
| purity | double | required | Tumor purity ρ | (0, 1] |
| tumorCopyNumber | int | required | Local tumor copy number N_T | ≥ 1 |
| multiplicity | int | required | Mutated copies per cancer cell m | [1, tumorCopyNumber] |
| ccfValues | IReadOnlyList&lt;double&gt; | required | CCF values to cluster | non-empty, finite |
| clusterCount | int | required (fixed-k overload) | Number of clusters k | [1, count] |
| minClusters, maxClusters | int | required (BIC overload) | Range of k searched by BIC (R `k = c(kmin, kmax)`) | 1 ≤ minClusters ≤ maxClusters; capped at #distinct values as in R |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| CcfEstimate.Ccf | double | CCF capped to [0, 1] |
| CcfEstimate.RawCcf | double | Uncapped formula value (may exceed 1) |
| CcfClustering.Centroids | IReadOnlyList&lt;double&gt; | Cluster centroids, ascending |
| CcfClustering.Assignments | IReadOnlyList&lt;int&gt; | Per-value cluster index (input order) |
| CcfClustering.ClonalClusterIndex | int | Index of the highest-centroid (clonal) cluster |

### 3.3 Preconditions and Validation

`EstimateCcf` throws `ArgumentOutOfRangeException` for vaf ∉ [0,1], purity ∉ (0,1], or tumorCopyNumber < 1, and
`ArgumentException` for multiplicity ∉ [1, tumorCopyNumber]. `ClusterCcfValues` throws `ArgumentNullException`
for null input, `ArgumentException` for an empty list or a NaN/infinite value, and `ArgumentOutOfRangeException`
for clusterCount ∉ [1, count]; the BIC overload throws `ArgumentOutOfRangeException` for minClusters < 1 or
maxClusters < minClusters (maxClusters above the number of values is allowed: R caps it at the number of distinct
values; if minClusters exceeds that number both bounds become it). Either overload rejects a DP matrix above 10⁸
cells. All indices are 0-based.

## 4. Algorithm

### 4.A CCF point estimation

#### High-Level Steps

1. Validate inputs.
2. Compute total DNA per cell `D = ρ·N_T + 2(1−ρ)`.
3. Raw CCF = VAF·D / (ρ·m); reported CCF = min(1, raw).

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| EstimateCcf | O(1) | O(1) | closed form |

### 4.B 1D CCF clustering

#### High-Level Steps

1. Validate; stable-sort values carrying original indices; k ← min(k, number of distinct values) [6].
2. Median-shifted prefix sums Σx, Σx² (Ckmeans.1d.dp `EWL2::fill_dp_matrix`).
3. Fill the DP rows q = 1..k−1 with the O(n) SMAWK row fill (`fill_row_q_SMAWK`, R's default
   `method = "linear"`, incl. its tie-breaking among equal-WCSS splits); only two rows of D are kept, the backtrack
   matrix J is k × n.
4. Backtrack cluster boundaries from J; centroid = block mean; clusters are ascending, the last is clonal.
5. BIC overload: adjust [k_min, k_max] to the number of distinct values (R `cluster.1d.dp`), fill one DP for
   K = k_max, compute BIC(K) for K = k_min..k_max by `select_levels`, backtrack with the selected K (all values
   equal → one cluster).

#### Decision Rules / Reference Tables

Port of Ckmeans.1d.dp 4.3.6 C++ (`EWL2_dynamic_prog.cpp`, `EWL2_fill_SMAWK.cpp`, `EWL2_within_cluster.h`,
`dynamic_prog.cpp::backtrack`); `ldouble` is `double` there, so arithmetic is identical. The R defaults are
reproduced: `method = "linear"` (SMAWK), unweighted (`y = 1`), criterion `"L2"`. Inputs whose effective
k × n exceeds 10⁸ backtrack cells are rejected (`ArgumentOutOfRangeException`).

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ClusterCcfValues | O(k·n) | O(k·n) | SMAWK DP row fill [6]; sort O(n log n) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.Clonality.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.Clonality.cs)

- `OncologyAnalyzer.EstimateCcf(vaf, purity, tumorCopyNumber, multiplicity)`: point CCF estimate (capped + raw).
- `OncologyAnalyzer.ClusterCcfValues(ccfValues, clusterCount)`: optimal 1D k-means (Ckmeans.1d.dp) + clonal-cluster id.
- `OncologyAnalyzer.ClusterCcfValues(ccfValues, minClusters, maxClusters)`: same, with k chosen by BIC
  (R `Ckmeans.1d.dp(x, k = c(minClusters, maxClusters))`); chosen k = `Centroids.Count`.

### 5.2 Current Behavior

Multiplicity is an integer input of `EstimateCcf` (multi-region/PICTograph convention); it can be inferred from VAF,
purity and local copy number with `OncologyAnalyzer.DeriveMultiplicity` (ONCO-ASCAT-001; facets-suite
`expected_mutant_copies`, McGranahan 2016). Clustering is the exact Ckmeans.1d.dp optimum (no seeding), so output is identical
across runs and independent of input order. Cross-checked against R Ckmeans.1d.dp 4.3.6 `Ckmeans.1d.dp(x, k)`
(defaults): 6000/6000 random inputs (n 2–60, k 1–8, 5 of 6 generators tie-heavy) bit-identical labels and centroids
(F46; the former log-linear fill, = R `method = "loglinear"`, differed on 21 equal-WCSS ties). The BIC overload matches
R `Ckmeans.1d.dp(x, c(kmin, kmax))` on 3000/3000 random inputs (kmin 1–3, kmax kmin..kmin+9; selected k 1–10):
labels, centers and the `$BIC` vector bit-identical (F47). No substring/pattern search is
involved, so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

#### 5.3.A CCF point estimation

**Implemented (verbatim from the cited theory/spec):**

- CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m) exactly as in McGranahan 2016, Tarabichi 2021 Box 1, and Zheng 2022 [1][2][3].

**Intentionally simplified:**

- Reported CCF capped at 1; **consequence:** raw values >1 from noise are clipped (the uncapped value is exposed as `RawCcf`).

**Not implemented:**

- Posterior/uncertainty modeling of CCF; **users should rely on:** ONCO-CLONAL-001 `ClassifyClonality` (Bayesian grid posterior) when probabilistic CCF is needed.
- Multiplicity inference inside `EstimateCcf`; **users should rely on:** `OncologyAnalyzer.DeriveMultiplicity(vaf, purity, totalCopyNumber, majorCopyNumber)` (ONCO-ASCAT-001, facets-suite `expected_mutant_copies`) to derive m, then pass it to `EstimateCcf`.

#### 5.3.B 1D CCF clustering

**Implemented (verbatim from the cited theory/spec):**

- Exact minimum-WCSS k-means partition by Ckmeans.1d.dp dynamic programming [5][6]; clonal cluster = highest centroid [1].

**Intentionally simplified:**

- ~~k must be supplied~~ — **Resolved (F47):** automatic k by the Ckmeans.1d.dp BIC (`select_levels`) via
  `ClusterCcfValues(ccfValues, minClusters, maxClusters)`; R-locked bit for bit (labels, centers, BIC vector).

**Not implemented:**

- Dirichlet-process / beta-binomial mixture clustering; **users should rely on:** external tools (PyClone, SciClone, DPClust) for full probabilistic subclone inference.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| VAF = 0 | CCF = 0 | linear formula (INV-CCF-03) |
| raw CCF > 1 | reported 1.0, RawCcf > 1 | INV-CCF-01 / CNAqc [4] |
| multi-copy locus (N_T>2, m>1) | uses supplied m | multiplicity definition [1] |
| k = 1 | single cluster at the global mean; clonal index 0 | trivial partition |
| BIC overload, all values equal | one cluster; BIC vector all 0 | Ckmeans.1d.dp `nUnique == 1` branch [6] |
| BIC overload, minClusters > distinct values | both bounds set to the number of distinct values | R `cluster.1d.dp` [6] |
| k > distinct values | k reduced to the number of distinct values (no empty clusters) | Ckmeans.1d.dp [6] |
| empty values / null / k out of range | exception | validation |

### 6.2 Limitations

Point CCF carries no uncertainty; clustering assumes clones are 1D-separable and that k is known. Aneuploid
matched normals violate the CN_n = 2 assumption. Not a replacement for full probabilistic subclone callers.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var est = OncologyAnalyzer.EstimateCcf(vaf: 0.20, purity: 0.80, tumorCopyNumber: 2, multiplicity: 1);
// est.Ccf == 0.5  (0.20 · (0.8·2 + 2·0.2) / (0.8·1) = 0.20·2.0/0.8)

var clustering = OncologyAnalyzer.ClusterCcfValues(
    new[] { 1.0, 0.98, 0.96, 0.50, 0.48, 0.52 }, clusterCount: 2);
// Centroids ≈ {0.50, 0.98}; ClonalClusterIndex == 1
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_EstimateCcf_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_EstimateCcf_Tests.cs) — covers `INV-CCF-01..03`, `INV-CL-01..04`
- Evidence: [ONCO-CCF-001-Evidence.md](../../../docs/Evidence/ONCO-CCF-001-Evidence.md)
- Related algorithms: [Clonal_Subclonal_Classification](Clonal_Subclonal_Classification.md), [Tumor_Phylogeny_Reconstruction](Tumor_Phylogeny_Reconstruction.md)

## 8. References

1. Tarabichi M, Salcedo A, Deshwar AG, et al. 2021. A practical guide to cancer subclonal reconstruction from DNA sequencing. Nature Methods 18:144–155. https://pmc.ncbi.nlm.nih.gov/articles/PMC7867630/
2. Zheng J, et al. 2022. Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors. Bioinformatics 38(15):3677–3683. https://academic.oup.com/bioinformatics/article/38/15/3677/6596597
3. McGranahan N, Furness AJS, Rosenthal R, et al. 2016. Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. Science 351(6280):1463–1469. https://www.science.org/doi/10.1126/science.aaf1490
4. Caravagna G, et al. CNAqc — Computation of Cancer Cell Fractions. https://caravagnalab.github.io/CNAqc/articles/a4_ccf_computation.html
5. Lloyd SP. 1982. Least squares quantization in PCM. IEEE Trans. Inf. Theory 28(2):129–137. https://doi.org/10.1109/TIT.1982.1056489
6. Wang H, Song M. 2011. Ckmeans.1d.dp: Optimal k-means clustering in one dimension by dynamic programming. The R Journal 3(2):29–33. Reference code: Ckmeans.1d.dp 4.3.6 (github.com/cran/Ckmeans.1d.dp, `src/EWL2_*`, `R/Ckmeans.1d.dp.R`).
