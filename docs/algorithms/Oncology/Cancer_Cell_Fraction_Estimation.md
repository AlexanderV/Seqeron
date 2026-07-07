# Cancer Cell Fraction Estimation and CCF Clustering

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / tumor subclonal reconstruction |
| Test Unit ID | ONCO-CCF-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-15 |

## 1. Overview

Cancer cell fraction (CCF) is the fraction of cancer cells in a sequenced sample that carry a given somatic
mutation. It is not observed directly; this unit computes the standard point estimate of CCF from a mutation's
variant allele fraction (VAF), the tumor purity, the local tumor copy number, and the mutation multiplicity, then
clusters a set of CCF values into clones/subclones and identifies the clonal cluster. `EstimateCcf` is an exact
closed-form (deterministic) estimate; `ClusterCcfValues` is a deterministic one-dimensional k-means partition.
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

Lloyd's k-means [5] partitions values into k clusters minimizing the within-cluster sum of squares
`Σ_j Σ_{x∈S_j} (x − μ_j)²` by alternating an assignment step (each value to the nearest centroid by squared
distance) and an update step (each centroid = mean of its members) until assignments stabilize [5]. In one
dimension this is exact and, with deterministic seeding, fully reproducible.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-CL-01 | The number of clones k is supplied | Wrong k over-/under-splits populations |
| ASM-CL-02 | Clones are separable by 1D CCF distance | Overlapping CCF distributions merge distinct clones |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-CL-01 | Every value assigned to exactly one cluster in [0, k) | k-means produces a partition [5] |
| INV-CL-02 | Each centroid equals the mean of its members | Update step [5] |
| INV-CL-03 | Output is deterministic for a given (values, k) | Quantile seeding uses no RNG; sorted, order-independent |
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
| clusterCount | int | required | Number of clusters k | [1, count] |

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
for clusterCount ∉ [1, count]. All indices are 0-based.

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

1. Validate; sort values carrying original indices.
2. Seed k centroids at evenly-spaced quantiles of the sorted data (deterministic).
3. Iterate assignment + update steps until assignments stop changing.
4. Relabel clusters by ascending centroid; the last (highest) is clonal.

#### Decision Rules / Reference Tables

Centroid j is seeded at the sorted value at quantile (j + 0.5)/k. Ties in nearest-centroid go to the lower index.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ClusterCcfValues | O(n·k·i) | O(n+k) | i ≤ n+1 iterations; sort O(n log n) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.EstimateCcf(vaf, purity, tumorCopyNumber, multiplicity)`: point CCF estimate (capped + raw).
- `OncologyAnalyzer.ClusterCcfValues(ccfValues, clusterCount)`: deterministic 1D k-means + clonal-cluster id.

### 5.2 Current Behavior

Multiplicity is a caller-supplied integer (multi-region/PICTograph convention); automatic multiplicity inference
from VAF is out of scope. Clustering uses no random seeding — centroids are seeded at quantiles of the sorted
input — so output is identical across runs and independent of input order. No substring/pattern search is
involved, so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

#### 5.3.A CCF point estimation

**Implemented (verbatim from the cited theory/spec):**

- CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m) exactly as in McGranahan 2016, Tarabichi 2021 Box 1, and Zheng 2022 [1][2][3].

**Intentionally simplified:**

- Reported CCF capped at 1; **consequence:** raw values >1 from noise are clipped (the uncapped value is exposed as `RawCcf`).

**Not implemented:**

- Posterior/uncertainty modeling of CCF; **users should rely on:** ONCO-CLONAL-001 `ClassifyClonality` (Bayesian grid posterior) when probabilistic CCF is needed.
- Automatic multiplicity inference; **users should rely on:** caller-supplied integer multiplicity.

#### 5.3.B 1D CCF clustering

**Implemented (verbatim from the cited theory/spec):**

- Lloyd assignment/update steps minimizing WCSS [5]; clonal cluster = highest centroid [1].

**Intentionally simplified:**

- k must be supplied; **consequence:** no automatic model selection of the number of clones.

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
