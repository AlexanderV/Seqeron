# Tumor Heterogeneity Analysis (MATH, Shannon Diversity, Subclone Count)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-HETERO-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-15 |

## 1. Overview

Quantifies intratumour heterogeneity (ITH) from somatic-mutation evidence using established, retrievable metrics. The MATH (Mutant-Allele Tumour Heterogeneity) score measures the spread of the variant-allele-fraction (VAF) distribution as `100·MAD/median` [1][2]; the Shannon diversity index `H = −Σ pᵢ ln pᵢ` measures clonal diversity over CCF-cluster fractions [4][5]; the subclone count is the number of occupied CCF clusters [4]; and the subclonal fraction is the proportion of mutations with CCF below the clonal threshold [6]. All metrics are exact deterministic statistics over the inputs.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A tumour is a mixture of genetically distinct cell populations (clones/subclones). Bulk sequencing reports each somatic mutation's VAF; correcting VAF for purity and copy number yields the cancer cell fraction (CCF, ONCO-CCF-001). The dispersion of VAFs and the diversity of CCF clusters are surrogates for the degree of intratumour heterogeneity [1][4].

### 2.2 Core Model

**MATH score** [1][2][3]:

```
MATH = 100 × MAD / median(f)
MAD  = 1.4826 × median(|fᵢ − median(f)|)
```

where `fᵢ` are the mutant-allele fractions. The 1.4826 factor scales the median absolute deviation so that, for a normally distributed variable, its expected MAD equals its standard deviation [2]. maftools `mathScore.R` computes the identical quantity: `pat.math = (median(abs(vaf − median(vaf)))·100)·1.4826 / median(vaf)` [3].

**Shannon diversity** [4][5]:

```
H = −Σᵢ pᵢ · ln(pᵢ)      (natural logarithm)
```

where `pᵢ` is the fraction of mutations assigned to clone/cluster `i` and richness `R` is the number of occupied clusters [4].

**Subclonal fraction** [6]: a mutation is subclonal when `CCF < 0.95`; the fraction is `#(CCF < 0.95)/n`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | MATH ≥ 0 | MAD ≥ 0 and median > 0 (enforced) [1] |
| INV-02 | MATH = 0 ⇔ MAD = 0 (all VAFs equal the median) | numerator 0 [1] |
| INV-03 | H ≥ 0; H = 0 ⇔ one occupied clone | −p ln p ≥ 0 for p ∈ (0,1] [5] |
| INV-04 | H ≤ ln(richness), equality for equal clones | maximum entropy of a uniform distribution [5] |
| INV-05 | 1 ≤ subclone count ≤ k; 0 ≤ subclonal fraction ≤ 1 | counts over n mutations [4][6] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| ccfDistribution | IReadOnlyList\<double\> | required | mutant-allele fractions for MATH | each finite ∈ [0,1]; median > 0 |
| ccfClusters | CcfClustering | required | clustering from `ClusterCcfValues` | ≥ 1 centroid and 1 assignment |
| variantAlleleFractions | IReadOnlyList\<double\> | required | per-mutation VAFs | each ∈ [0,1]; count = ccfValues count |
| ccfValues | IReadOnlyList\<double\> | required | per-mutation CCFs | each ∈ [0,1] |
| clusterCount | int | required | number of clones k | [1, count] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CalculateITH` return | double | MATH score (≥ 0) |
| `InferSubclones` return | int | number of occupied CCF clusters |
| HeterogeneityResult.MathScore | double | MATH over VAFs |
| HeterogeneityResult.ShannonDiversity | double | H = −Σ pᵢ ln pᵢ |
| HeterogeneityResult.SubcloneCount | int | occupied clusters |
| HeterogeneityResult.SubclonalFraction | double | fraction with CCF < 0.95 |

### 3.3 Preconditions and Validation

Null lists throw `ArgumentNullException`. Empty lists, non-finite or out-of-[0,1] values, mismatched VAF/CCF lengths, and a zero median (MATH division by zero) throw `ArgumentException`; `clusterCount` outside `[1, count]` throws `ArgumentOutOfRangeException`. Inputs are not mutated.

## 4. Algorithm

### 4.1 High-Level Steps

1. **MATH:** compute `median(f)`; reject median = 0; compute raw MAD = `median(|fᵢ − median|)`; `MATH = 100·1.4826·MAD/median`.
2. **Subclones:** cluster CCFs (ONCO-CCF-001); count distinct occupied cluster labels.
3. **Shannon:** clone fractions `pᵢ = sizeᵢ/n`; `H = −Σ pᵢ ln pᵢ`.
4. **Subclonal fraction:** count CCF < 0.95, divide by n.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- MAD consistency constant `1.4826 = 1/Φ⁻¹(3/4)` [2][3].
- Percentage scale `100` [1].
- Clonal CCF threshold `0.95` (reused from ONCO-CLONAL-001, Landau et al. 2013) [6].
- Median for even counts = mean of the two central order statistics (R/maftools convention) [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateITH | O(n log n) | O(n) | two median sorts |
| InferSubclones | O(n) | O(k) | hash set of labels |
| AnalyzeHeterogeneity | O(n log n) | O(n) | MATH + clustering (ONCO-CCF-001) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CalculateITH(IReadOnlyList<double>)`: MATH score.
- `OncologyAnalyzer.InferSubclones(CcfClustering)`: occupied-cluster count.
- `OncologyAnalyzer.AnalyzeHeterogeneity(IReadOnlyList<double>, IReadOnlyList<double>, int)`: aggregate ITH metrics.

### 5.2 Current Behavior

Reuses `ClusterCcfValues` (ONCO-CCF-001) for CCF clustering and the existing `ClonalCcfThreshold` constant (ONCO-CLONAL-001) for the subclonal cutoff. No search/matching is involved, so the repository suffix tree is not applicable. The median helper clones the input before sorting, so callers' arrays are never mutated.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `MATH = 100·1.4826·median(|f − median|)/median` exactly as Mroz & Rocco (2013) / Mroz et al. (2015) and maftools `mathScore.R` [1][2][3].
- `H = −Σ pᵢ ln pᵢ` with natural log, over clone fractions [4][5].
- Subclonal cutoff CCF < 0.95 [6].

**Intentionally simplified:**

- Clone fractions for Shannon use per-cluster mutation counts (cluster sizes), not posterior cluster weights; **consequence:** H reflects mutation-count diversity rather than CCF-weighted diversity.

**Not implemented:**

- Probabilistic/Bayesian subclone inference (e.g., PyClone/SciClone posterior clustering); **users should rely on:** dedicated tools — clustering here is the deterministic k-means of ONCO-CCF-001.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Even-count median = mean of central two | Assumption | matches maftools `median` | accepted | R convention [3] |
| 2 | Shannon over cluster sizes | Assumption | mutation-count diversity | accepted | see 5.3 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| All VAFs identical | MATH = 0 | MAD = 0 [1] / INV-02 |
| Single VAF | MATH = 0 | median = value, MAD = 0 |
| Median VAF = 0 | ArgumentException | division by zero [3] |
| Single clone | H = 0 | INV-03 |
| k equal clones | H = ln k | INV-04 |
| Null / empty / out-of-range | exception | §3.3 |

### 6.2 Limitations

MATH is sensitive to mutation calling and lacks strong clinical predictive power across some datasets (Heng et al. 2018, Sci Rep). The metrics summarise dispersion/diversity, not the clonal tree (see ONCO-PHYLO-001 for phylogeny).

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through (MATH, odd count):**

VAFs {0.10,0.20,0.30,0.40,0.50}: median = 0.30; abs deviations {0.20,0.10,0,0.10,0.20}; raw MAD = 0.10; scaled MAD = 1.4826·0.10 = 0.14826; MATH = 100·0.14826/0.30 = **49.42** [1].

**API usage example:**

```csharp
double math = OncologyAnalyzer.CalculateITH(new[] { 0.10, 0.20, 0.30, 0.40, 0.50 }); // 49.42
var ccf = new[] { 0.20, 0.25, 0.95, 1.00 };
var result = OncologyAnalyzer.AnalyzeHeterogeneity(new[] { 0.1, 0.12, 0.45, 0.50 }, ccf, clusterCount: 2);
// result.SubcloneCount, result.ShannonDiversity, result.SubclonalFraction
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs) — covers INV-01..INV-05
- Evidence: [ONCO-HETERO-001-Evidence.md](../../../docs/Evidence/ONCO-HETERO-001-Evidence.md)
- Related algorithms: [Cancer_Cell_Fraction](../Oncology/Cancer_Cell_Fraction_Estimation.md)

## 8. References

1. Mroz EA, Rocco JW. 2013. MATH, a novel measure of intratumor genetic heterogeneity, is high in poor-outcome classes of head and neck squamous cell carcinoma. Oral Oncology 49(3):211–215. https://pubmed.ncbi.nlm.nih.gov/23079694/
2. Mroz EA, Tward AD, Hammon RJ, Ren Y, Rocco JW. 2015. Intra-tumor genetic heterogeneity and mortality in head and neck cancer. PLOS Medicine 12(2):e1001786. https://doi.org/10.1371/journal.pmed.1001786
3. Mayakonda A et al. maftools `mathScore.R`. https://github.com/PoisonAlien/maftools/blob/master/R/mathScore.R
4. Liu Z, Zhang S. 2017. Quantification of within-sample genetic heterogeneity from SNP-array data. BMC Genomics 18:457. https://pmc.ncbi.nlm.nih.gov/articles/PMC5468233/
5. Shannon CE. 1948. A mathematical theory of communication. Bell System Technical Journal 27:379–423. https://en.wikipedia.org/wiki/Diversity_index#Shannon_index
6. Landau DA et al. 2013. Evolution and impact of subclonal mutations in chronic lymphocytic leukemia. Cell 152(4):714–726. https://doi.org/10.1016/j.cell.2013.01.019
