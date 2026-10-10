# Tumor Heterogeneity Analysis (MATH, Shannon Diversity, Subclone Count)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-HETERO-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-10-10 (review-2026-09 B24: F20–F22, F56–F57) |

## 1. Overview

Quantifies intratumour heterogeneity (ITH) from somatic-mutation evidence using established, retrievable metrics. The MATH (Mutant-Allele Tumour Heterogeneity) score measures the spread of the variant-allele-fraction (VAF) distribution as `100·MAD/median` [1][2]; the Shannon diversity index `H = −Σ pᵢ ln pᵢ` measures clonal diversity over CCF-cluster fractions [4][5]; the subclone count is the number of occupied CCF clusters [4]; and the subclonal fraction is the proportion of mutations that are not clonal (clonal ⇔ CCF > 0.95) [6]. All metrics are exact deterministic statistics over the inputs.

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

Martinez et al. [4] compute the reference Shannon diversity of each clone mixture "using the clonal frequencies of each mixture", i.e. `pᵢ` = cellular frequency of clone `i` (prevalence-weighted); richness `R` = number of clones present [4]. `CalculateCloneShannonDiversity(cloneFractions)` implements exactly this (`pᵢ = fᵢ/Σf`, = `scipy.stats.entropy`). `AnalyzeHeterogeneity` — which receives per-mutation CCFs, not clone frequencies — uses `pᵢ` = fraction of mutations assigned to CCF cluster `i` (count-weighted; see §5.3).

**maftools `math.score` pre-filters** [3] (optional overload `CalculateITH(vafs, vafCutOff, minMutations = 5)`): keep VAFs with `!(vaf < vafCutOff)` (default `0.075`; a VAF equal to the cutoff is kept); if fewer than 5 remain (`length(vaf) < 5`, hard-coded) the sample is skipped (no row → `null`); otherwise MATH over the retained VAFs. maftools derives a missing `t_vaf` as `t_alt_count/(t_ref_count + t_alt_count)` and divides all VAFs by 100 when their maximum exceeds 1; here VAFs must already be fractions in [0, 1].

**Subclonal fraction** [6]: Landau et al. "classified a mutation as clonal if the CCF harboring it was >0.95 [with probability >0.5] and subclonal otherwise"; for CCF point estimates a mutation is subclonal when `CCF ≤ 0.95` (so exactly 0.95 is subclonal), and the fraction is `#(CCF ≤ 0.95)/n = 1 − |IdentifyClonalMutations|/n`.

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
| vafCutOff | double | required in the filtered overload (maftools 0.075 = `MaftoolsMathVafCutOff`) | drop VAF < cutoff | finite ∈ [0,1] |
| minMutations | int | 5 (`MaftoolsMathMinMutations`) | skip sample (null) if fewer retained | ≥ 1 |
| cloneFractions | IReadOnlyList\<double\> | required | cellular frequency of each clone (`CalculateCloneShannonDiversity`) | each finite ∈ [0,1]; sum > 0 (normalised) |
| ccfClusters | CcfClustering | required | clustering from `ClusterCcfValues` | ≥ 1 centroid and 1 assignment |
| variantAlleleFractions | IReadOnlyList\<double\> | required | per-mutation VAFs | each ∈ [0,1]; count = ccfValues count |
| ccfValues | IReadOnlyList\<double\> | required | per-mutation CCFs | each ∈ [0,1] |
| clusterCount | int | required | number of clones k | [1, count] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CalculateITH` return | double | MATH score (≥ 0) |
| `CalculateITH(vafs, vafCutOff, minMutations)` return | double? | MATH over VAFs ≥ cutoff; null when < minMutations retained (maftools skip) |
| `CalculateCloneShannonDiversity` return | double | H over clone cellular frequencies (≥ 0) |
| `InferSubclones` return | int | number of occupied CCF clusters |
| HeterogeneityResult.MathScore | double | MATH over VAFs |
| HeterogeneityResult.ShannonDiversity | double | H = −Σ pᵢ ln pᵢ |
| HeterogeneityResult.SubcloneCount | int | occupied clusters |
| HeterogeneityResult.SubclonalFraction | double | fraction with CCF ≤ 0.95 (not clonal) |

### 3.3 Preconditions and Validation

Null lists throw `ArgumentNullException`. Empty lists, non-finite or out-of-[0,1] values, mismatched VAF/CCF lengths, a `CcfClustering` assignment label outside `[0, centroid count)`, and a zero median (MATH division by zero) throw `ArgumentException`; `clusterCount` outside `[1, count]` throws `ArgumentOutOfRangeException`. Inputs are not mutated.

## 4. Algorithm

### 4.1 High-Level Steps

1. **MATH:** compute `median(f)`; reject median = 0; compute raw MAD = `median(|fᵢ − median|)`; `MATH = ((MAD·100)·1.4826)/median` — maftools' operation order, bit-identical to R (0/19 999 random vectors differ; the former `100·(1.4826·MAD)/median` order differed by 1 ulp on 6 912).
2. **Subclones:** cluster CCFs (ONCO-CCF-001); count distinct occupied cluster labels.
3. **Shannon:** clone fractions `pᵢ = sizeᵢ/n`; `H = −Σ pᵢ ln pᵢ` via canonical `StatisticsHelper.ShannonIndex` (= `scipy.stats.entropy`).
4. **Subclonal fraction:** `1 − |IdentifyClonalMutations(ccf)|/n` (canonical Landau rule, CCF > 0.95 clonal).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- MAD consistency constant `1.4826 = 1/Φ⁻¹(3/4)` [2][3].
- Percentage scale `100` [1].
- Clonal CCF threshold `0.95` (reused from ONCO-CLONAL-001, Landau et al. 2013) [6].
- Median = canonical `StatisticsHelper.Median` (R `median`: mean of the two central order statistics for even counts) [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateITH | O(n log n) | O(n) | two median sorts |
| InferSubclones | O(n) | O(k) | hash set of labels |
| AnalyzeHeterogeneity | O(n log n) | O(n) | MATH + clustering (ONCO-CCF-001) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.PhylogenyHeterogeneity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.PhylogenyHeterogeneity.cs); helpers `StatisticsHelper.Median` / `StatisticsHelper.ShannonIndex` in [StatisticsHelper.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/StatisticsHelper.cs)

- `OncologyAnalyzer.CalculateITH(IReadOnlyList<double>)`: MATH score (no filtering).
- `OncologyAnalyzer.CalculateITH(IReadOnlyList<double>, double vafCutOff, int minMutations = 5)`: maftools `math.score` per-sample pre-filters, then MATH; `null` = skipped.
- `OncologyAnalyzer.CalculateCloneShannonDiversity(IReadOnlyList<double>)`: Shannon over clone cellular frequencies (Martinez et al. 2017) via canonical `StatisticsHelper.ShannonIndexOfWeights`.
- `OncologyAnalyzer.InferSubclones(CcfClustering)`: occupied-cluster count.
- `OncologyAnalyzer.AnalyzeHeterogeneity(IReadOnlyList<double>, IReadOnlyList<double>, int)`: aggregate ITH metrics.

### 5.2 Current Behavior

Reuses `ClusterCcfValues` (ONCO-CCF-001, Ckmeans.1d.dp) for CCF clustering and `IdentifyClonalMutations` (ONCO-CLONAL-001) for the subclonal cutoff. No search/matching is involved, so the repository suffix tree is not applicable. The median helper clones the input before sorting, so callers' arrays are never mutated.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `MATH = 100·1.4826·median(|f − median|)/median` exactly as Mroz & Rocco (2013) / Mroz et al. (2015) and maftools `mathScore.R` [1][2][3].
- `H = −Σ pᵢ ln pᵢ` with natural log, over clone fractions [4][5]; `CalculateCloneShannonDiversity` takes the clones' cellular frequencies exactly as Martinez et al. [4] (F57).
- maftools `math.score` pre-filters (`vafCutOff` = 0.075 kept-if-equal, skip if < 5 retained) via the `CalculateITH(vafs, vafCutOff, minMutations)` overload, bit-identical to R `math.score` (F56) [3]. The single-argument `CalculateITH` stays the unfiltered MATH definition [1][2]; the filters are maftools' caller-side wrapper.
- Subclonal ⇔ not clonal, clonal ⇔ CCF > 0.95 [6].

**Intentionally simplified:**

- `AnalyzeHeterogeneity.ShannonDiversity` uses per-cluster mutation counts (cluster sizes) because its inputs are per-mutation CCFs, which do not give clone cellular frequencies directly (nested subclones' CCFs overlap); **consequence:** H reflects mutation-count diversity. The source's prevalence-weighted index [4] is available as `CalculateCloneShannonDiversity(cloneFractions)` (F57).

**Not implemented:**

- Probabilistic/Bayesian subclone inference (e.g., PyClone/SciClone posterior clustering); **users should rely on:** dedicated tools — clustering here is the exact 1-D k-means (Ckmeans.1d.dp) of ONCO-CCF-001.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Even-count median = mean of central two | Assumption | matches maftools `median` | accepted | R convention [3] |
| 2 | Shannon over cluster sizes in `AnalyzeHeterogeneity` | Assumption | mutation-count diversity | accepted | see 5.3; source form = `CalculateCloneShannonDiversity` |
| 3 | `CalculateITH` filtered overload requires VAFs in [0,1] | Assumption | no maftools ×1/100 rescale of percent VAFs | accepted | maftools rescales only when max > 1 [3] |

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
| Filtered overload: VAF = vafCutOff | kept | maftools `!t_vaf < vafCutOff` [3] |
| Filtered overload: < minMutations retained | null (skipped) | maftools `length(vaf) < 5` [3] |

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

- Tests: [OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs) — covers INV-01..INV-05
- Evidence: [ONCO-HETERO-001-Evidence.md](../../../docs/Evidence/ONCO-HETERO-001-Evidence.md)
- Related algorithms: [Cancer_Cell_Fraction](../Oncology/Cancer_Cell_Fraction_Estimation.md)

## 8. References

1. Mroz EA, Rocco JW. 2013. MATH, a novel measure of intratumor genetic heterogeneity, is high in poor-outcome classes of head and neck squamous cell carcinoma. Oral Oncology 49(3):211–215. https://pubmed.ncbi.nlm.nih.gov/23079694/
2. Mroz EA, Tward AD, Hammon RJ, Ren Y, Rocco JW. 2015. Intra-tumor genetic heterogeneity and mortality in head and neck cancer. PLOS Medicine 12(2):e1001786. https://doi.org/10.1371/journal.pmed.1001786
3. Mayakonda A et al. maftools `mathScore.R`. https://github.com/PoisonAlien/maftools/blob/master/R/mathScore.R
4. Martinez P et al. 2017. Quantification of within-sample genetic heterogeneity from SNP-array data. Scientific Reports 7:3248. https://doi.org/10.1038/s41598-017-03496-0 (PMC5468233; earlier revisions mis-cited this as "Liu & Zhang, BMC Genomics 18:457").
5. Shannon CE. 1948. A mathematical theory of communication. Bell System Technical Journal 27:379–423. https://en.wikipedia.org/wiki/Diversity_index#Shannon_index
6. Landau DA et al. 2013. Evolution and impact of subclonal mutations in chronic lymphocytic leukemia. Cell 152(4):714–726. https://doi.org/10.1016/j.cell.2013.01.019
