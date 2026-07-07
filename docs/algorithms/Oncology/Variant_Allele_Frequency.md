# Variant Allele Frequency Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-VAF-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Variant allele frequency (VAF) analysis quantifies, at a single genomic locus, the fraction of sequencing reads that support the alternate allele, and expresses the statistical uncertainty of that estimate. This unit implements three exact, specification-driven computations: the empirical VAF (alt-supporting reads / total reads), a Wilson score confidence interval for the underlying allele proportion, and a purity/ploidy correction that recovers the per-tumour-copy mutant fraction from an observed VAF. All three are closed-form and deterministic; none are heuristic or probabilistic in the sense of a caller model.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In tumor sequencing, a somatic variant is observed as a mixture of reference- and alternate-supporting reads at a locus. The empirical allele fraction (alt AD / Σ AD in a VCF, or read counts from samtools mpileup) is the model-free estimate of how prevalent the allele is among observed reads [2]. It differs from a caller's modelled allele fraction such as Mutect2's `AF`, which is a Bayesian estimate marginalised over allele fractions rather than the raw ratio [2]. Because the VAF is a binomial proportion (alt reads out of n covering reads), its sampling uncertainty is summarised by a binomial proportion confidence interval [1]. Observed VAF is further distorted by non-tumour (normal) cell admixture and by copy-number changes; correcting for tumor purity and segment ploidy recovers the biologically meaningful mutant fraction [3][4].

### 2.2 Core Model

**Empirical VAF** [2]:

```
VAF = altReads / totalReads
```

**Wilson score interval** for the proportion p̂ = altReads / n (n = totalReads) at standard-normal quantile z [1]:

```
center = (p̂ + z²/(2n)) / (1 + z²/n)
margin = (z / (1 + z²/n)) · √( p̂(1−p̂)/n + z²/(4n²) )
interval = center ± margin
```

For a 95% interval, z = 1.96 [1]. The Wilson interval is asymmetric and is contained in [0, 1] (no overshoot) with non-zero width even at p̂ = 0 or 1, in contrast to the Wald interval p̂ ± (z/√n)·√(p̂(1−p̂)) which can fall outside [0, 1] [1].

**Expected VAF as a function of purity and copy number** [3][4]:

```
v = (m · π) / ( 2(1−π) + π · n_tot )
```

where m = mutation multiplicity, π = tumor purity, n_tot = tumor total copy number, and the normal contribution is fixed at 2 (autosomal diploid). For a clonal heterozygous SNV in a diploid region (m = 1, n_tot = 2) this reduces to v = π/2; e.g. 80% purity ⇒ expected VAF 0.4 [3].

**Purity/ploidy correction** is the inversion of the expected-VAF relation, recovering m·CCF [4]:

```
adjusted = vaf · ( 2(1−π) + π · n_tot ) / π
```

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ VAF ≤ 1 | altReads ≤ totalReads by precondition; ratio of non-negatives [2] |
| INV-02 | lower ≤ center ≤ upper | margin ≥ 0 by construction (square root of a non-negative term) [1] |
| INV-03 | 0 ≤ lower and upper ≤ 1 | Wilson interval is contained in [0, 1] (no overshoot) [1] |
| INV-04 | AdjustVAFForPurity(π/2, π, 2) = 1 | diploid heterozygous round-trip of the CNAqc relation [3][4] |

### 2.5 Comparison with Related Methods

| Aspect | Wilson score interval | Wald (normal approximation) |
|--------|------------------------|------------------------------|
| Bounds | always within [0, 1] | can overshoot below 0 / above 1 [1] |
| Width at p̂ = 0 or 1 | non-zero | zero (degenerate) [1] |
| Symmetry | asymmetric | symmetric about p̂ |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| altReads | int | required | Alt-supporting reads | ≥ 0 and ≤ totalReads |
| totalReads | int | required | Total covering reads | ≥ 0 (> 0 for a CI) |
| confidence | double | 0.95 | Two-sided confidence level | (0, 1); only 0.95 supported (z = 1.96) |
| vaf | double | required | Observed VAF | [0, 1] |
| purity | double | required | Tumor purity π | (0, 1] |
| ploidy | double | required | Tumor total copy number n_tot | > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| CalculateVAF | double | Empirical VAF in [0, 1] |
| VafConfidenceInterval.Vaf | double | VAF point estimate |
| VafConfidenceInterval.Lower / Upper | double | Wilson interval bounds in [0, 1] |
| VafConfidenceInterval.Confidence | double | Confidence level used |
| AdjustVAFForPurity | double | Purity/ploidy-corrected mutant fraction (m·CCF) |

### 3.3 Preconditions and Validation

Read counts are 0-based non-negative integers; `altReads > totalReads` (VAF > 1 from alignment artifacts) and negative counts throw `ArgumentOutOfRangeException`. `CalculateVAF` returns 0 for `totalReads == 0` (no coverage); `CalculateVAFConfidenceInterval` throws for `totalReads == 0` (an interval is undefined with no trials) and for confidence outside (0, 1) or ≠ 0.95. `AdjustVAFForPurity` throws for `vaf ∉ [0, 1]`, `purity ∉ (0, 1]` (correction divides by purity), or `ploidy ≤ 0`.

## 4. Algorithm

### 4.1 High-Level Steps

1. **CalculateVAF:** validate counts; return totalReads == 0 ? 0 : altReads/totalReads.
2. **CalculateVAFConfidenceInterval:** compute VAF; require coverage; resolve z from confidence; compute Wilson center and margin; clamp bounds to [0, 1].
3. **AdjustVAFForPurity:** validate inputs; compute average copies per cell 2(1−π) + π·n_tot; return vaf · avg / π.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- z(95%) = 1.96 [1] — the only supported confidence level; named constant `ZScore95`.
- Normal copy number = 2 [3][4] — named constant `NormalDiploidCopyNumber`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| All three methods | O(1) | O(1) | per-locus closed-form arithmetic |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CalculateVAF(int, int)`: empirical VAF.
- `OncologyAnalyzer.CalculateVAFConfidenceInterval(int, int, double)`: VAF + Wilson score interval.
- `OncologyAnalyzer.AdjustVAFForPurity(double, double, double)`: purity/ploidy correction.

### 5.2 Current Behavior

The Wilson bounds are clamped to [0, 1] to absorb floating-point drift at the exact boundaries p̂ = 0 / p̂ = 1; the unclamped values are already within [0, 1] mathematically. `CalculateVAF` reuses the existing private `CalculateVaf` validation shared with the somatic-calling path. **Search reuse:** N/A — this unit performs no substring/pattern search, so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Empirical VAF = altReads/totalReads [2].
- Wilson score interval center/margin with z = 1.96 for 95% [1].
- Purity/ploidy correction adjusted = vaf·(2(1−π) + π·n_tot)/π [3][4].

**Intentionally simplified:**

- Confidence level: only 0.95 (z = 1.96) is supported; **consequence:** other levels throw rather than returning an interval, because no other z value was retrieved from an authoritative source.
- Purity correction fixes normal ploidy at 2 and does not model multiplicity m or sub-clonality separately; **consequence:** the returned value is m·CCF for the diploid-normal case, not a decomposed CCF.

**Not implemented:**

- Mutect2-style Bayesian `AF` estimate; **users should rely on:** an external caller (GATK Mutect2) for the modelled allele fraction — this unit computes the empirical ratio by design [2].

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| totalReads = 0 (VAF) | returns 0 | no coverage ⇒ allele absent [2] |
| totalReads = 0 (CI) | ArgumentOutOfRangeException | interval undefined with no trials [1] |
| altReads > totalReads | ArgumentOutOfRangeException | invalid (VAF > 1 artifact) [2] |
| p̂ = 0 | Wilson lower = 0, upper > 0 | no overshoot, non-zero width [1] |
| p̂ = 1 | Wilson lower < 1, upper = 1 | no overshoot [1] |
| purity = 0 | ArgumentOutOfRangeException | correction divides by purity [3][4] |

### 6.2 Limitations

Sex chromosomes and non-diploid normal backgrounds are not modelled (normal ploidy fixed at 2). The correction returns m·CCF, not a separated multiplicity and cancer cell fraction. Only the 95% confidence level is supported.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
double vaf = OncologyAnalyzer.CalculateVAF(25, 100);                 // 0.25
var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(25, 100);   // ~[0.1755, 0.3430]
double ccf = OncologyAnalyzer.AdjustVAFForPurity(0.40, 0.80, 2);     // 1.0 (clonal het, diploid)
```

**Numerical walk-through:** for 25/100 at 95% (z = 1.96): p̂ = 0.25, denom = 1 + 1.96²/100 = 1.038416; center = (0.25 + 1.96²/200)/denom = 0.2592487; margin = (1.96/denom)·√(0.25·0.75/100 + 1.96²/40000) = 0.0837978 ⇒ [0.1754509, 0.3430465].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_CalculateVAF_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_CalculateVAF_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-VAF-001-Evidence.md](../../../docs/Evidence/ONCO-VAF-001-Evidence.md)
- Related algorithms: [Somatic_Mutation_Calling](../Oncology/Somatic_Mutation_Calling.md)

## 8. References

1. Wilson, E.B. 1927. Probable inference, the law of succession, and statistical inference. Journal of the American Statistical Association 22(158):209–212. https://doi.org/10.1080/01621459.1927.10502953 (formula retrieved via https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval)
2. GATK team / Broad Institute. FAQ for Mutect2; AlleleFraction. https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2
3. Tarabichi, M., Salcedo, A., Deshwar, A.G., et al. 2017. Principles of Reconstructing the Subclonal Architecture of Cancers. Cold Spring Harb Perspect Med. https://pmc.ncbi.nlm.nih.gov/articles/PMC5538405/
4. Househam, J., Caravagna, G., et al. 2024. Computational validation of clonal and subclonal copy number alterations from bulk tumor sequencing using CNAqc. Genome Biology 25:38. https://doi.org/10.1186/s13059-024-03170-5
