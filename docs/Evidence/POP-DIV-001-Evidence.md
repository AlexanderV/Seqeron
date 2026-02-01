# Evidence Document: POP-DIV-001 (Diversity Statistics)

**Test Unit ID:** POP-DIV-001
**Algorithm:** Nucleotide Diversity, Watterson's Theta, Tajima's D
**Date:** 2026-02-01
**Status:** Evidence Gathered

---

## 1. Sources Consulted

### Primary Sources

| # | Source | Type | URL | Accessed |
|---|--------|------|-----|----------|
| 1 | Wikipedia: Nucleotide diversity | Encyclopedia | https://en.wikipedia.org/wiki/Nucleotide_diversity | 2026-02-01 |
| 2 | Wikipedia: Watterson estimator | Encyclopedia | https://en.wikipedia.org/wiki/Watterson_estimator | 2026-02-01 |
| 3 | Wikipedia: Tajima's D | Encyclopedia | https://en.wikipedia.org/wiki/Tajima%27s_D | 2026-02-01 |
| 4 | Wikipedia: Zygosity (Heterozygosity) | Encyclopedia | https://en.wikipedia.org/wiki/Zygosity | 2026-02-01 |
| 5 | Nei & Li (1979) - PNAS | Original Paper | doi:10.1073/pnas.76.10.5269 | 2026-02-01 |
| 6 | Watterson (1975) - TPB | Original Paper | doi:10.1016/0040-5809(75)90020-9 | 2026-02-01 |
| 7 | Tajima (1989) - Genetics | Original Paper | doi:10.1093/genetics/123.3.585 | 2026-02-01 |
| 8 | Hartl & Clark (2007) - Principles of Population Genetics | Textbook | ISBN 978-0-87893-308-2 | 2026-02-01 |

---

## 2. Algorithm Definitions

### 2.1 Nucleotide Diversity (π)

**Definition (Wikipedia - Nei & Li 1979):**
> "Nucleotide diversity is a concept in molecular genetics which is used to measure the degree of polymorphism within a population. One commonly used measure of nucleotide diversity was first introduced by Nei and Li in 1979. This measure is defined as the average number of nucleotide differences per site between two DNA sequences in all possible pairs in the sample population, and is denoted by π."

**Formula (Wikipedia):**
$$\hat{\pi} = \frac{n}{n-1} \sum_{i} \sum_{j} x_i x_j \pi_{ij}$$

Where:
- $x_i$ and $x_j$ are the respective frequencies of the $i$-th and $j$-th sequences
- $\pi_{ij}$ is the number of nucleotide differences per nucleotide site between sequences $i$ and $j$
- $n$ is the number of sequences in the sample

**Simplified for equal sample sizes:**
$$\pi = \frac{\sum_{i<j} d_{ij}}{\binom{n}{2} \cdot L}$$

Where:
- $d_{ij}$ = number of differences between sequences $i$ and $j$
- $L$ = sequence length
- $\binom{n}{2}$ = number of pairwise comparisons = $n(n-1)/2$

**Key Properties:**
- Measures polymorphism at nucleotide level
- Related to expected heterozygosity
- **Invariant:** $\pi \geq 0$

### 2.2 Watterson's Theta (θ_W)

**Definition (Wikipedia - Watterson 1975):**
> "The Watterson estimator is a method for describing the genetic diversity in a population. It is estimated by counting the number of polymorphic sites."

**Formula (Wikipedia):**
$$\hat{\theta}_W = \frac{K}{a_n}$$

Where:
- $K$ = number of segregating sites (polymorphic positions)
- $a_n = \sum_{i=1}^{n-1} \frac{1}{i}$ (the $(n-1)$-th harmonic number)

**Per-site Watterson's theta:**
$$\theta_W = \frac{S}{a_n \cdot L}$$

Where:
- $S$ = segregating sites
- $L$ = sequence length

**Key Properties:**
- Based on coalescent theory
- Unbiased estimator when assumptions are met
- Can be biased by population structure or exponential growth
- **Invariant:** $\theta_W \geq 0$

### 2.3 Tajima's D

**Definition (Wikipedia - Tajima 1989):**
> "Tajima's D is a population genetic test statistic created by and named after the Japanese researcher Fumio Tajima. Tajima's D is computed as the difference between two measures of genetic diversity: the mean number of pairwise differences and the number of segregating sites, each scaled so that they are expected to be the same in a neutrally evolving population of constant size."

**Purpose:**
> "The purpose of Tajima's D test is to distinguish between a DNA sequence evolving randomly ('neutrally') and one evolving under a non-random process, including directional selection or balancing selection, demographic expansion or contraction, genetic hitchhiking, or introgression."

**Formula (Wikipedia):**
$$D = \frac{d}{\sqrt{\hat{V}(d)}} = \frac{\hat{k} - \frac{S}{a_1}}{\sqrt{e_1 S + e_2 S(S-1)}}$$

Where:
- $\hat{k}$ = average number of pairwise differences (π × L)
- $S$ = number of segregating sites
- $a_1 = \sum_{i=1}^{n-1} \frac{1}{i}$
- $a_2 = \sum_{i=1}^{n-1} \frac{1}{i^2}$
- $e_1, e_2$ = constants derived from sample size

**Variance components:**
- $b_1 = \frac{n+1}{3(n-1)}$
- $b_2 = \frac{2(n^2 + n + 3)}{9n(n-1)}$
- $c_1 = b_1 - \frac{1}{a_1}$
- $c_2 = b_2 - \frac{n+2}{a_1 \cdot n} + \frac{a_2}{a_1^2}$
- $e_1 = \frac{c_1}{a_1}$
- $e_2 = \frac{c_2}{a_1^2 + a_2}$

**Interpretation (Wikipedia):**

| Value | Mathematical reason | Biological interpretation |
|-------|---------------------|---------------------------|
| D = 0 | π ≈ θ (Observed ≈ Expected) | Neutral evolution, mutation-drift equilibrium |
| D < 0 | π < θ (excess rare alleles) | Recent selective sweep, population expansion after bottleneck |
| D > 0 | π > θ (lack of rare alleles) | Balancing selection, sudden population contraction |

**Rule of thumb (Wikipedia):**
> "A very rough rule of thumb to significance is that values greater than +2 or less than -2 are likely to be significant."

### 2.4 Segregating Sites (S)

**Definition:**
A segregating site is a position in the sequence alignment where not all sequences have the same nucleotide. Also called a polymorphic site or SNP.

**Calculation:**
- For each position in the alignment, check if any sequence differs from the first (or any other reference)
- Count the number of such positions

### 2.5 Heterozygosity

**Observed Heterozygosity (H_o) - Wikipedia:**
$$H_o = \frac{\sum_{i=1}^{n} (1 \text{ if } a_{i1} \neq a_{i2})}{n}$$

For sequence data (adapted): fraction of polymorphic sites

**Expected Heterozygosity (H_e) - Wikipedia:**
$$H_e = 1 - \sum_{i=1}^{m} f_i^2$$

Where:
- $m$ = number of alleles at the locus
- $f_i$ = frequency of the $i$-th allele

For sequence data: average expected heterozygosity across all positions

---

## 3. Test Cases from Sources

### 3.1 Wikipedia Tajima's D Example

**Source:** Wikipedia Tajima's D - Example section

**Dataset:**
```
Position:  1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
Person Y:  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Person A:  0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0
Person B:  0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0
Person C:  0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0
Person D:  0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0
```

**Parameters (n=5):**
- Segregating sites: S = 4 (positions 3, 7, 13, 19)
- $a_1 = 1 + 1/2 + 1/3 + 1/4 = 2.083...$

**Pairwise differences:**
- Y vs A: 3, Y vs B: 2, Y vs C: 2, Y vs D: 3
- A vs B: 1, A vs C: 3, A vs D: 2
- B vs C: 2, B vs D: 1
- C vs D: 1

- Total: 20 differences
- Comparisons: 10
- Average: $\hat{k}$ = 2.0

**Watterson's estimate:**
- M = S/a₁ = 4/2.083 = 1.92

**Tajima's d (unnormalized):**
- d = 2.0 - 1.92 = 0.08

### 3.2 Nucleotide Diversity Test Cases

| # | Test Case | Source | Input | Expected |
|---|-----------|--------|-------|----------|
| ND-1 | Identical sequences | Definition | All same sequences | π = 0 |
| ND-2 | All different (two seqs) | Definition | "AAAA", "TTTT" | π = 1.0 |
| ND-3 | Single sequence | Edge case | n = 1 | π = 0 (undefined) |
| ND-4 | Wikipedia example | Wikipedia | n=5, S=4, L=20 | k̂ = 2.0, π = 0.1 |

### 3.3 Watterson's Theta Test Cases

| # | Test Case | Source | Input | Expected |
|---|-----------|--------|-------|----------|
| WT-1 | Wikipedia harmonic example | Wikipedia | S=10, n=10, L=1000 | θ ≈ 0.00353 |
| WT-2 | No segregating sites | Definition | S=0 | θ = 0 |
| WT-3 | Sample size n=2 | Edge case | S=5, n=2, L=100 | θ = 5/(1×100) = 0.05 |
| WT-4 | Wikipedia Tajima example | Wikipedia | S=4, n=5, L=20 | θ = 4/(2.083×20) ≈ 0.096 |

### 3.4 Tajima's D Test Cases

| # | Test Case | Source | Input | Expected |
|---|-----------|--------|-------|----------|
| TD-1 | Neutral evolution (π ≈ θ) | Wikipedia | π = θ | D ≈ 0 |
| TD-2 | Positive selection (π < θ) | Wikipedia | π << θ | D < 0 |
| TD-3 | Balancing selection (π > θ) | Wikipedia | π >> θ | D > 0 |
| TD-4 | No segregating sites | Edge case | S = 0 | D = 0 |
| TD-5 | Minimum sample size | Edge case | n < 3 | D = 0 (undefined) |

---

## 4. Edge Cases and Boundary Conditions

### From Wikipedia/Standard Definitions

| Edge Case | Expected Behavior | Source |
|-----------|-------------------|--------|
| n = 0 (empty input) | Return zeros / handle gracefully | ASSUMPTION |
| n = 1 (single sequence) | π = 0, θ undefined (return 0) | Definition |
| n = 2 (minimum for comparison) | π calculable, Tajima's D undefined | Tajima (1989) |
| n < 3 for Tajima's D | D = 0 or undefined | Wikipedia |
| S = 0 (monomorphic) | π = 0, θ = 0, D = 0 | Definition |
| All sequences identical | S = 0, π = 0, θ = 0 | Definition |
| Sequences of different lengths | Error or handle first length | ASSUMPTION |

---

## 5. Implementation Notes

### Current Implementation (PopulationGeneticsAnalyzer.cs)

**CalculateNucleotideDiversity:**
- Returns 0 for n < 2 ✓
- Uses pairwise comparison ✓
- Divides by comparisons × length ✓

**CalculateWattersonTheta:**
- Returns 0 for n < 2 or L ≤ 0 ✓
- Uses harmonic number a₁ ✓
- Divides by a₁ × L ✓

**CalculateTajimasD:**
- Returns 0 for S = 0 or n < 3 ✓
- Implements full Tajima (1989) formula ✓
- Uses harmonic numbers a₁, a₂ ✓

**CalculateDiversityStatistics:**
- Aggregates all metrics ✓
- Includes heterozygosity ✓

### Deviations from Standard

1. **Heterozygosity calculation:** Implementation uses polymorphic sites for observed heterozygosity rather than diploid genotype data. This is an adaptation for haploid sequence data. (ASSUMPTION: This is acceptable for the intended use case)

---

## 6. References

1. Nei, M.; Li, W.-H. (1979). "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". PNAS. 76(10): 5269-73.

2. Watterson, G.A. (1975). "On the number of segregating sites in genetical models without recombination". Theoretical Population Biology. 7(2): 256-276.

3. Tajima, F. (1989). "Statistical method for testing the neutral mutation hypothesis by DNA polymorphism". Genetics. 123(3): 585-95.

4. Hartl, D.L.; Clark, A.G. (2007). Principles of Population Genetics (4th ed.). Sinauer Associates. ISBN 978-0-87893-308-2.
