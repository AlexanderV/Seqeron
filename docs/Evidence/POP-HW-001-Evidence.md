# Evidence Document: POP-HW-001 (Hardy-Weinberg Equilibrium Test)

**Test Unit ID:** POP-HW-001  
**Algorithm:** Hardy-Weinberg Equilibrium Chi-Square Test  
**Date:** 2026-02-01  
**Status:** Evidence Gathered  

---

## 1. Sources Consulted

### Primary Sources

| # | Source | Type | URL | Accessed |
|---|--------|------|-----|----------|
| 1 | Wikipedia: Hardy-Weinberg principle | Encyclopedia | https://en.wikipedia.org/wiki/Hardy-Weinberg_principle | 2026-02-01 |
| 2 | Wikipedia: Chi-squared test | Encyclopedia | https://en.wikipedia.org/wiki/Chi-squared_test | 2026-02-01 |
| 3 | Hardy (1908) - Science | Original Paper | doi:10.1126/science.28.706.49 | 2026-02-01 |
| 4 | Weinberg (1908) - Jahreshefte | Original Paper | Historical reference | 2026-02-01 |
| 5 | Ford (1971) - Ecological Genetics | Textbook Example | ISBN: reference in Wikipedia | 2026-02-01 |
| 6 | Emigh (1980) - Biometrics | Chi-square test comparison | doi:10.2307/2556115 | 2026-02-01 |
| 7 | Wigginton et al. (2005) - AJHG | Exact test method | doi:10.1086/429864 | 2026-02-01 |

---

## 2. Algorithm Definitions

### 2.1 Hardy-Weinberg Equilibrium Principle

**Definition (Wikipedia - Hardy 1908, Weinberg 1908):**
> "In population genetics, the Hardy-Weinberg principle states that allele and genotype frequencies in a population will remain constant from generation to generation in the absence of other evolutionary influences."

**Expected Genotype Frequencies:**
For a biallelic locus with allele frequencies p (major) and q (minor) where p + q = 1:

$$f(AA) = p^2$$
$$f(Aa) = 2pq$$
$$f(aa) = q^2$$

Where:
- $p$ = frequency of allele A (major allele)
- $q$ = frequency of allele a (minor allele)
- $p + q = 1$

### 2.2 Allele Frequency Calculation from Genotypes

**Formula (Wikipedia):**
$$p = \frac{2 \times n_{AA} + n_{Aa}}{2n}$$
$$q = 1 - p = \frac{2 \times n_{aa} + n_{Aa}}{2n}$$

Where:
- $n_{AA}$ = observed count of homozygous major genotype
- $n_{Aa}$ = observed count of heterozygous genotype  
- $n_{aa}$ = observed count of homozygous minor genotype
- $n = n_{AA} + n_{Aa} + n_{aa}$ = total sample size

### 2.3 Expected Genotype Counts

**Formula:**
$$E(AA) = p^2 \times n$$
$$E(Aa) = 2pq \times n$$
$$E(aa) = q^2 \times n$$

### 2.4 Chi-Square Test for HWE

**Definition (Wikipedia - Pearson Chi-Square Test):**
> "Pearson's chi-squared test is used to determine whether there is a statistically significant difference between the expected frequencies and the observed frequencies."

**Formula (Wikipedia - HWE example):**
$$\chi^2 = \sum \frac{(O - E)^2}{E} = \frac{(n_{AA} - E_{AA})^2}{E_{AA}} + \frac{(n_{Aa} - E_{Aa})^2}{E_{Aa}} + \frac{(n_{aa} - E_{aa})^2}{E_{aa}}$$

**Degrees of Freedom (Wikipedia):**
> "There is 1 degree of freedom (degrees of freedom for test for Hardy-Weinberg proportions are # genotypes − # alleles)."

For biallelic case: df = 3 genotypes - 2 alleles = 1

**Significance (Wikipedia):**
> "The 5% significance level for 1 degree of freedom is 3.84"

---

## 3. Published Test Datasets

### 3.1 Ford's Scarlet Tiger Moth Data (Wikipedia Example)

**Source:** Ford (1971) Ecological Genetics, cited in Wikipedia HWE article

| Phenotype | White-spotted (AA) | Intermediate (Aa) | Little spotting (aa) | Total |
|-----------|-------------------|-------------------|----------------------|-------|
| Count | 1469 | 138 | 5 | 1612 |

**Calculated Values:**
- $p = \frac{2 \times 1469 + 138}{2 \times 1612} = \frac{3076}{3224} = 0.954$
- $q = 1 - 0.954 = 0.046$
- $E(AA) = 0.954^2 \times 1612 = 1467.4$
- $E(Aa) = 2 \times 0.954 \times 0.046 \times 1612 = 141.2$
- $E(aa) = 0.046^2 \times 1612 = 3.4$
- $\chi^2 = \frac{(1469-1467.4)^2}{1467.4} + \frac{(138-141.2)^2}{141.2} + \frac{(5-3.4)^2}{3.4}$
- $\chi^2 = 0.001 + 0.073 + 0.756 = 0.83$
- **Result:** p-value > 0.05, population IS in HWE

### 3.2 Perfect HWE (Synthetic Dataset)

**When allele frequencies exactly match genotype expectations:**
For p = 0.5, n = 100:
- $E(AA) = 0.25 \times 100 = 25$
- $E(Aa) = 0.50 \times 100 = 50$
- $E(aa) = 0.25 \times 100 = 25$
- If observed = expected, $\chi^2 = 0$

### 3.3 Excess Heterozygotes (Deviation from HWE)

**Example indicating non-random mating or selection:**
- Observed: AA=10, Aa=80, aa=10 (n=100)
- $p = \frac{20 + 80}{200} = 0.5$
- Expected: AA=25, Aa=50, aa=25
- $\chi^2 = \frac{(10-25)^2}{25} + \frac{(80-50)^2}{50} + \frac{(10-25)^2}{25}$
- $\chi^2 = 9 + 18 + 9 = 36$ (very high, significant deviation)

---

## 4. Edge Cases (from Wikipedia and Implementation Analysis)

### 4.1 Zero Sample Size
- **Input:** n = 0 (no observations)
- **Expected:** Return InEquilibrium = true, ChiSquare = 0, PValue = 1
- **Rationale:** No data to test, default to equilibrium (ASSUMPTION)

### 4.2 Fixed Allele (Monomorphic)
- **Input:** All AA (e.g., AA=100, Aa=0, aa=0)
- **Expected:** p = 1.0, q = 0.0
- **Expected counts:** E(AA)=100, E(Aa)=0, E(aa)=0
- **ChiSquare:** 0 (observed = expected)
- **Result:** InEquilibrium = true

### 4.3 All Heterozygotes
- **Input:** AA=0, Aa=100, aa=0
- **Expected:** p = 0.5, q = 0.5
- **Expected counts:** E(AA)=25, E(Aa)=50, E(aa)=25
- **ChiSquare:** Very high (deviation from HWE)
- **Interpretation:** Impossible under random mating, indicates heterozygote advantage or assortative mating

### 4.4 Single Individual
- **Input:** n = 1 (e.g., AA=1, Aa=0, aa=0)
- **Expected:** Test still runs, but statistically meaningless
- **Chi-square calculation:** Division might have issues if expected = 0

### 4.5 Division by Zero Protection
- **When:** Expected count = 0 for a genotype class
- **Implementation:** Skip that term in chi-square sum (Wikipedia formula uses conditional)

---

## 5. Deviations from HWE (Wikipedia)

**Causes of deviation:**
1. **Non-random mating** (inbreeding → excess homozygotes)
2. **Selection** (directional, balancing)
3. **Mutation**
4. **Migration/gene flow**
5. **Small population size** (genetic drift)
6. **Genotyping errors** (real-world issue in bioinformatics)

**Wikipedia Quote:**
> "In real world genotype data, deviations from Hardy-Weinberg Equilibrium may be a sign of genotyping error."

---

## 6. Significance Testing (Wikipedia)

### 6.1 Chi-Square Critical Values

| df | α = 0.10 | α = 0.05 | α = 0.01 |
|----|----------|----------|----------|
| 1 | 2.706 | 3.841 | 6.635 |

### 6.2 Interpretation

| P-value | Decision | Interpretation |
|---------|----------|----------------|
| p ≥ α | Fail to reject H₀ | Population is consistent with HWE |
| p < α | Reject H₀ | Significant deviation from HWE |

Default significance level: α = 0.05

---

## 7. Implementation Notes

### Current Implementation Analysis

The `TestHardyWeinberg` method in `PopulationGeneticsAnalyzer.cs`:

1. Calculates allele frequencies from observed genotype counts
2. Computes expected counts under HWE
3. Performs chi-square goodness-of-fit test
4. Returns p-value calculated via chi-square CDF (df=1)
5. Uses default significance level of 0.05

**Chi-Square CDF Implementation:**
- Uses lower incomplete gamma function approximation
- Standard approach for calculating p-values from chi-square distribution

---

## 8. Test Case Summary (Evidence-Based)

| ID | Test Case | Source | Expected Result |
|----|-----------|--------|-----------------|
| HW-E01 | Ford's moth data (1469, 138, 5) | Wikipedia | χ² ≈ 0.83, InEquilibrium = true |
| HW-E02 | Perfect HWE (25, 50, 25) | Mathematical definition | χ² = 0, InEquilibrium = true |
| HW-E03 | Excess heterozygotes (10, 80, 10) | Wikipedia deviation | χ² >> 3.84, InEquilibrium = false |
| HW-E04 | Zero samples (0, 0, 0) | Edge case | InEquilibrium = true, PValue = 1 |
| HW-E05 | Fixed allele (100, 0, 0) | Definition | InEquilibrium = true |
| HW-E06 | All heterozygotes (0, 100, 0) | Extreme deviation | InEquilibrium = false |
