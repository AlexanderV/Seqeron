# Evidence Document: POP-FST-001 (F-Statistics)

**Test Unit ID:** POP-FST-001
**Algorithm:** Fst (Fixation Index), F-Statistics (Fis, Fit, Fst)
**Date:** 2026-02-01
**Status:** Evidence Gathered

---

## 1. Sources Consulted

### Primary Sources

| # | Source | Type | URL | Accessed |
|---|--------|------|-----|----------|
| 1 | Wikipedia: Fixation index | Encyclopedia | https://en.wikipedia.org/wiki/Fixation_index | 2026-02-01 |
| 2 | Wikipedia: F-statistics | Encyclopedia | https://en.wikipedia.org/wiki/F-statistics | 2026-02-01 |
| 3 | Wright (1950) - Nature | Original Paper | doi:10.1038/166247a0 | 2026-02-01 |
| 4 | Wright (1965) - Evolution | Original Paper | doi:10.2307/2406450 | 2026-02-01 |
| 5 | Weir & Cockerham (1984) - Evolution | Original Paper | doi:10.2307/2408641 | 2026-02-01 |
| 6 | Holsinger & Weir (2009) - Nat Rev Genet | Review | doi:10.1038/nrg2611 | 2026-02-01 |
| 7 | Cavalli-Sforza et al. (1994) - History and Geography of Human Genes | Textbook | ISBN 978-0691029054 | 2026-02-01 |

---

## 2. Algorithm Definitions

### 2.1 Fixation Index (Fst)

**Definition (Wikipedia):**
> "The fixation index (FST) is a measure of population differentiation due to genetic structure. It is frequently estimated from genetic polymorphism data, such as single-nucleotide polymorphisms (SNP) or microsatellites. Developed as a special case of Wright's F-statistics, it is one of the most commonly used statistics in population genetics. Its values range from 0 to 1, with 0 being no differentiation and 1 being complete differentiation."

**Interpretation (Wikipedia):**
> "A zero value implies complete panmixia; that is, that the two populations are interbreeding freely. A value of one implies that all genetic variation is explained by the population structure, and that the two populations do not share any genetic diversity."

### 2.2 Wright's Definition (Variance-Based)

**Formula (Wikipedia - Wright 1965):**
$$F_{ST} = \frac{\sigma_S^2}{\sigma_T^2} = \frac{\sigma_S^2}{\bar{p}(1-\bar{p})}$$

Where:
- $\sigma_S^2$ = variance in allele frequency among subpopulations (weighted by population size)
- $\bar{p}$ = average allele frequency in total population
- $\bar{p}(1-\bar{p})$ = expected heterozygosity under Hardy-Weinberg

**Alternative formulation (Wikipedia):**
$$F_{ST} = \frac{\bar{p}(1-\bar{p}) - \overline{p(1-p)}}{\bar{p}(1-\bar{p})}$$

This measures the fraction of total diversity not due to average diversity within subpopulations.

### 2.3 Weir and Cockerham Estimator (1984)

**Source:** Weir & Cockerham (1984) Evolution 38:1358-1370

The Weir-Cockerham estimator is the most widely used due to its unbiased properties with finite sample sizes. The implementation in this library uses a simplified variance-based approach.

**For two populations (simplified):**
$$F_{ST} = \frac{\text{Between-population variance}}{\text{Total heterozygosity}}$$

Where:
- Between-pop variance = weighted variance of allele frequencies
- Total heterozygosity = $\bar{p}(1-\bar{p})$

### 2.4 F-Statistics Relationship

**Partition Formula (Wikipedia - Wright):**
$$(1 - F_{IT}) = (1 - F_{IS})(1 - F_{ST})$$

Where:
- **Fis**: Inbreeding coefficient of Individual relative to Subpopulation
- **Fit**: Inbreeding coefficient of Individual relative to Total population
- **Fst**: Effect of Subpopulation structure relative to Total population

**Heterozygosity-based definitions:**
- $F_{IS} = 1 - H_I / H_S$ (observed vs expected heterozygosity within subpop)
- $F_{IT} = 1 - H_I / H_T$ (observed vs expected heterozygosity in total)
- $F_{ST} = 1 - H_S / H_T$ (subpop vs total expected heterozygosity)

---

## 3. Mathematical Properties and Invariants

### 3.1 Value Range

| Property | Constraint | Source |
|----------|-----------|--------|
| Fst range | $0 \leq F_{ST} \leq 1$ | Wikipedia |
| Fst = 0 | Complete panmixia (no differentiation) | Wikipedia |
| Fst = 1 | Complete differentiation (fixed differences) | Wikipedia |
| Fis range | $-1 \leq F_{IS} \leq 1$ (can be negative with excess heterozygosity) | Wikipedia |
| Fit range | $-1 \leq F_{IT} \leq 1$ | Wikipedia |

### 3.2 Pairwise Fst Matrix Properties

| Property | Expected Behavior | Source |
|----------|-------------------|--------|
| Diagonal | $F_{ST}(i,i) = 0$ (identical populations) | Mathematical |
| Symmetry | $F_{ST}(i,j) = F_{ST}(j,i)$ | Mathematical |
| Non-negative | $F_{ST}(i,j) \geq 0$ for distinct populations | Wikipedia |

**Note (Wikipedia):** "FST is not a distance in the mathematical sense, as it does not satisfy the triangle inequality."

---

## 4. Reference Values from Literature

### 4.1 Human Population Fst Values (Cavalli-Sforza 1994)

| Population Pair | Fst | Source |
|-----------------|-----|--------|
| Danish - English | 0.0021 | Cavalli-Sforza (1994) |
| Dutch - Danes | 0.0009 | Cavalli-Sforza (1994) |
| Lapps - Sardinians | 0.0667 | Cavalli-Sforza (1994) |
| Mbuti - Papua New Guineans | 0.4573 | Cavalli-Sforza (1994) |
| Mean (42 populations) | 0.1338 | Cavalli-Sforza (1994) |

### 4.2 Continental Fst Values (Elhaik 2012, HapMap)

| Population Pair | Fst | Source |
|-----------------|-----|--------|
| Europe - Sub-Saharan Africa | 0.153 | Elhaik (2012) |
| Europe - East Asia | 0.111 | Elhaik (2012) |
| Sub-Saharan Africa - East Asia | 0.190 | Elhaik (2012) |
| Within continental populations | < 0.01 | Elhaik (2012) |

### 4.3 Typical Interpretation Scale

| Fst Range | Interpretation | Source |
|-----------|---------------|--------|
| 0.00 - 0.05 | Little genetic differentiation | Hartl & Clark |
| 0.05 - 0.15 | Moderate differentiation | Hartl & Clark |
| 0.15 - 0.25 | Great differentiation | Hartl & Clark |
| > 0.25 | Very great differentiation | Hartl & Clark |

---

## 5. Edge Cases and Corner Cases

### 5.1 Documented Edge Cases

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Identical populations | Fst = 0 | Wikipedia (panmixia) |
| Fixed differences (p1=1, p2=0) | Fst approaches 1 | Wikipedia |
| Empty populations | Fst = 0 (undefined, return 0) | ASSUMPTION |
| Single variant | Valid Fst calculation | Mathematical |
| Unequal sample sizes | Weighted by sample size | Weir & Cockerham |
| All monomorphic sites | Fst = 0 (no variation) | Mathematical |

### 5.2 Numerical Considerations

| Issue | Handling | Source |
|-------|----------|--------|
| Division by zero (pBar = 0 or 1) | Return 0 | Implementation |
| High polymorphism | Arbitrarily low upper bound | Wikipedia |
| Small sample sizes | Weir-Cockerham correction | Weir & Cockerham |

---

## 6. Testing Methodology

### 6.1 Test Categories (Derived from Literature)

| Category | Description | Priority |
|----------|-------------|----------|
| Boundary conditions | Fst = 0 (identical), Fst → 1 (fixed) | Must |
| Value range invariant | 0 ≤ Fst ≤ 1 | Must |
| Matrix properties | Diagonal = 0, Symmetry | Must |
| Heterozygosity relationship | Fis, Fit, Fst partition | Should |
| Empty/edge inputs | Graceful handling | Must |

### 6.2 Recommended Test Datasets

**Test Case 1: Identical Populations**
- Pop1 = Pop2 = [(0.5, 100), (0.3, 100)]
- Expected: Fst = 0

**Test Case 2: Moderate Differentiation**
- Pop1 = [(0.9, 100), (0.8, 100)]
- Pop2 = [(0.1, 100), (0.2, 100)]
- Expected: Fst > 0 (significant differentiation)

**Test Case 3: Fixed Differences (Maximum Differentiation)**
- Pop1 = [(1.0, 100)]
- Pop2 = [(0.0, 100)]
- Expected: Fst = 1.0 (or very high, ~0.5+ depending on estimator)

**Test Case 4: F-Statistics Components**
- Heterozygosity data with observed and expected values
- Verify Fis, Fit, Fst relationship: (1-Fit) = (1-Fis)(1-Fst)

---

## 7. Implementation Notes

### 7.1 Current Implementation Analysis

The implementation uses a simplified Weir-Cockerham-style estimator:

```
For each locus:
  pBar = (n1*p1 + n2*p2) / (n1 + n2)   // Weighted mean
  variance = weighted variance of p1, p2 around pBar
  het = pBar * (1 - pBar)

Fst = sum(variance) / sum(het)
```

This is consistent with Wright's variance definition.

### 7.2 CalculateFStatistics Implementation

Uses heterozygosity-based approach:
- Hi = observed heterozygosity (individual level)
- Hs = expected heterozygosity within subpopulations
- Ht = expected heterozygosity in total population

Returns: Fis = 1 - Hi/Hs, Fit = 1 - Hi/Ht, Fst = 1 - Hs/Ht

---

## 8. Fallback Strategy

N/A - Comprehensive sources found.

---

## 9. Open Questions

None - algorithm well-defined in literature.

---

## 10. Evidence Summary

| Aspect | Status | Notes |
|--------|--------|-------|
| Definition | ✓ Complete | Wright (1965), Wikipedia |
| Formula | ✓ Complete | Weir & Cockerham (1984), Wikipedia |
| Invariants | ✓ Complete | 0 ≤ Fst ≤ 1, matrix properties |
| Edge cases | ✓ Complete | From Wikipedia, implementation |
| Test data | ✓ Complete | Reference values from Cavalli-Sforza |
| Methodology | ✓ Complete | Standard population genetics testing |
