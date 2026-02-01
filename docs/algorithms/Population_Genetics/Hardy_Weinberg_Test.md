# Hardy-Weinberg Equilibrium Test

**Algorithm Group:** Population Genetics
**Test Unit:** POP-HW-001
**Last Updated:** 2026-02-01

---

## Overview

The Hardy-Weinberg equilibrium (HWE) test is a fundamental statistical test in population genetics used to determine whether genotype frequencies in a population conform to expected frequencies under the assumptions of random mating, no selection, no mutation, no migration, and infinite population size.

The test uses a chi-square goodness-of-fit approach to compare observed genotype counts against expected counts calculated from allele frequencies.

---

## Theory

### Hardy-Weinberg Principle

**Historical Background:**
The principle was independently derived by G. H. Hardy (1908) and Wilhelm Weinberg (1908). Hardy's paper was focused on refuting the incorrect notion that a dominant allele would automatically increase in frequency.

**Core Principle:**
In a population at Hardy-Weinberg equilibrium, allele and genotype frequencies remain constant across generations. For a biallelic locus with allele frequencies p and q (where p + q = 1), the expected genotype frequencies are:

$$f(AA) = p^2$$
$$f(Aa) = 2pq$$
$$f(aa) = q^2$$

### Allele Frequency Estimation

From observed genotype counts, allele frequencies are calculated as:

$$p = \frac{2n_{AA} + n_{Aa}}{2n}$$

$$q = 1 - p = \frac{2n_{aa} + n_{Aa}}{2n}$$

Where:
- $n_{AA}, n_{Aa}, n_{aa}$ = observed counts of each genotype
- $n = n_{AA} + n_{Aa} + n_{aa}$ = total sample size

### Chi-Square Test

**Test Statistic:**
$$\chi^2 = \sum_{i} \frac{(O_i - E_i)^2}{E_i}$$

Expanded for HWE:
$$\chi^2 = \frac{(n_{AA} - p^2 n)^2}{p^2 n} + \frac{(n_{Aa} - 2pqn)^2}{2pqn} + \frac{(n_{aa} - q^2 n)^2}{q^2 n}$$

**Degrees of Freedom:**
For biallelic HWE test: df = (# genotypes) - (# alleles) = 3 - 2 = 1

**Critical Values (df = 1):**

| Significance (α) | Critical χ² |
|------------------|-------------|
| 0.10 | 2.706 |
| 0.05 | 3.841 |
| 0.01 | 6.635 |

### Interpretation

| Chi-Square Result | P-value | Decision |
|-------------------|---------|----------|
| χ² < critical value | p ≥ α | Fail to reject HWE |
| χ² ≥ critical value | p < α | Reject HWE |

---

## Assumptions

The HWE model requires:

1. **Random mating** - no mate selection based on genotype
2. **No selection** - all genotypes have equal fitness
3. **No mutation** - allele frequencies unchanged by new mutations
4. **No migration** - closed population
5. **Infinite population size** - no genetic drift
6. **Equal allele frequencies in both sexes**

Violation of any assumption can cause deviation from HWE.

---

## Causes of HWE Deviation

**Biological causes:**
- Inbreeding (excess homozygotes)
- Population structure/stratification
- Selection (directional, balancing)
- Non-random mating

**Technical causes:**
- Genotyping errors
- Null alleles
- Sample contamination

---

## Implementation

### Method Signature

```csharp
public static HardyWeinbergResult TestHardyWeinberg(
    string variantId,
    int observedAA,
    int observedAa,
    int observedaa,
    double significanceLevel = 0.05)
```

### Return Type

```csharp
public readonly record struct HardyWeinbergResult(
    string VariantId,
    int ObservedAA,
    int ObservedAa,
    int Observedaa,
    double ExpectedAA,
    double ExpectedAa,
    double Expectedaa,
    double ChiSquare,
    double PValue,
    bool InEquilibrium);
```

### Algorithm Steps

1. Calculate total sample size: $n = n_{AA} + n_{Aa} + n_{aa}$
2. Handle edge case: if n = 0, return equilibrium with p-value = 1
3. Calculate allele frequencies: $p = (2n_{AA} + n_{Aa}) / 2n$
4. Calculate expected counts: $E_{AA} = p^2 n$, $E_{Aa} = 2pqn$, $E_{aa} = q^2 n$
5. Calculate chi-square statistic (skip terms where expected = 0)
6. Calculate p-value from chi-square CDF with df = 1
7. Determine equilibrium status: InEquilibrium = (p-value ≥ significanceLevel)

### Edge Cases

| Case | Handling |
|------|----------|
| n = 0 | Return equilibrium (χ² = 0, p = 1) |
| Fixed allele (all one genotype) | χ² = 0, in equilibrium |
| Expected count = 0 | Skip that term in χ² sum |

---

## Sources

1. Hardy, G. H. (1908). "Mendelian Proportions in a Mixed Population." Science, 28(706), 49-50.
2. Weinberg, W. (1908). "Über den Nachweis der Vererbung beim Menschen." Jahreshefte des Vereins für vaterländische Naturkunde in Württemberg, 64, 368-382.
3. Emigh, T. H. (1980). "A Comparison of Tests for Hardy-Weinberg Equilibrium." Biometrics, 36(4), 627-642.
4. Wikipedia: Hardy-Weinberg principle. https://en.wikipedia.org/wiki/Hardy-Weinberg_principle
5. Wigginton, J. E., Cutler, D. J., & Abecasis, G. R. (2005). "A Note on Exact Tests of Hardy-Weinberg Equilibrium." American Journal of Human Genetics, 76(5), 887-893.
