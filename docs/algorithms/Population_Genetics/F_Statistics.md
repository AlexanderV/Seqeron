# F-Statistics (Fst, Fis, Fit)

**Algorithm Group:** Population Genetics
**Test Unit:** POP-FST-001
**Last Updated:** 2026-02-01

---

## Overview

F-statistics (fixation indices) are a set of measures that describe the partitioning of genetic variation within and among populations. Developed by Sewall Wright in the 1920s-1960s, they are fundamental tools in population genetics for quantifying population structure.

---

## Mathematical Definitions

### Fst (Fixation Index)

Fst measures population differentiation due to genetic structure.

**Wright's variance-based definition:**

$$F_{ST} = \frac{\sigma_S^2}{\bar{p}(1-\bar{p})}$$

Where:
- $\sigma_S^2$ = variance in allele frequency among subpopulations
- $\bar{p}$ = mean allele frequency across populations

**Heterozygosity-based definition:**

$$F_{ST} = \frac{H_T - H_S}{H_T}$$

Where:
- $H_T$ = expected heterozygosity in total population
- $H_S$ = mean expected heterozygosity within subpopulations

### F-Statistics Partition

The three F-statistics are related by:

$$(1 - F_{IT}) = (1 - F_{IS})(1 - F_{ST})$$

| Statistic | Measures | Range |
|-----------|----------|-------|
| Fis | Inbreeding within subpopulation | -1 to 1 |
| Fit | Inbreeding in total population | -1 to 1 |
| Fst | Population differentiation | 0 to 1 |

---

## Interpretation

| Fst Value | Interpretation |
|-----------|----------------|
| 0 | No differentiation (panmixia) |
| 0.00 - 0.05 | Little genetic differentiation |
| 0.05 - 0.15 | Moderate differentiation |
| 0.15 - 0.25 | Great differentiation |
| > 0.25 | Very great differentiation |
| 1 | Complete differentiation (fixed differences) |

---

## Implementation

### CalculateFst

Calculates pairwise Fst between two populations using a Weir-Cockerham style variance estimator.

**Signature:**
```csharp
public static double CalculateFst(
    IEnumerable<(double AlleleFreq, int SampleSize)> population1,
    IEnumerable<(double AlleleFreq, int SampleSize)> population2)
```

**Algorithm:**
1. For each variant locus:
   - Calculate weighted mean frequency: $\bar{p} = (n_1 p_1 + n_2 p_2) / (n_1 + n_2)$
   - Calculate variance: $\sigma^2 = [(p_1 - \bar{p})^2 n_1 + (p_2 - \bar{p})^2 n_2] / (n_1 + n_2)$
   - Calculate expected heterozygosity: $H = \bar{p}(1 - \bar{p})$
2. Sum across loci: $F_{ST} = \sum \sigma^2 / \sum H$

**Complexity:** O(n) where n = number of loci

### CalculateFStatistics

Calculates all three F-statistics (Fis, Fit, Fst) from heterozygosity data.

**Signature:**
```csharp
public static FStatistics CalculateFStatistics(
    string pop1Name,
    string pop2Name,
    IEnumerable<(int HetObs1, int N1, int HetObs2, int N2,
                 double AlleleFreq1, double AlleleFreq2)> variantData)
```

**Algorithm:**
1. For each variant:
   - Accumulate observed heterozygosity (Hi)
   - Calculate expected heterozygosity within subpops (Hs)
   - Calculate total expected heterozygosity (Ht)
2. Compute:
   - $F_{IS} = 1 - H_I / H_S$
   - $F_{IT} = 1 - H_I / H_T$
   - $F_{ST} = 1 - H_S / H_T$

### CalculatePairwiseFst

Calculates an Fst matrix for multiple populations.

**Signature:**
```csharp
public static double[,] CalculatePairwiseFst(
    IEnumerable<(string PopulationId,
                 IReadOnlyList<(double AlleleFreq, int SampleSize)> Variants)> populations)
```

**Properties:**
- Diagonal entries are 0
- Matrix is symmetric: Fst(i,j) = Fst(j,i)
- All values ≥ 0

**Complexity:** O(k² × n) where k = populations, n = loci

---

## Invariants

| Invariant | Description |
|-----------|-------------|
| Value range | 0 ≤ Fst ≤ 1 |
| Identical populations | Fst = 0 |
| Self-comparison | Fst(pop, pop) = 0 |
| Symmetry | Fst(A, B) = Fst(B, A) |
| Partition | (1-Fit) = (1-Fis)(1-Fst) |

---

## Edge Cases

| Case | Behavior |
|------|----------|
| Empty populations | Returns 0 |
| Single locus | Valid calculation |
| Unequal sample sizes | Weighted by n |
| Monomorphic loci | No contribution (het = 0) |
| Fixed differences | Fst approaches 1 |

---

## References

1. Wright, S. (1950). Genetical structure of populations. Nature 166:247-249.
2. Wright, S. (1965). The interpretation of population structure by F-statistics with special regard to systems of mating. Evolution 19:395-420.
3. Weir, B.S. & Cockerham, C.C. (1984). Estimating F-statistics for the analysis of population structure. Evolution 38:1358-1370.
4. Holsinger, K.E. & Weir, B.S. (2009). Genetics in geographically structured populations: defining, estimating and interpreting Fst. Nat Rev Genet 10:639-650.

---

## See Also

- [Allele Frequency](Allele_Frequency.md)
- [Hardy-Weinberg Test](Hardy_Weinberg_Test.md)
- [Diversity Statistics](Diversity_Statistics.md)
