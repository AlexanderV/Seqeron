# Alpha Diversity

## Overview

| Property | Value |
|----------|-------|
| **Algorithm** | Alpha Diversity Calculation |
| **Category** | Metagenomics |
| **Time Complexity** | O(n) where n = number of taxa |
| **Space Complexity** | O(n) |

## Description

Alpha diversity (α-diversity) measures the species diversity within a single sample or site at a local scale. The term was introduced by R. H. Whittaker (1960, 1972) as part of the alpha-beta-gamma diversity framework. Alpha diversity quantifies both the number of species present (richness) and their relative abundances (evenness).

## Indices Computed

### 1. Observed Species (Richness)

The simplest diversity measure: count of distinct species with non-zero abundance.

$$S_{obs} = \text{count}(\{i : p_i > 0\})$$

**Source:** Wikipedia — Species Richness

### 2. Shannon Index (H)

Quantifies uncertainty in predicting species identity of a randomly chosen individual.

$$H = -\sum_{i=1}^{S} p_i \ln(p_i)$$

Where:
- $S$ = number of species with non-zero abundance
- $p_i$ = proportional abundance of species $i$

**Properties:**
- $H = 0$ when only one species (certainty)
- $H = \ln(S)$ when all species equally abundant (maximum entropy)
- Higher values indicate higher diversity

**Source:** Shannon (1948), Wikipedia — Diversity Index

### 3. Simpson Index (λ)

Probability that two randomly chosen individuals belong to the same species.

$$\lambda = \sum_{i=1}^{S} p_i^2$$

**Properties:**
- $\lambda = 1$ when only one species
- $\lambda = 1/S$ when all species equally abundant
- Higher values indicate lower diversity (more dominance)

**Source:** Simpson (1949), Wikipedia — Diversity Index

### 4. Inverse Simpson Index

The effective number of equally abundant species (true diversity of order 2).

$$D = \frac{1}{\lambda} = \frac{1}{\sum_{i=1}^{S} p_i^2}$$

**Source:** Hill (1973), Wikipedia — Diversity Index

### 5. Pielou's Evenness (J)

Measures how evenly individuals are distributed among species.

$$J = \frac{H}{\ln(S)}$$

Where:
- $H$ = Shannon index
- $S$ = observed species count

**Properties:**
- $J \in [0, 1]$
- $J = 1$ when all species equally abundant
- Undefined when $S \leq 1$ (implementation returns 0)

**Source:** Pielou (1966), Wikipedia — Species Evenness

### 6. Chao1 Estimator

Estimates true species richness accounting for unobserved species.

$$\hat{S}_{Chao1} = S_{obs} + \frac{f_1^2}{2f_2}$$

Where:
- $f_1$ = number of singletons
- $f_2$ = number of doubletons

**Note:** Implementation uses simplified version ($Chao1 = S_{obs}$) as abundance proportions lack singleton/doubleton counts.

**Source:** Chao (1984)

## Algorithm

### Input

- `IReadOnlyDictionary<string, double> abundances`: Map of taxon names to relative abundances

### Output

- `AlphaDiversity` record containing all indices

### Steps

1. **Handle empty/null input**: Return all-zero result
2. **Filter zero abundances**: Remove taxa with $p_i \leq 0$
3. **Count observed species**: $S = $ filtered count
4. **Normalize abundances**: Ensure $\sum p_i = 1$
5. **Calculate Shannon**: $H = -\sum p_i \ln(p_i)$
6. **Calculate Simpson**: $\lambda = \sum p_i^2$
7. **Calculate Inverse Simpson**: $1/\lambda$ (or 0 if $\lambda = 0$)
8. **Calculate Pielou**: $J = H / \ln(S)$ (or 0 if $S \leq 1$)
9. **Set Chao1**: Simplified to $S_{obs}$

## Implementation

### Location

`Seqeron.Genomics.MetagenomicsAnalyzer.CalculateAlphaDiversity()`

### Return Type

```csharp
public readonly record struct AlphaDiversity(
    double ShannonIndex,
    double SimpsonIndex,
    double InverseSimpson,
    double Chao1Estimate,
    double ObservedSpecies,
    double PielouEvenness);
```

### Implementation Notes

1. Uses natural logarithm (ln) for Shannon index
2. Normalizes input abundances before calculation
3. Filters zero/negative abundances
4. Returns zero for all indices on empty input
5. Chao1 is simplified (equals ObservedSpecies)

## Edge Cases

| Scenario | Shannon | Simpson | InvSimpson | Pielou | ObservedSpecies |
|----------|---------|---------|------------|--------|-----------------|
| Empty input | 0 | 0 | 0 | 0 | 0 |
| Single species | 0 | 1.0 | 1.0 | 0 | 1 |
| Two equal (0.5, 0.5) | ln(2) ≈ 0.693 | 0.5 | 2.0 | 1.0 | 2 |
| Four equal (0.25 each) | ln(4) ≈ 1.386 | 0.25 | 4.0 | 1.0 | 4 |
| Highly uneven (0.99, 0.01) | ≈ 0.056 | ≈ 0.98 | ≈ 1.02 | < 1 | 2 |

## References

1. Shannon, C. E. (1948). A Mathematical Theory of Communication. Bell System Technical Journal.
2. Simpson, E. H. (1949). Measurement of Diversity. Nature, 163(4148), 688.
3. Hill, M. O. (1973). Diversity and Evenness: A Unifying Notation. Ecology, 54(2), 427-432.
4. Chao, A. (1984). Non-parametric estimation of the number of classes in a population. Scandinavian Journal of Statistics, 11, 265-270.
5. Whittaker, R. H. (1960). Vegetation of the Siskiyou Mountains. Ecological Monographs, 30, 279-338.
6. Wikipedia — Diversity Index: https://en.wikipedia.org/wiki/Diversity_index
7. Wikipedia — Alpha Diversity: https://en.wikipedia.org/wiki/Alpha_diversity

## Related Test Units

- [META-PROF-001](../../TestSpecs/META-PROF-001.md) — Taxonomic Profile (uses Shannon/Simpson)
- [META-BETA-001](./Beta_Diversity.md) — Beta Diversity
