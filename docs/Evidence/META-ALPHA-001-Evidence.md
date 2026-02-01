# META-ALPHA-001: Alpha Diversity — Evidence Document

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | META-ALPHA-001 |
| **Area** | Metagenomics |
| **Canonical Methods** | `MetagenomicsAnalyzer.CalculateAlphaDiversity` |
| **Date** | 2026-02-02 |

## Evidence Sources

### Primary Sources

1. **Wikipedia — Diversity Index**
   - URL: https://en.wikipedia.org/wiki/Diversity_index
   - Accessed: 2026-02-02
   - Content: Comprehensive coverage of Shannon, Simpson, and other diversity indices with formulas

2. **Wikipedia — Alpha Diversity**
   - URL: https://en.wikipedia.org/wiki/Alpha_diversity
   - Accessed: 2026-02-02
   - Content: Definition of alpha diversity as mean species diversity at local scale

3. **Wikipedia — Species Richness**
   - URL: https://en.wikipedia.org/wiki/Species_richness
   - Accessed: 2026-02-02
   - Content: Species richness definition and relationship to diversity

4. **Wikipedia — Species Evenness**
   - URL: https://en.wikipedia.org/wiki/Species_evenness
   - Accessed: 2026-02-02
   - Content: Pielou evenness definition and Shannon formula

### Foundational References (cited in Wikipedia)

5. **Shannon (1948)**: "A Mathematical Theory of Communication"
   - Shannon entropy formula: H' = −Σ pᵢ ln(pᵢ)

6. **Simpson (1949)**: "Measurement of Diversity"
   - Simpson index formula: λ = Σ pᵢ²

7. **Chao (1984)**: "Non-parametric estimation of the number of classes in a population"
   - Chao1 richness estimator (requires singleton/doubleton counts)

8. **Hill (1973)**: "Diversity and evenness: a unifying notation and its consequences"
   - Unified framework for diversity indices

9. **Pielou (1966)**: Evenness as J = H / ln(S)
   - Pielou's evenness index

## Extracted Formulas

### Shannon Index (H)

From Wikipedia Diversity Index:

$$H' = -\sum_{i=1}^{R} p_i \ln(p_i)$$

Where:
- R = species richness (number of species)
- pᵢ = proportional abundance of species i

**Properties:**
- H = 0 when only one species (no uncertainty)
- H = ln(R) when all species equally abundant (maximum entropy)
- Uses natural logarithm (nats)

### Simpson Index (λ)

From Wikipedia Diversity Index:

$$\lambda = \sum_{i=1}^{R} p_i^2$$

**Properties:**
- λ = 1 when only one species (all probability concentrated)
- λ = 1/R when all species equally abundant
- Represents probability that two randomly chosen individuals belong to same species

### Inverse Simpson Index

$$\frac{1}{\lambda} = \frac{1}{\sum_{i=1}^{R} p_i^2}$$

Equals true diversity of order 2 (effective number of species).

### Pielou's Evenness (J)

From Wikipedia Species Evenness and Diversity Index:

$$J = \frac{H}{\ln(S)}$$

Where:
- H = Shannon index
- S = species richness (observed species count)

**Properties:**
- J ∈ [0, 1]
- J = 1 when all species equally abundant
- J = 0 when one species completely dominates
- Undefined (or 0 by convention) when S ≤ 1

### Chao1 Estimator

From Chao (1984):

$$\hat{S}_{Chao1} = S_{obs} + \frac{f_1^2}{2f_2}$$

Where:
- S_obs = observed species count
- f₁ = number of singletons (species observed exactly once)
- f₂ = number of doubletons (species observed exactly twice)

**Note:** Implementation currently uses simplified version (Chao1 = ObservedSpecies) as singleton/doubleton counts are not available from abundance data.

## Edge Cases Documented

### From Wikipedia/Shannon Theory

1. **Empty input**: No species → all metrics = 0
2. **Single species**: H = 0, λ = 1, J undefined (or 0)
3. **Perfect evenness**: H = ln(S), λ = 1/S, J = 1
4. **Zero abundances**: Must be filtered out (ln(0) undefined)
5. **Unnormalized abundances**: Must normalize to sum = 1 before calculation

### Corner Cases

1. **All zero abundances**: Should behave like empty input
2. **Negative abundances**: Invalid input (undefined behavior)
3. **Null input**: Should return zero/empty result
4. **Single species with zero abundance filtered**: Returns S = 0

## Test Data from Sources

### From Wikipedia Diversity Index Examples

| Distribution | Shannon H (ln) | Simpson λ |
|--------------|----------------|-----------|
| Single species (p=1.0) | 0 | 1.0 |
| Two equal species (0.5, 0.5) | ln(2) ≈ 0.693 | 0.5 |
| Four equal species (0.25 each) | ln(4) ≈ 1.386 | 0.25 |

### Evenness Test Cases

| Distribution | Observed S | Shannon H | Pielou J |
|--------------|------------|-----------|----------|
| (0.5, 0.5) | 2 | ln(2) | 1.0 |
| (0.9, 0.1) | 2 | 0.325 | 0.469 |
| (0.25, 0.25, 0.25, 0.25) | 4 | ln(4) | 1.0 |

## Implementation Notes

### Current Implementation (MetagenomicsAnalyzer.cs)

The implementation:
1. Filters out zero/negative abundances
2. Normalizes abundances to sum = 1
3. Calculates Shannon using natural log
4. Calculates Simpson as Σpᵢ²
5. Calculates Inverse Simpson as 1/λ
6. Calculates Pielou's J = H / ln(S) for S > 1, else 0
7. Uses simplified Chao1 = ObservedSpecies (no singleton/doubleton data)

### Deviation from Theory

**Chao1 Simplification**: The implementation returns `Chao1 = ObservedSpecies` because the input is abundance proportions, not raw counts with singleton/doubleton information. This is explicitly noted in code comments.

## Assumptions

1. **ASSUMPTION**: Input abundances are relative proportions (sum to 1.0 or will be normalized)
2. **ASSUMPTION**: Zero abundances indicate absent species and should be filtered
3. **ASSUMPTION**: When S ≤ 1, Pielou evenness = 0 (not undefined)
4. **ASSUMPTION**: Null input treated as empty input (returns all zeros)
