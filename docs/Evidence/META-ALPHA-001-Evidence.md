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

When f₂ = 0, the bias-corrected form is used:

$$\hat{S}_{Chao1} = S_{obs} + \frac{f_1(f_1 - 1)}{2}$$

When f₁ = 0 (no singletons), Chao1 = S_obs.
For proportional data (non-integer values), singletons/doubletons are undefined, so Chao1 = S_obs.

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

### Implementation (MetagenomicsAnalyzer.cs)

The implementation matches all source formulas:
1. Filters out zero/negative abundances (ln(0) undefined — mathematical requirement)
2. Normalizes abundances to proportions summing to 1 (standard ecological practice)
3. Shannon H' = −Σ pᵢ ln(pᵢ) using natural logarithm — per Shannon (1948)
4. Simpson λ = Σ pᵢ² — per Simpson (1949)
5. Inverse Simpson = 1/λ = ²D (effective number of species) — per Hill (1973)
6. Pielou J = H / ln(S) for S > 1; J = 0 for S ≤ 1 — per Pielou (1966); ln(1) = 0 makes J mathematically undefined when S = 1, J = 0 is the standard ecological convention
7. Chao1 = S_obs + f₁²/(2·f₂) for integer count data with f₂ > 0; bias-corrected form S_obs + f₁·(f₁−1)/2 when f₂ = 0; S_obs for proportional (non-integer) data — per Chao (1984)
8. Null/empty input → all metrics = 0 (standard defensive boundary handling)

## Deviations and Assumptions

None. All formulas match external sources exactly.
