# META-BETA-001 Evidence

**Test Unit:** META-BETA-001 (Beta Diversity)  
**Algorithm:** CalculateBetaDiversity  
**Date Collected:** 2026-02-04  

## Authoritative Sources

### Primary Sources

1. **Whittaker (1960)** - "Vegetation of the Siskiyou Mountains, Oregon and California"
   - Ecological Monographs 30(3): 279-338
   - Original definition of beta diversity concept
   - Source: [Wikipedia Beta Diversity](https://en.wikipedia.org/wiki/Beta_diversity)

2. **Bray & Curtis (1957)** - "An Ordination of the Upland Forest Communities of Southern Wisconsin"
   - Ecological Monographs 27(4): 325-349
   - Original Bray-Curtis dissimilarity formula
   - Source: [Wikipedia Bray-Curtis dissimilarity](https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity)

3. **Jaccard (1901)** - "Étude comparative de la distribution florale dans une portion des Alpes et des Jura"
   - Bulletin de la Société vaudoise des sciences naturelles 37(142): 547-579
   - Original Jaccard index definition
   - Source: [Wikipedia Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index)

### Key Algorithmic Principles

#### Bray-Curtis Dissimilarity
- **Formula**: `BC_jk = 1 - (2 * C_jk) / (S_j + S_k)`
- **Where**: C_jk = sum of lesser counts for each species, S_j/S_k = total specimens in each site
- **Range**: [0, 1] where 0 = identical composition, 1 = no shared species
- **Properties**: NOT a true distance metric (violates triangle inequality)

#### Jaccard Distance  
- **Formula**: `J_distance = 1 - |A ∩ B| / |A ∪ B|`
- **Alternative**: `J_distance = (|A| + |B| - 2|A ∩ B|) / (|A| + |B| - |A ∩ B|)`
- **Range**: [0, 1] where 0 = identical sets, 1 = no shared elements
- **Properties**: True distance metric (satisfies triangle inequality)

## Edge Cases and Testing Requirements

### Critical Edge Cases (Evidence-Based)

1. **Identical Samples**
   - Both Bray-Curtis and Jaccard should equal 0
   - Source: Mathematical definition in Wikipedia sources

2. **Completely Disjoint Samples**
   - Both Bray-Curtis and Jaccard should equal 1
   - Source: Mathematical definition in Wikipedia sources

3. **Empty Sample Handling**
   - Mathematical formula undefined when both samples empty (division by zero)
   - Requires implementation-specific handling
   - Source: Mathematical analysis of formulas

4. **Single-Species Dominance**
   - Test behavior when one species has 100% abundance
   - Source: Ecological interpretation from Wikipedia

### Algorithmic Invariants

1. **Range Constraints**
   - Bray-Curtis ∈ [0, 1]
   - Jaccard Distance ∈ [0, 1]
   - Source: Mathematical definitions

2. **Symmetry Property**
   - Distance(A, B) = Distance(B, A)
   - Source: Mathematical definitions

3. **Identity Property**
   - Distance(A, A) = 0 for any valid sample A
   - Source: Mathematical definitions

## Implementation-Specific Notes

Based on current implementation in MetagenomicsAnalyzer.cs:
- Returns BetaDiversity record with both metrics
- Uses abundance data (continuous values, not just presence/absence)
- Includes UniFracDistance field (set to 0, requires phylogenetic tree)
- Handles shared species counting

## Test Dataset Requirements

**ASSUMPTION**: No published standard test datasets found in sources for beta diversity validation. Implementation testing should use constructed examples with known mathematical outcomes.

## References

1. Whittaker, R.H. (1960). Vegetation of the Siskiyou Mountains, Oregon and California. Ecological Monographs, 30(3), 279-338.
2. Bray, J.R. & Curtis, J.T. (1957). An ordination of the upland forest communities of southern Wisconsin. Ecological Monographs, 27(4), 325-349.
3. Jaccard, P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société vaudoise des sciences naturelles, 37(142), 547-579.
4. Wikipedia contributors. Beta diversity. Wikipedia. Retrieved 2026-02-04.
5. Wikipedia contributors. Bray-Curtis dissimilarity. Wikipedia. Retrieved 2026-02-04.  
6. Wikipedia contributors. Jaccard index. Wikipedia. Retrieved 2026-02-04.