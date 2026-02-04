# Beta Diversity Analysis

## Overview

Beta diversity measures the compositional dissimilarity between two ecological communities or samples. It quantifies species turnover and complementarity between sites, providing insight into biodiversity patterns across landscapes.

## Mathematical Foundation

### Bray-Curtis Dissimilarity

The Bray-Curtis dissimilarity coefficient (Bray & Curtis, 1957) is defined as:

```
BC_jk = 1 - (2 * C_jk) / (S_j + S_k)
```

Where:
- `C_jk` = sum of the lesser counts for each species found in both sites
- `S_j` = total number of specimens in site j  
- `S_k` = total number of specimens in site k

**Properties:**
- Range: [0, 1]
- 0 = identical composition
- 1 = no shared species
- Not a true distance metric (violates triangle inequality)

### Jaccard Distance

The Jaccard distance is the complement of the Jaccard similarity index:

```
J_distance = 1 - |A ∩ B| / |A ∪ B|
```

For presence/absence data:
```
J_distance = (unique_to_A + unique_to_B) / (shared + unique_to_A + unique_to_B)
```

**Properties:**
- Range: [0, 1]  
- 0 = identical species sets
- 1 = no shared species
- True distance metric (satisfies triangle inequality)

## Implementation Notes

### Current Implementation

The `MetagenomicsAnalyzer.CalculateBetaDiversity` method:

1. Accepts two abundance dictionaries (species → count)
2. Calculates species overlap statistics
3. Computes both Bray-Curtis and Jaccard distances
4. Returns structured result with metrics and species counts

### Key Considerations

- **Abundance vs. Presence**: Implementation uses abundance data for Bray-Curtis, converts to presence/absence for Jaccard
- **Zero Handling**: Species with zero abundance are treated as absent
- **Normalization**: No explicit normalization applied to input abundances
- **UniFrac Distance**: Field provided but requires phylogenetic tree data (not implemented)

## Applications

- Comparing microbial community composition between samples
- Assessing habitat similarity in ecological studies
- Monitoring temporal changes in community structure
- Evaluating treatment effects on species composition

## Limitations

- Sensitive to sampling effort differences
- May be dominated by rare species in presence/absence calculations
- Does not incorporate phylogenetic relationships (standard metrics)
- Assumes equal weighting of all species

## References

1. Bray, J.R. & Curtis, J.T. (1957). An ordination of the upland forest communities of southern Wisconsin. Ecological Monographs, 27(4), 325-349.
2. Jaccard, P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société vaudoise des sciences naturelles, 37(142), 547-579.
3. Whittaker, R.H. (1960). Vegetation of the Siskiyou Mountains, Oregon and California. Ecological Monographs, 30(3), 279-338.