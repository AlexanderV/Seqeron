# Taxonomic Profile Generation

## Overview

| Property | Value |
|----------|-------|
| **Algorithm** | Taxonomic Profile Generation |
| **Category** | Metagenomics |
| **Time Complexity** | O(n) where n = number of classifications |
| **Space Complexity** | O(t) where t = unique taxa |

## Description

Taxonomic profiling aggregates individual read classifications into a community-level summary showing the relative abundances of taxa at multiple taxonomic ranks. This is a fundamental analysis in metagenomics that transforms raw classification data into interpretable ecological metrics.

## Algorithm

### Input

- Collection of `TaxonomicClassification` records, each containing:
  - ReadId
  - Taxonomic assignments at each rank (kingdom, phylum, class, order, family, genus, species)
  - Confidence score

### Output

- `TaxonomicProfile` containing:
  - Abundance distributions at each taxonomic rank
  - Diversity metrics (Shannon, Simpson)
  - Read statistics (total, classified)

### Steps

1. **Count reads per taxon** at each taxonomic rank
2. **Exclude unclassified reads** from abundance denominators
3. **Compute relative abundances**: `abundance = count / total_classified`
4. **Calculate diversity metrics** from species-level abundances

## Mathematical Formulas

### Relative Abundance

For taxon $i$ with count $c_i$:

$$
\text{abundance}_i = \frac{c_i}{\sum_{j=1}^{n} c_j}
$$

where the sum excludes unclassified reads.

### Shannon Diversity Index

$$
H = -\sum_{i=1}^{S} p_i \ln(p_i)
$$

where:
- $S$ = number of species
- $p_i$ = relative abundance of species $i$

**Properties:**
- $H = 0$ when only one species present (no uncertainty)
- $H = \ln(S)$ when all species equally abundant (maximum entropy)

### Simpson Concentration Index

$$
D = \sum_{i=1}^{S} p_i^2
$$

**Properties:**
- $D = 1$ when only one species present
- $D = 1/S$ when all species equally abundant

## Implementation Notes

### Seqeron.Genomics Implementation

Located in `MetagenomicsAnalyzer.GenerateTaxonomicProfile()`:

```csharp
public static TaxonomicProfile GenerateTaxonomicProfile(
    IEnumerable<TaxonomicClassification> classifications)
```

**Key behaviors:**
1. Filters out reads where `Kingdom == "Unclassified"`
2. Empty taxonomic rank values excluded from rank-specific abundance maps
3. Diversity metrics computed from species-level abundances
4. Uses natural logarithm for Shannon index

### Edge Cases

| Case | Behavior |
|------|----------|
| Empty input | Returns profile with TotalReads=0, empty abundances |
| All unclassified | ClassifiedReads=0, empty abundance maps, Shannon=Simpson=0 |
| Single species | Shannon=0 (no uncertainty), Simpson=1.0 |
| Missing rank values | Excluded from that rank's abundance map |

## Sources

1. **Shannon CE (1948)**. "A Mathematical Theory of Communication." Bell System Technical Journal.
2. **Simpson EH (1949)**. "Measurement of Diversity." Nature.
3. **Segata N et al. (2012)**. "Metagenomic microbial community profiling using unique clade-specific marker genes." Nature Methods. DOI: 10.1038/nmeth.2066
4. **Wikipedia**. "Metagenomics" — https://en.wikipedia.org/wiki/Metagenomics

## Related Algorithms

- [Taxonomic Classification](./Taxonomic_Classification.md) — Upstream algorithm that produces classifications
- [Alpha Diversity](./Alpha_Diversity.md) — Detailed diversity analysis
- [Beta Diversity](./Beta_Diversity.md) — Cross-sample comparison
