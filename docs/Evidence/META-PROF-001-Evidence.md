# Evidence: META-PROF-001 — Taxonomic Profile

## Overview

| Field | Value |
|-------|-------|
| Test Unit ID | META-PROF-001 |
| Algorithm | Taxonomic Profile Generation |
| Area | Metagenomics |
| Methods | `GenerateTaxonomicProfile` |
| Date | 2026-02-01 |

---

## Authoritative Sources

### Primary Sources

1. **Wikipedia — Metagenomics**
   - URL: https://en.wikipedia.org/wiki/Metagenomics
   - Relevance: Defines taxonomic profiling, relative abundance computation
   - Key excerpts:
     - "MetaPhlAn and AMPHORA are methods based on unique clade-specific markers for estimating organismal relative abundances"
     - "Metagenomes are most commonly compared by analyzing taxonomic composition—such as differences in normalized species or genus abundance between groups"
     - Taxonomic profiling aggregates classification results into abundance distributions

2. **Wikipedia — Relative Abundance Distribution**
   - URL: https://en.wikipedia.org/wiki/Relative_abundance_distribution
   - Relevance: Defines relative abundance in ecological/metagenomic contexts
   - Key excerpts:
     - "The SAD (Species Abundance Distribution) describes the relationship between the number of species observed as a function of their observed abundance"
     - Abundance is typically normalized as fractions summing to 1.0

3. **MetaPhlAn Documentation**
   - URL: https://huttenhower.sph.harvard.edu/metaphlan/
   - Relevance: Canonical metagenomic profiler
   - Key concepts:
     - Generates taxonomic profiles at multiple taxonomic levels (kingdom, phylum, class, order, family, genus, species)
     - Relative abundance expressed as percentages or fractions

4. **Segata et al. (2012)** — Primary Reference
   - Citation: Segata N, et al. Metagenomic microbial community profiling using unique clade-specific marker genes. Nature Methods 2012.
   - DOI: 10.1038/nmeth.2066
   - Relevance: MetaPhlAn paper describing taxonomic profile methodology

---

## Algorithm Characteristics

### Taxonomic Profile Definition

A taxonomic profile aggregates classified reads into:
1. **Counts per taxon** at each taxonomic rank (kingdom, phylum, genus, species)
2. **Relative abundances** (counts normalized to sum to 1.0 for classified reads)
3. **Diversity metrics** (Shannon, Simpson) computed from abundance distribution

### Abundance Calculation

Per Wikipedia/MetaPhlAn:
- **Absolute abundance**: Count of reads classified to each taxon
- **Relative abundance**: Fraction = count / total_classified_reads
- "Unclassified" reads excluded from abundance denominators

### Key Formulas

| Metric | Formula | Source |
|--------|---------|--------|
| Relative abundance | count(taxon) / Σcount(all classified taxa) | MetaPhlAn |
| Shannon diversity | H = -Σ(pᵢ × ln(pᵢ)) | Shannon (1948), Wikipedia |
| Simpson diversity | D = Σpᵢ² | Simpson (1949), Wikipedia |

---

## Test Datasets from Sources

### Documented Test Cases

| Case | Source | Expected Behavior |
|------|--------|-------------------|
| Empty classification list | Standard robustness | TotalReads=0, ClassifiedReads=0, empty abundances |
| All Unclassified | MetaPhlAn output format | ClassifiedReads=0, no abundance entries |
| Single taxon | Trivial case | Abundance = 1.0 for that taxon |
| Uniform distribution | Shannon/Simpson definition | Shannon = ln(n), Simpson = 1/n |
| Highly skewed distribution | Diversity theory | Low Shannon, high Simpson |

---

## Edge Cases and Corner Cases

### From Literature

1. **Empty Input**
   - Empty classification list → empty profile
   - TotalReads = 0, ClassifiedReads = 0

2. **All Unclassified**
   - No reads classified → empty abundance maps
   - ClassifiedReads = 0

3. **Missing Taxonomic Ranks**
   - Some reads may lack genus/species information
   - Empty strings filtered from abundance maps

4. **Diversity with Single Species**
   - Shannon = 0 (no uncertainty)
   - Simpson = 1.0 (certainty of picking same species)

5. **Diversity Normalization**
   - Abundances should sum to 1.0 before diversity calculation
   - Handle unnormalized input gracefully

---

## Implementation-Specific Notes

### Current Seqeron.Genomics Implementation

From `MetagenomicsAnalyzer.cs` (lines 229-280):

1. **Input**: `IEnumerable<TaxonomicClassification>`
2. **Output**: `TaxonomicProfile` record containing:
   - Abundance maps at four ranks (kingdom, phylum, genus, species)
   - Diversity metrics (Shannon, Simpson)
   - Read counts (TotalReads, ClassifiedReads)

3. **Behavior Notes**:
   - Excludes "Unclassified" from abundance calculation
   - Filters empty strings from rank-specific abundances
   - Diversity computed from species-level abundances

---

## Invariants

| Invariant | Description |
|-----------|-------------|
| Sum invariant | Σ(abundance values) ≈ 1.0 at each rank (for classified reads) |
| Count invariant | ClassifiedReads ≤ TotalReads |
| Shannon bounds | Shannon ≥ 0 |
| Simpson bounds | 0 ≤ Simpson ≤ 1.0 |
| Consistency | ClassifiedReads = Σ(counts at any rank) |

---

## Testing Methodology

### Recommended Approach

Per diversity metrics literature:

1. **Verify abundance normalization**: Sum of abundances = 1.0
2. **Verify count consistency**: ClassifiedReads matches abundance denominators
3. **Verify diversity formulas**: Use known distributions with calculable expected values
4. **Edge case coverage**: Empty input, single taxon, all unclassified

---

## ASSUMPTIONS

1. **Shannon uses natural log**: Implementation uses `Math.Log()` (natural log), consistent with ecological convention
2. **Simpson is concentration index**: Implementation returns Σpᵢ² (concentration), not 1-D (diversity)
3. **Empty abundances yield Shannon=0**: No species → no entropy (ASSUMPTION based on limit behavior)

---

## References

1. Shannon CE (1948). "A Mathematical Theory of Communication." Bell System Technical Journal.
2. Simpson EH (1949). "Measurement of Diversity." Nature.
3. Segata N et al. (2012). "Metagenomic microbial community profiling..." Nature Methods.
4. Wikipedia contributors. "Metagenomics." Wikipedia, The Free Encyclopedia.
5. Wikipedia contributors. "Relative abundance distribution." Wikipedia, The Free Encyclopedia.
