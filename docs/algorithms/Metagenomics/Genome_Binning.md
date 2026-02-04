# Genome Binning

## Overview

Genome binning is the computational process of grouping assembled contigs from metagenomic data and assigning them to their separate genomes of origin. The resulting grouped contigs are called **Metagenome Assembled Genomes (MAGs)**.

## Algorithm

### Compositional Features

Binning methods exploit organism-specific characteristics of DNA:

1. **GC Content**: The ratio of guanine and cytosine bases to total bases. Different organisms have characteristic GC content signatures that remain relatively constant across their genomes.

2. **Tetranucleotide Frequencies**: There are 256 (4⁴) possible 4-mer combinations. The frequency distribution of these tetramers is species-specific and can be used to cluster contigs from the same organism.

### Coverage Features

Read depth across samples provides additional signal:
- Contigs from the same genome should have correlated coverage patterns
- Differential coverage across samples improves separation

### Clustering Approach

Common approaches include:
- **K-means clustering**: Partitions contigs into k clusters based on feature similarity
- **DBSCAN**: Density-based clustering that can identify arbitrary-shaped clusters
- **Hierarchical clustering**: Builds nested clusters

## Quality Metrics

### Completeness

Proportion of expected single-copy marker genes present in a bin. Real implementations (CheckM) use:
- Sets of universal single-copy genes
- Lineage-specific marker gene sets
- Completeness = (observed markers / expected markers) × 100

### Contamination

Indicator of multiple genomes in a single bin:
- Measured by duplicated single-copy marker genes
- Contamination = (duplicated markers / total markers) × 100

### MIMAG Standards

Minimum Information about a Metagenome-Assembled Genome:
- **High-quality**: >90% complete, <5% contamination
- **Medium-quality**: ≥50% complete, <10% contamination
- **Low-quality**: <50% complete, <10% contamination

## Complexity

- **Time Complexity**: O(n × k × i) where n = contigs, k = bins, i = iterations
- **Space Complexity**: O(n × f) where f = feature dimensions

## Implementation Notes

### Current Implementation

The `MetagenomicsAnalyzer.BinContigs` method provides a simplified binning implementation:

1. Calculates GC content and tetranucleotide frequencies for each contig
2. Assigns contigs to bins based on GC content (simplified k-means)
3. Filters bins below minimum size threshold
4. Estimates completeness from total length (approximation)
5. Estimates contamination from GC variance within bin

### Simplifications

1. **Completeness**: Estimated by ratio of bin length to 2Mbp (average bacterial genome)
   - Real tools use marker gene detection

2. **Contamination**: Estimated from GC variance
   - Real tools detect duplicated single-copy marker genes

3. **Clustering**: Uses only GC content
   - Real tools use multi-dimensional feature space

## References

1. Teeling H, et al. (2004). "TETRA: a web-service and a stand-alone program for the analysis and comparison of tetranucleotide usage patterns in DNA sequences." BMC Bioinformatics. doi:10.1186/1471-2105-5-163

2. Parks DH, et al. (2014). "Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome Research 25:1043-1055.

3. Maguire F, et al. (2020). "Metagenome-assembled genome binning methods with short reads disproportionately fail for plasmids and genomic Islands." Microbial Genomics 6(10). doi:10.1099/mgen.0.000436

4. Wikipedia: Binning (metagenomics). https://en.wikipedia.org/wiki/Binning_(metagenomics)

---

**Document Version:** 1.0
**Related Test Unit:** META-BIN-001
