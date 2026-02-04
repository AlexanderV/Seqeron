# META-BIN-001: Genome Binning - Evidence Document

## Sources

### Primary Sources

1. **Wikipedia - Binning (metagenomics)**
   - URL: https://en.wikipedia.org/wiki/Binning_(metagenomics)
   - Key concepts:
     - Metagenomic binning groups assembled contigs to their separate genomes of origin
     - Results are called Metagenome Assembled Genomes (MAGs)
     - Methods use compositional features (GC-content, tetranucleotide frequencies)
     - Methods also use coverage across samples
     - Can employ supervised or unsupervised classification
   - Referenced algorithms: TETRA, MEGAN, Phylopythia, SOrt-ITEMS, MetaBAT2, MaxBin2, CONCOCT, DAS Tool

2. **TETRA (Teeling et al., 2004)**
   - URL: https://doi.org/10.1186/1471-2105-5-163
   - Tetranucleotide usage patterns for genomic fragment classification
   - 256 possible tetramers (4^4)
   - Z-scores indicate over/under-representation vs expected frequencies

3. **CheckM (Parks et al., 2014)**
   - URL: http://genome.cshlp.org/content/25/7/1043
   - Tool for assessing genome quality from metagenomes
   - Estimates completeness using single-copy marker genes
   - Estimates contamination using duplicated marker genes
   - Reference: Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW

4. **Maguire et al. (2020) - MAG binning limitations**
   - URL: https://doi.org/10.1099/mgen.0.000436
   - Key findings:
     - 82-94% chromosomes recovered correctly
     - Only 38-44% of GIs recovered
     - Only 1-29% of plasmid sequences recovered
     - Variable copy number and composition problematic for binning

### Documented Binning Features

1. **Compositional Features**
   - GC content: genome-specific signature
   - Tetranucleotide frequencies: species-specific patterns
   - Codon usage patterns

2. **Coverage Features**
   - Read depth across samples
   - Differential coverage binning

3. **Quality Metrics (from CheckM)**
   - Completeness: proportion of expected marker genes present
   - Contamination: proportion of duplicated marker genes
   - Strain heterogeneity: indicator of multiple strains

### Edge Cases from Literature

1. **Empty input**: Should return empty results
2. **Single contig**: May form its own bin or be filtered by size
3. **Minimum bin size**: Filter bins below threshold
4. **Similar GC content across organisms**: Can cause incorrect binning
5. **Variable coverage regions**: Cause assembly and binning issues
6. **High contamination bins**: Multiple genomes in one bin

### Test Dataset Design (Based on Literature)

1. **Pure genomes** (high GC vs low GC): Should separate into different bins
2. **Similar GC content**: May cluster together
3. **Variable coverage**: Test coverage feature contribution
4. **Minimum size filtering**: Bins below threshold excluded

## Implementation Analysis

### Current Implementation: `MetagenomicsAnalyzer.BinContigs`

```csharp
public static IEnumerable<GenomeBin> BinContigs(
    IEnumerable<(string ContigId, string Sequence, double Coverage)> contigs,
    int numBins = 10,
    double minBinSize = 500000)
```

**Algorithm (simplified k-means based on GC):**
1. Calculate features: GC content, tetranucleotide frequency per contig
2. Assign to bins based on GC content: `binIndex = (int)(GcContent * numBins) % numBins`
3. Filter bins below minimum size
4. Estimate completeness based on total length (simplified)
5. Estimate contamination based on GC variance

**Output: GenomeBin record**
- BinId, ContigIds, TotalLength, GcContent, Coverage
- Completeness, Contamination, PredictedTaxonomy

### Deviations from Literature

1. **ASSUMPTION**: Completeness estimated by length ratio to 2Mbp genome
   - Real tools use marker gene detection (CheckM)
   - Simplified for this implementation

2. **ASSUMPTION**: Contamination estimated from GC variance
   - Real tools detect duplicated single-copy marker genes
   - High GC variance suggests multiple organisms

3. **Tetranucleotide frequency calculated but not used in binning**
   - Implementation calculates it but bins only on GC
   - Full implementation would use multi-feature clustering

## Invariants

1. Result bin IDs must be unique
2. Completeness must be in range [0, 100]
3. Contamination must be in range [0, 100]
4. GC content must be in range [0, 1]
5. Total length must equal sum of contig lengths in bin
6. Empty input returns empty output
7. Bins below minBinSize are excluded

## Test Categories

### Must Tests (Evidence-backed)
1. Empty input handling
2. Basic binning returns non-null
3. GC-based separation (high vs low GC contigs)
4. MinBinSize filtering
5. Completeness range [0, 100]
6. Contamination range [0, 100]
7. Unique bin IDs
8. Correct contig grouping
9. TotalLength accuracy

### Should Tests
1. Coverage values preserved
2. GC content averaging in bins
3. Multiple bins from diverse input

### Could Tests
1. Large dataset performance
2. Edge case with all contigs same GC

---

**Document Version:** 1.0
**Created:** 2026-02-04
**Test Unit:** META-BIN-001
