# Genome Binning

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-BIN-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Genome binning groups assembled metagenomic contigs into bins intended to represent their genomes of origin, producing candidate metagenome-assembled genomes (MAGs). Classical binning methods use compositional signals such as GC content and oligonucleotide frequencies, often combined with coverage information, to cluster contigs from the same organism. The repository implementation follows that general pattern with a deterministic k-means workflow over GC content, tetranucleotide frequencies, and coverage. It also reports completeness and contamination as proxy scores rather than marker-gene-based estimates.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Metagenomic assemblies interleave contigs from many organisms, so binning attempts to recover organism-level groupings after assembly. Composition-based methods rely on the observation that genomes tend to maintain characteristic GC content and oligonucleotide usage signatures, while coverage adds an orthogonal abundance signal when contigs from the same organism share similar read depth patterns (Teeling et al., 2004).

### 2.2 Core Model

The theoretical basis for the current feature set is a joint representation of each contig by composition and abundance-related signals:

- GC content as the fraction of `G` and `C` bases in the contig.
- Tetranucleotide frequency (TNF) as the normalized frequency vector over the `4^4 = 256` possible 4-mers (Teeling et al., 2004).
- Coverage as a depth-derived abundance feature.

The implementation combines these signals with a composite distance described in the source code:

$$
d(a, b) = |GC_a - GC_b| + |Cov_a - Cov_b| + \frac{1 - r_{TNF}(a, b)}{2}
$$

where $r_{TNF}(a, b)$ is the Pearson correlation between the two TNF vectors. The repository then applies hard clustering by k-means over that feature space.

The document also cites MIMAG quality categories for interpreting MAG quality:

- High-quality: `> 90%` complete and `< 5%` contamination.
- Medium-quality: `>= 50%` complete and `< 10%` contamination.
- Low-quality: `< 50%` complete and `< 10%` contamination.

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Contigs from the same genome have similar GC content and tetranucleotide composition | Composition-based clustering may merge unrelated contigs or split one genome across bins |
| ASM-02 | Coverage carries organism-specific signal that helps separate bins | Coverage contributes noise instead of separation when samples or contigs do not preserve comparable abundance patterns |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each reported bin is a hard cluster of contigs rather than an overlapping assignment | The implementation uses k-means assignments, so each contig index belongs to one cluster |
| INV-02 | Reported completeness and contamination proxy values are clamped to `[0, 100]` | Both formulas use `Math.Min(..., 100)` in the implementation |
| INV-03 | Bins smaller than the configured minimum size are not emitted | The code skips clusters whose `totalLength` is below `minBinSize` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `contigs` | `IEnumerable<(string ContigId, string Sequence, double Coverage)>` | required | Input contigs with sequence and coverage data | Empty input yields no output bins |
| `numBins` | `int` | `10` | Maximum number of bins (`k`) for k-means clustering | Effective `k` becomes `min(numBins, contigCount)` |
| `minBinSize` | `double` | `500000` | Minimum total assembled length required for a bin to be reported | Clusters below this threshold are skipped |
| `expectedGenomeSize` | `double` | `4000000` | Expected genome size in base pairs for completeness estimation | Used only in the completeness proxy formula |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `BinId` | `string` | Generated bin identifier of the form `bin.<n>` |
| `ContigIds` | `IReadOnlyList<string>` | Contig identifiers assigned to the bin |
| `TotalLength` | `double` | Sum of sequence lengths of all contigs in the bin |
| `GcContent` | `double` | Mean GC content over the bin's contigs |
| `Coverage` | `double` | Mean raw coverage over the bin's contigs |
| `Completeness` | `double` | Proxy completeness percentage based on total length versus `expectedGenomeSize` |
| `Contamination` | `double` | Proxy contamination percentage based on within-bin GC standard deviation |
| `PredictedTaxonomy` | `string` | Taxonomy placeholder; the current implementation returns an empty string |

### 3.3 Preconditions and Validation

An empty input collection produces no bins. Coverage values are normalized internally to `[0, 1]` using the maximum observed coverage before clustering. Empty sequences produce zero GC content and no TNF entries through the helper methods. After clustering, any cluster whose total assembled length is less than `minBinSize` is filtered out before output.

## 4. Algorithm

### 4.1 High-Level Steps

1. Materialize the input contigs.
2. Compute GC content, TNF vector, and raw coverage for each contig.
3. Normalize coverage values by the maximum coverage observed in the input set.
4. Initialize `k` centroids by spreading selections across GC-sorted contigs.
5. Run k-means assignment and centroid updates until assignments stabilize or `50` iterations are reached.
6. Group contigs by final cluster assignment.
7. Discard clusters whose total sequence length is below `minBinSize`.
8. Emit `GenomeBin` records with mean GC content, mean coverage, completeness proxy, and contamination proxy.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The source code defines the composite distance as the sum of GC difference, normalized coverage difference, and TNF Pearson distance. TNF Pearson distance is `(1 - r) / 2`, which maps Pearson correlation into `[0, 1]`. Completeness is computed as `min(totalLength / expectedGenomeSize * 100, 100)`. Contamination is computed as `min(gcStdDev / 0.5 * 100, 100)`, where `0.5` is treated as the theoretical maximum standard deviation for GC fractions in `[0, 1]`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `BinContigs` | `O(n * k * i)` | `O(n * f)` | `n` = contigs, `k` = bins, `i` = k-means iterations, `f` = feature dimensions |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs)

- `MetagenomicsAnalyzer.BinContigs(IEnumerable<(string ContigId, string Sequence, double Coverage)>, int, double, double)`: Performs the repository's contig binning workflow and emits `GenomeBin` results.

### 5.2 Current Behavior

The repository implementation fixes the clustering method to k-means and uses deterministic centroid initialization by spreading across GC-sorted contigs. Coverage is normalized by the maximum input coverage, but average coverage in the output bin record is reported in the original scale. Completeness and contamination are proxies rather than marker-gene-derived scores. The `PredictedTaxonomy` field is currently emitted as an empty string.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Use of GC content, tetranucleotide frequency, and coverage as binning features.
- Pearson-correlation-based TNF comparison in the spirit of TETRA-style composition analysis.
- Hard clustering of contigs into bins.

**Intentionally simplified:**

- Bin quality uses `totalLength / expectedGenomeSize` as a completeness proxy instead of single-copy marker genes; **consequence:** completeness tracks assembled length rather than marker-gene recovery.
- Contamination uses within-bin GC standard deviation instead of duplicated marker genes; **consequence:** contamination is a composition-based proxy and not a CheckM-equivalent estimate.
- The implementation fixes the clustering method to k-means; **consequence:** density-based or hierarchical alternatives discussed in the document are not available in this API.

**Not implemented:**

- Marker-gene-based completeness and contamination estimation; **users should rely on:** external MAG quality tools such as CheckM rather than this class.
- Predicted taxonomy assignment for output bins; **users should rely on:** no current alternative in this method.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Completeness proxy by expected genome size | Deviation | Length-rich but incomplete bins can appear more complete than marker-gene-based tools would report | accepted | Controlled by the caller-supplied `expectedGenomeSize` parameter |
| 2 | Contamination proxy by GC variance | Assumption | Low-variance mixed bins and high-variance single-genome bins can be misinterpreted relative to marker-gene methods | accepted | The code documents this as a proxy tied to within-bin GC dispersion |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input collection | No bins are returned | The method exits early when the contig list is empty |
| `numBins` greater than the number of contigs | Effective `k` equals the number of contigs | The implementation sets `k = min(numBins, contigCount)` |
| Cluster total length below `minBinSize` | Cluster is omitted from the output | Post-clustering filter on `totalLength` |
| Bin with fewer than two contigs | `Contamination = 0` | The GC-variance helper returns `0` when fewer than two values are present |

### 6.2 Limitations

This implementation does not include marker-gene databases, lineage-specific models, or plasmid/genomic-island-specific handling. The cited literature notes that short-read metagenome binning can fail disproportionately for plasmids and genomic islands, so output bins should be interpreted as an algorithmic approximation rather than production-grade MAG curation.

## 7. Examples and Related Material

- [META-BIN-001](../../../tests/TestSpecs/META-BIN-001.md) documents the repository's genome-binning test specification.

## 8. References

1. Teeling, H., et al. 2004. TETRA: a web-service and a stand-alone program for the analysis and comparison of tetranucleotide usage patterns in DNA sequences. BMC Bioinformatics 5:163. doi:10.1186/1471-2105-5-163.
2. Parks, D. H., et al. 2014. Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research 25:1043-1055.
3. Maguire, F., et al. 2020. Metagenome-assembled genome binning methods with short reads disproportionately fail for plasmids and genomic Islands. Microbial Genomics 6(10). doi:10.1099/mgen.0.000436.
4. Wikipedia contributors. Binning (metagenomics). Wikipedia. https://en.wikipedia.org/wiki/Binning_(metagenomics)

