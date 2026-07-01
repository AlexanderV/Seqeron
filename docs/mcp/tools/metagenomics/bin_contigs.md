# bin_contigs

Cluster contigs into metagenome-assembled genomes (MAGs).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `bin_contigs` |
| **Method ID** | `MetagenomicsAnalyzer.BinContigs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Bins contigs by k-means clustering over three compositional/coverage features:
GC content, normalized coverage, and tetranucleotide frequency (TETRA, Teeling et al. 2004,
compared via Pearson distance). Centroids are seeded deterministically by spreading across
GC-sorted contigs, so the result is reproducible.

For each bin with total length `≥ minBinSize` it reports:

- **completeness** `= min(totalLength / expectedGenomeSize × 100, 100)`
- **contamination** `= min(stddev(GC) / 0.5 × 100, 100)` (within-bin GC variance, Parks et al. 2014)
- mean GC content, mean coverage, and member contig ids.

Bins below `minBinSize` are dropped.

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L855](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L855)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `contigs` | ContigInput[] | Yes | `{ contigId, sequence, coverage }` entries. |
| `numBins` | integer | No | Maximum bins / k (default 10). |
| `minBinSize` | number | No | Minimum reported bin length in bp (default 500000). |
| `expectedGenomeSize` | number | No | Expected genome size in bp (default 4000000). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].binId` | string | Bin id (`bin.1`, `bin.2`, …). |
| `items[].contigIds` | string[] | Member contig ids. |
| `items[].totalLength` | number | Total bin length in bp. |
| `items[].gcContent` | number | Mean within-bin GC fraction [0,1]. |
| `items[].coverage` | number | Mean within-bin coverage. |
| `items[].completeness` | number | Completeness percentage. |
| `items[].contamination` | number | Contamination percentage. |
| `items[].predictedTaxonomy` | string | Reserved (empty). |

## Errors

None. Empty input or all-below-`minBinSize` returns an empty list.

## Example

**User Prompt:**
> Bin these 10 contigs (200 kb each, ~50% GC, coverage 15) into a single MAG.

With `numBins = 1`, `minBinSize = 100000`:

**Response:**
```json
{
  "items": [
    { "binId": "bin.1", "totalLength": 2000000, "gcContent": 0.5, "coverage": 15.0, "completeness": 50.0, "contamination": 0.0, "predictedTaxonomy": "" }
  ]
}
```

Completeness = 2 Mb / 4 Mb × 100 = 50 %; uniform GC → contamination 0.

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `BinContigs`
- Teeling H. et al. (2004) BMC Bioinformatics 5:163 (TETRA); Parks D.H. et al. (2014) Genome Res. 25:1043 (CheckM).

## See Also

- [taxonomic_profile](taxonomic_profile.md)
- [predict_functions](predict_functions.md)
