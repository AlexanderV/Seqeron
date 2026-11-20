# nucmer_clusters

MUMmer nucmer anchor clustering (mummer | mgaps).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `nucmer_clusters` |
| **Method ID** | `GenomicAnalyzer.FindNucmerClusters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The anchor-clustering stage of MUMmer 3.23 nucmer = `mummer -mum|-mumreference|-maxmatch -b -l L -n | mgaps -l c -s g -d D -f F` (700 related pairs, both strands, all three match modes equal; 758 further mgaps inputs equal): the exact matches (mode 'mum' = maximal unique matches in both sequences, 'mumreference' = unique in the reference, 'maxmatch' = all maximal matches; minMatchLength = -l; only A/C/G/T match) are clustered: matches are sorted by their start in the query, matches inside a repeat are removed and overlapping matches of one diagonal merged, neighbours whose diagonals differ by at most max(diagonalDifference, diagonalFraction * separation) within maxSeparation are joined, and each component is cut into chains of maximum score (match length - overlap - diagonal difference); a chain is reported when its length (sum of the match lengths, or the extent with useExtent) reaches minClusterLength. Coordinates are 1-based as printed by mummer / mgaps; on the reverse strand the query start is a position in the reverse-complemented query; a match's overlap is the number of bases cut off its start because it overlaps the previous match. This is the stage that decides which anchors belong together; the banded extension (postnuc) and delta-filter are not part of this tool.

## Core Documentation Reference

- Source: [GenomicAnalyzer.NucmerClusters.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.NucmerClusters.cs)
- Clustering port: [NucmerClusters.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/NucmerClusters.cs)
- Algorithm doc: [Common_Region_Detection.md](../../../algorithms/Sequence_Comparison/Common_Region_Detection.md)
- Oracle: [generate_nucmer_clusters_oracle.py](../../../../scripts/alignment/generate_nucmer_clusters_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference DNA sequence. |
| `query` | string | Yes | Query DNA sequence. |
| `minMatchLength` | integer | No | Minimum exact match length (nucmer -l, default 20). Default: `20`. |
| `mode` | string | No | Match mode: 'mum' (default), 'mumreference' or 'maxmatch'. Default: `"mum"`. |
| `bothStrands` | boolean | No | Also cluster the reverse-complement strand of the query (nucmer default true). Default: `true`. |
| `maxSeparation` | integer | No | Largest separation between matches of a cluster (nucmer -g, default 90). Default: `90`. |
| `minClusterLength` | integer | No | Smallest reported cluster length (nucmer -c, default 65). Default: `65`. |
| `diagonalDifference` | integer | No | Fixed diagonal difference that joins matches (nucmer -D, default 5). Default: `5`. |
| `diagonalFraction` | number | No | Fraction of the separation added to the allowed diagonal difference (nucmer -d, default 0.12). Default: `0.12`. |
| `useExtent` | boolean | No | Cluster length = extent in the reference instead of the sum of the match lengths (mgaps -e, default false). Default: `false`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `clusters` | array | Clusters `{ strand (forward|reverse), matches: [{ start1, start2, length, overlap }] }`, forward strand first, in the order mgaps prints them. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | mode must be 'mum', 'mumreference' or 'maxmatch' |

## Examples

### Example 1: Two segments split by a 4 base substitution block

**User Prompt:**
> Cluster the nucmer anchors of these two sequences with minimum match length 15.

**Expected Tool Call:**
```json
{"tool": "nucmer_clusters", "arguments": {"reference": "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA", "query": "AGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA", "minMatchLength": 15}}
```

**Response:**
```json
{"clusters": [{"strand": "forward", "matches": [{"start1": 6, "start2": 1, "length": 65, "overlap": 0}, {"start1": 75, "start2": 70, "length": 66, "overlap": 0}]}]}
```

mummer -mum -b -l 15 | mgaps prints "6 1 65 none - -" and "75 70 66 none 4 4": the two unique matches are one cluster of length 131.

## See Also

- [find_common_regions](find_common_regions.md) - the exact matches (MUMs / MEMs) themselves
- [find_longest_common_region](find_longest_common_region.md) - longest common region
- [local_align_top_n](local_align_top_n.md) - non-intersecting local alignments
