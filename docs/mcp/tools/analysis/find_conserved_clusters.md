# find_conserved_clusters

Find conserved gene clusters (common intervals) across genomes.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_conserved_clusters` |
| **Method ID** | `ComparativeGenomics.FindConservedClusters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds conserved gene clusters as **common intervals** of the ortholog-group permutations: a set
of ortholog-group labels that occupies a contiguous window in *every* genome (Uno & Yagiura 2000;
Heber & Stoye 2001). Genes are mapped to groups via `orthologGroups`; genes with no group break
windows. Requires at least two genomes (fewer makes the question vacuous → empty). Results are
sorted by size, then lexicographically. `maxGap` is retained for compatibility but the model is
strict (gap-free).

## Core Documentation Reference

- Source: [ComparativeGenomics.cs#L915](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs#L915)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genomes` | object[][] | Yes | One gene array per genome |
| `orthologGroups` | object | Yes | Map gene id → ortholog-group id |
| `minClusterSize` | integer | No | Minimum distinct groups per cluster (default 3, ≥ 2) |
| `maxGap` | integer | No | Retained for compatibility (default 2) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `clusters` | string[][] | Each cluster is a sorted list of ortholog-group ids |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one genome is required |
| 1001 | Ortholog groups map cannot be null |

## Examples

### Example 1: Two identical genomes

Genes A,B,C (groups g1,g2,g3) in the same order in both genomes, `minClusterSize` 2:

**Response:**
```json
{ "clusters": [ ["g1","g2"], ["g2","g3"], ["g1","g2","g3"] ] }
```

### Example 2: A single genome

One genome → **Response:** `{ "clusters": [] }` (a common interval requires ≥ 2 genomes).

## Performance

- **Time Complexity:** O(g · n²) over g genomes of n genes. **Space Complexity:** O(n²) candidate sets.

## See Also

- [find_syntenic_blocks](find_syntenic_blocks.md)
- [find_orthologs](find_orthologs.md)
