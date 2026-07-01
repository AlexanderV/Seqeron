# cluster_genes_by_expression

Cluster genes by expression-profile correlation.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `cluster_genes_by_expression` |
| **Method ID** | `TranscriptomeAnalyzer.ClusterGenesByExpression` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Clusters genes by similarity of their expression profiles with a k-means-like procedure: the first
`numClusters` genes seed the clusters, and each remaining gene is assigned to the cluster with the highest
mean Pearson correlation to its members. Each cluster reports its member genes, mean within-cluster
correlation, a representative gene, and (placeholder) enriched functions.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L1140](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L1140)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `geneProfiles` | array | Yes | Per-gene `{ geneId, expression[] }` |
| `numClusters` | integer | No | Number of clusters to form (default 5) |
| `correlationThreshold` | number | No | Informational min within-cluster correlation (default 0.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `clusters[].clusterId` | integer | Cluster id |
| `clusters[].genes` | string[] | Member gene ids |
| `clusters[].meanCorrelation` | number | Mean within-cluster correlation |
| `clusters[].representativeGene` | string | Representative gene |
| `clusters[].enrichedFunctions` | string[] | Enriched functions (placeholder) |

## Errors

| Code | Message |
|------|---------|
| 1001 | geneProfiles cannot be null |

## Examples

### Example 1: Two correlation clusters

G1/G3 (positively correlated) cluster together; G2/G4 form the other cluster:

**Response:**
```json
{ "clusters": [ { "genes": ["G1", "G3"] }, { "genes": ["G2", "G4"] } ] }
```

## Performance

- **Time Complexity:** O(g² · s) for g genes, s samples
- **Space Complexity:** O(g)

## See Also

- [pearson_correlation](pearson_correlation.md) — pairwise correlation
- [build_coexpression_network](build_coexpression_network.md) — correlation-thresholded network
