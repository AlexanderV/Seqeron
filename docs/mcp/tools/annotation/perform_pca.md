# perform_pca

Project samples onto the first two principal components (approximate PCA).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `perform_pca` |
| **Method ID** | `TranscriptomeAnalyzer.PerformPCA` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Selects the `topGenes` most variable genes across the samples, then approximates each sample's first two
principal-component scores as the sum of the first half (`PC1`) and second half (`PC2`) of its selected gene
values. With fewer than two samples every sample projects to the origin (0, 0). This is a lightweight
approximation, not a full SVD-based PCA.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L1243](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L1243)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `samples` | array | Yes | Per-sample `{ sampleId, expression[] }` (equal length) |
| `topGenes` | integer | No | Number of top-variable genes to use (default 500) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `points[].sampleId` | string | Sample id |
| `points[].pc1` | number | PC1 score |
| `points[].pc2` | number | PC2 score |

## Errors

| Code | Message |
|------|---------|
| 1001 | samples cannot be null |
| 1001 | samples cannot be empty. |
| 1001 | All sample expression vectors must have the same length |

## Examples

### Example 1: Two contrasting samples

All four genes are selected; A → PC1 = 4+4 = 8, PC2 = 8; B → 0, 0:

**Response:**
```json
{ "points": [ { "sampleId": "A", "pc1": 8, "pc2": 8 }, { "sampleId": "B", "pc1": 0, "pc2": 0 } ] }
```

### Example 2: Single sample

A single sample projects to the origin:

**Response:**
```json
{ "points": [ { "sampleId": "A", "pc1": 0, "pc2": 0 } ] }
```

## Performance

- **Time Complexity:** O(s · g) for s samples and g genes
- **Space Complexity:** O(g)

## See Also

- [cluster_genes_by_expression](cluster_genes_by_expression.md) — correlation-based gene clustering
- [rnaseq_quality_metrics](rnaseq_quality_metrics.md) — library QC metrics
