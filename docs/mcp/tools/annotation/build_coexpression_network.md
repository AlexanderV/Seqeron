# build_coexpression_network

Build a gene co-expression network from per-gene expression profiles.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `build_coexpression_network` |
| **Method ID** | `TranscriptomeAnalyzer.BuildCoExpressionNetwork` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For every unordered pair of genes (upper triangle, `i < j`), computes the Pearson correlation of
their expression profiles across samples and emits an edge `{ gene1, gene2, correlation }` when
`|correlation| >= correlationThreshold`. The `correlation` value is signed, so both strong positive
and strong negative relationships are retained. Genes with fewer than two samples correlate to `0`
and are therefore never linked.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L1117](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L1117)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `geneProfiles` | array | Yes | Per-gene profiles `{ geneId, expression[] }` |
| `correlationThreshold` | number | No | Minimum `|Pearson correlation|` for an edge (default `0.7`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `edges` | array | Co-expression edges `{ gene1, gene2, correlation }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | geneProfiles cannot be null |

## Examples

### Example 1: Perfectly (anti)correlated genes

`GENE1 = [1,2,3,4,5]`, `GENE2 = [2,4,6,8,10]` (corr `+1`), `GENE3 = [5,4,3,2,1]` (corr `-1` with
GENE1). All three pairs have magnitude `1`, so all edges are emitted at threshold `0.7`.

**Response:**
```json
{ "edges": [
  { "gene1": "GENE1", "gene2": "GENE2", "correlation": 1.0 },
  { "gene1": "GENE1", "gene2": "GENE3", "correlation": -1.0 },
  { "gene1": "GENE2", "gene2": "GENE3", "correlation": -1.0 }
] }
```

### Example 2: Single gene

**Response:**
```json
{ "edges": [] }
```

## Performance

- **Time Complexity:** O(g² · s) — gene pairs times samples per correlation
- **Space Complexity:** O(e) in the number of emitted edges

## See Also

- [pearson_correlation](pearson_correlation.md) - Pairwise correlation of two vectors
- [cluster_genes_by_expression](cluster_genes_by_expression.md) - Cluster genes by profile similarity
