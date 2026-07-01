# over_representation_analysis

Pathway / gene-set over-representation analysis (ORA).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `over_representation_analysis` |
| **Method ID** | `TranscriptomeAnalyzer.PerformOverRepresentationAnalysis` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Tests each pathway for over-representation among the differentially expressed (DE) genes using a
hypergeometric / Fisher's-exact approximation over a background of `backgroundGeneCount` genes. Pathways with
no DE overlap are omitted. Each reported pathway records its gene count, the number of overlapping DE genes,
an enrichment score (> 1 = enriched), a p-value, and the overlapping gene IDs.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L638](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L638)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `differentiallyExpressedGenes` | string[] | Yes | DE gene IDs |
| `pathways` | array | Yes | Pathway definitions `{ pathwayId, pathwayName, genes[] }` |
| `backgroundGeneCount` | integer | Yes | Universe size (> 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `results[].pathwayId` | string | Pathway id |
| `results[].pathwayName` | string | Pathway name |
| `results[].genesInPathway` | integer | Pathway size |
| `results[].overlappingGenes` | integer | Overlap with DE set |
| `results[].enrichmentScore` | number | > 1 = enriched |
| `results[].pValue` | number | Over-representation p-value |
| `results[].genes` | string[] | Overlapping gene IDs |

## Errors

| Code | Message |
|------|---------|
| 1001 | differentiallyExpressedGenes cannot be null |
| 1001 | pathways cannot be null |
| 1001 | backgroundGeneCount must be positive. |

## Examples

### Example 1: Enriched pathway

DE `{A,B,C,D,E}` against pathway `{A,B,C,X,Y}` (background 1000) → 3 overlapping genes, enriched.

### Example 2: No overlap

A pathway sharing no DE gene is omitted:

**Response:**
```json
{ "results": [] }
```

## Performance

- **Time Complexity:** O(p · g) for p pathways and gene-set sizes g
- **Space Complexity:** O(p)

## See Also

- [differential_expression](differential_expression.md) — produce the DE gene set
- [enrichment_score](enrichment_score.md) — ranked-list GSEA enrichment
