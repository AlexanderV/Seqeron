# enrichment_score

Compute a GSEA-like running-sum enrichment score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `enrichment_score` |
| **Method ID** | `TranscriptomeAnalyzer.CalculateEnrichmentScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Walks the ranked gene list, adding `1/hits` when a gene is in the gene set and subtracting `1/misses`
otherwise, and returns the maximum-magnitude running-sum deviation. A positive score means the gene set is
enriched near the top of the ranking; a negative score means it is enriched near the bottom. The score is 0
when the ranked list contains all or none of the set (no meaningful walk).

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L697](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L697)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rankedGenes` | string[] | Yes | Gene IDs ordered by ranking metric |
| `geneSet` | string[] | Yes | Gene IDs in the tested set |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | Max running-sum deviation |

## Errors

| Code | Message |
|------|---------|
| 1001 | rankedGenes cannot be null |
| 1001 | geneSet cannot be null |
| 1001 | rankedGenes cannot be empty. |
| 1001 | geneSet cannot be empty. |

## Examples

### Example 1: Set enriched at the top

`{A,B}` at the top of `[A,B,C,D]`:

**Response:**
```json
{ "score": 1.0 }
```

### Example 2: Set enriched at the bottom

`{C,D}`:

**Response:**
```json
{ "score": -1.0 }
```

## Performance

- **Time Complexity:** O(n) for n ranked genes
- **Space Complexity:** O(1)

## See Also

- [over_representation_analysis](over_representation_analysis.md) — hypergeometric ORA
- [differential_expression](differential_expression.md) — produce the ranking metric
