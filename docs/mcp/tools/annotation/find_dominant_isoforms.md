# find_dominant_isoforms

Identify each gene's dominant transcript isoform.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_dominant_isoforms` |
| **Method ID** | `TranscriptomeAnalyzer.FindDominantIsoforms` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Groups isoforms by gene and reports, for each gene, the isoform with the highest expression together with
its **dominance ratio** — its expression divided by the total expression of all the gene's isoforms.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L1005](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L1005)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `isoforms` | array | Yes | Transcript isoforms with expression values |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `dominants[].geneId` | string | Gene id |
| `dominants[].dominantIsoform` | object | Highest-expression isoform |
| `dominants[].dominanceRatio` | number | dominant expression / total gene expression |

## Errors

| Code | Message |
|------|---------|
| 1001 | isoforms cannot be null |

## Examples

### Example 1: Dominant isoform of a gene

GENE1 with isoform expressions 100/50/25 → dominant TX1, dominance ratio `100 / 175 ≈ 0.571`:

**Response:**
```json
{ "dominants": [ { "geneId": "GENE1", "dominantIsoform": { "transcriptId": "TX1" }, "dominanceRatio": 0.5714285714285714 } ] }
```

## Performance

- **Time Complexity:** O(n) for n isoforms
- **Space Complexity:** O(g) for g genes

## See Also

- [detect_isoform_switching](detect_isoform_switching.md) — dominant-isoform switching across conditions
