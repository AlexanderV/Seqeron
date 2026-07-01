# pearson_correlation

Pearson product-moment correlation of two expression vectors.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `pearson_correlation` |
| **Method ID** | `TranscriptomeAnalyzer.CalculatePearsonCorrelation` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the Pearson correlation coefficient `r ∈ [-1, 1]` between two equal-length expression vectors. `r`
is 1 for a perfect positive linear relationship, −1 for a perfect negative one, and near 0 when the vectors
are uncorrelated.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L1084](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L1084)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `expression1` | number[] | Yes | First expression vector |
| `expression2` | number[] | Yes | Second vector (same length) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `correlation` | number | Pearson r in [-1, 1] |

## Errors

| Code | Message |
|------|---------|
| 1001 | expression1 cannot be null |
| 1001 | expression2 cannot be null |
| 1001 | Expression vectors cannot be empty. |
| 1001 | Expression vectors must have equal length |

## Examples

### Example 1: Perfect positive

`[1,2,3,4,5]` vs `[2,4,6,8,10]`:

**Response:**
```json
{ "correlation": 1.0 }
```

### Example 2: Perfect negative

`[1,2,3,4,5]` vs `[5,4,3,2,1]`:

**Response:**
```json
{ "correlation": -1.0 }
```

## Performance

- **Time Complexity:** O(n)
- **Space Complexity:** O(1)

## See Also

- [build_coexpression_network](build_coexpression_network.md) — correlation-thresholded gene network
- [cluster_genes_by_expression](cluster_genes_by_expression.md) — correlation-based clustering
