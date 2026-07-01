# quantile_normalize

Quantile-normalize multiple expression vectors.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `quantile_normalize` |
| **Method ID** | `TranscriptomeAnalyzer.QuantileNormalize` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Applies Bolstad (2003) quantile normalization: for each rank the mean of the sorted values across samples
is computed, then each value is replaced by the rank mean corresponding to its within-sample rank (ties get
the average of the spanning rank means). After normalization every sample shares an identical value
distribution. All input vectors must be the same length.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L175](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L175)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `samples` | number[][] | Yes | Per-sample expression vectors (all equal length) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `normalized` | number[][] | Quantile-normalized vectors (same shape as input) |

## Errors

| Code | Message |
|------|---------|
| 1001 | samples cannot be null |
| 1001 | samples cannot be empty. |
| 1001 | All samples must have the same length |

## Examples

### Example 1: Two-sample normalization

Rank means of `[2,3,4]` and `[8,6,5]` are `[3.5, 4.5, 6.0]`:

**Response:**
```json
{ "normalized": [ [3.5, 4.5, 6.0], [6.0, 4.5, 3.5] ] }
```

## Performance

- **Time Complexity:** O(k · n log n) for k samples of length n
- **Space Complexity:** O(k · n)

## See Also

- [log2_transform](log2_transform.md) — variance stabilization
- [differential_expression](differential_expression.md) — group comparison
