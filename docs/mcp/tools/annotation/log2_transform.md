# log2_transform

Apply log2(x + pseudocount) to each value.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `log2_transform` |
| **Method ID** | `TranscriptomeAnalyzer.Log2Transform` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Applies `log2(value + pseudocount)` to each input value. The pseudocount (default 1) avoids `log(0)` and
stabilizes the variance of count-like expression data.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L237](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L237)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `values` | number[] | Yes | Values to transform |
| `pseudocount` | number | No | Pseudocount added before log2 (default 1.0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `transformed` | number[] | `log2(value + pseudocount)` per input |

## Errors

| Code | Message |
|------|---------|
| 1001 | values cannot be null |

## Examples

### Example 1: Default pseudocount

`[0, 1, 3, 7]` → `log2([1, 2, 4, 8])`:

**Response:**
```json
{ "transformed": [0.0, 1.0, 2.0, 3.0] }
```

### Example 2: Custom pseudocount

`[1, 5]` with pseudocount 3 → `log2([4, 8])`:

**Response:**
```json
{ "transformed": [2.0, 3.0] }
```

## Performance

- **Time Complexity:** O(n)
- **Space Complexity:** O(n)

## See Also

- [quantile_normalize](quantile_normalize.md) — cross-sample normalization
- [calculate_tpm](calculate_tpm.md) — TPM abundance
