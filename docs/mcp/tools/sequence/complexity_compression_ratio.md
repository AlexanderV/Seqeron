# complexity_compression_ratio

Estimate sequence complexity using compression ratio.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complexity_compression_ratio` |
| **Method ID** | `SequenceComplexity.EstimateCompressionRatio` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Estimates sequence complexity using an LZ77-like compression approach. The algorithm counts unique substrings at various lengths and compares them to the expected number for a random sequence. Lower ratios indicate more repetitive/less complex sequences, while higher ratios indicate more complex sequences with more unique patterns.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L391](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceComplexity.cs#L391)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `compressionRatio` | number | Estimated compression ratio (0 to 1) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Complex DNA sequence

**User Prompt:**
> What's the compression ratio for "ATGCGATCGATCG"?

**Expected Tool Call:**
```json
{
  "tool": "complexity_compression_ratio",
  "arguments": {
    "sequence": "ATGCGATCGATCG"
  }
}
```

**Response:**
```json
{
  "compressionRatio": 0.85
}
```

### Example 2: Highly repetitive sequence

**User Prompt:**
> Calculate complexity for "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

**Expected Tool Call:**
```json
{
  "tool": "complexity_compression_ratio",
  "arguments": {
    "sequence": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  }
}
```

**Response:**
```json
{
  "compressionRatio": 0.15
}
```

## Performance

- **Time Complexity:** O(n × m) where n is sequence length and m is max substring length (10)
- **Space Complexity:** O(n × m) for unique substring storage

## See Also

- [complexity_dust_score](complexity_dust_score.md) - DUST algorithm complexity
- [shannon_entropy](shannon_entropy.md) - Information-theoretic complexity
- [linguistic_complexity](linguistic_complexity.md) - K-mer diversity measure
