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

Returns the normalized Lempel–Ziv (1976) complexity c / (n / log_b n) (Zhang et al. 2009), where c is the number of components in the LZ76 exhaustive history (Kaspar–Schuster scan), n the length and b the number of distinct symbols (clamped to ≥ 2). Values near 1 indicate random-like sequences; lower values indicate more repetitive/compressible sequences. Finite sequences can exceed 1. See [Lempel_Ziv_Complexity.md](../../../algorithms/Complexity/Lempel_Ziv_Complexity.md).

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L1631](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L1631)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `compressionRatio` | number | Normalized Lempel–Ziv complexity (≥ 0; ≈ 1 for random sequences, may exceed 1) |

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
  "compressionRatio": 0.9962722318072171
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
  "compressionRatio": 0.26609640474436813
}
```

## Performance

- **Time Complexity:** O(n²) worst case (Kaspar–Schuster scan of the LZ76 exhaustive history)
- **Space Complexity:** O(σ) (distinct-symbol count only)

## See Also

- [complexity_dust_score](complexity_dust_score.md) - DUST algorithm complexity
- [shannon_entropy](shannon_entropy.md) - Information-theoretic complexity
- [linguistic_complexity](linguistic_complexity.md) - K-mer diversity measure
