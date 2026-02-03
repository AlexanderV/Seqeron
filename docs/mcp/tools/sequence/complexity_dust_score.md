# complexity_dust_score

Calculate DUST score for low-complexity filtering.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complexity_dust_score` |
| **Method ID** | `SequenceComplexity.CalculateDustScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the DUST score for a DNA sequence, which is used for low-complexity filtering in BLAST and other sequence analysis tools. The DUST algorithm identifies simple/repetitive regions by counting triplet word frequencies. Higher scores indicate lower complexity (more repetitive sequences).

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L296](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceComplexity.cs#L296)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to analyze (min length: 1) |
| `wordSize` | integer | No | Word size for triplet counting (default: 3, minimum: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `dustScore` | number | DUST score (higher values = lower complexity) |
| `wordSize` | integer | Word size used for calculation |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Word size must be at least 1 |

## Examples

### Example 1: Complex DNA sequence

**User Prompt:**
> What's the DUST score for "ATGCGATCGATCG"?

**Expected Tool Call:**
```json
{
  "tool": "complexity_dust_score",
  "arguments": {
    "sequence": "ATGCGATCGATCG",
    "wordSize": 3
  }
}
```

**Response:**
```json
{
  "dustScore": 0.45,
  "wordSize": 3
}
```

### Example 2: Low complexity (repetitive) sequence

**User Prompt:**
> Calculate DUST score for "AAAAAAAAAAAA"

**Expected Tool Call:**
```json
{
  "tool": "complexity_dust_score",
  "arguments": {
    "sequence": "AAAAAAAAAAAA"
  }
}
```

**Response:**
```json
{
  "dustScore": 4.5,
  "wordSize": 3
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(4^wordSize) for triplet storage

## See Also

- [complexity_mask_low](complexity_mask_low.md) - Mask low-complexity regions
- [shannon_entropy](shannon_entropy.md) - Information-theoretic complexity
- [linguistic_complexity](linguistic_complexity.md) - K-mer diversity measure
