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

Calculates the DUST score for a DNA sequence, which is used for low-complexity filtering in BLAST and other sequence analysis tools. The DUST algorithm identifies simple/repetitive regions by counting triplet word frequencies: score = Σ_t c_t·(c_t − 1)/2 / (ℓ − 1), where ℓ = L − 2 is the number of overlapping triplets (Morgulis et al. 2006; the normaliser used by NCBI dustmasker and lh3/sdust). Higher scores indicate lower complexity (more repetitive sequences); fewer than two triplets yield 0. DUST is defined for triplets only, so `wordSize` must be 3 (kept for compatibility; other values are rejected); the sourced k-mer generalisation is the C# API `SequenceComplexity.FindLongdustRegions` / `CalculateLongdustScore` (longdust, Li & Li 2025).

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L723](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L723)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to analyze (min length: 1) |
| `wordSize` | integer | No | Word size; must be 3 (default: 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `dustScore` | number | DUST score (higher values = lower complexity) |
| `wordSize` | integer | Word size used for calculation |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | The DUST score is defined for triplets only (word size 3) |

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
  "dustScore": 0.4,
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
  "dustScore": 5.0,
  "wordSize": 3
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(64) for triplet storage

## See Also

- [complexity_mask_low](complexity_mask_low.md) - Mask low-complexity regions
- [shannon_entropy](shannon_entropy.md) - Information-theoretic complexity
- [linguistic_complexity](linguistic_complexity.md) - K-mer diversity measure
