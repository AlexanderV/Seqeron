# complexity_mask_low

Mask low-complexity regions using DUST algorithm.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complexity_mask_low` |
| **Method ID** | `SequenceComplexity.MaskLowComplexity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Masks low-complexity regions in a DNA sequence using the symmetric DUST (SDUST) algorithm of Morgulis et al. (2006) — output identical to lh3/sdust and NCBI dustmasker: every perfect interval of at most `windowSize` bases whose DUST score `Σ c(c−1)/2 / (ℓ−1)` exceeds `threshold` is masked. This is commonly used as a preprocessing step before BLAST searches or other sequence analyses to prevent spurious matches caused by simple/repetitive sequences. Low-complexity regions are replaced with a mask character (typically 'N').

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L346](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L346)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to mask (must be valid DNA) |
| `windowSize` | integer | No | Window size for analysis (default: 64; values 1–2 are rejected by the core, which requires ≥ 3) |
| `threshold` | number | No | DUST threshold above which to mask (default: 2.0) |
| `maskChar` | string | No | Character to use for masking (default: 'N') |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `maskedSequence` | string | Sequence with low-complexity regions masked |
| `originalLength` | integer | Original sequence length |
| `maskChar` | string | Character used for masking |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Window size must be at least 1 |
| 2001 | Invalid DNA sequence |

## Examples

### Example 1: Mask repetitive region

**User Prompt:**
> Mask the low-complexity regions in this sequence

**Expected Tool Call:**
```json
{
  "tool": "complexity_mask_low",
  "arguments": {
    "sequence": "ATGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGC",
    "windowSize": 64,
    "threshold": 2.0,
    "maskChar": "N"
  }
}
```

**Response:**
```json
{
  "maskedSequence": "ATGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATGC",
  "originalLength": 68,
  "maskChar": "N"
}
```

### Example 2: Custom threshold

**User Prompt:**
> Mask low-complexity regions with a stricter threshold of 1.0

**Expected Tool Call:**
```json
{
  "tool": "complexity_mask_low",
  "arguments": {
    "sequence": "ATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGC",
    "threshold": 1.0
  }
}
```

## Performance

- **Time Complexity:** O(n × w) typical; O(n × w³) worst case (SDUST perfect-interval search), w = window size
- **Space Complexity:** O(n) for the masked sequence

## See Also

- [complexity_dust_score](complexity_dust_score.md) - Calculate DUST score
- [shannon_entropy](shannon_entropy.md) - Information-theoretic complexity
