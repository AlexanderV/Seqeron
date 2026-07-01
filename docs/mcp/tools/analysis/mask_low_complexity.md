# mask_low_complexity

Mask low-complexity windows of a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mask_low_complexity` |
| **Method ID** | `SequenceComplexity.MaskLowComplexity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Masks **low-complexity windows** of a DNA sequence with a chosen character. Each
sliding window whose DUST score exceeds `threshold` has all of its positions replaced
with `maskChar`. The result is the same length as the input. Sequences shorter than
`windowSize` are returned unchanged.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L411](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L411)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `windowSize` | integer | No | Window size (default 64) |
| `threshold` | number | No | DUST threshold above which to mask (default 2.0) |
| `maskChar` | string | No | Mask character (default `N`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `masked` | string | Masked sequence (same length as input) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: Fully masked poly-A

**User Prompt:**
> Mask low-complexity regions of a 100-nt poly-A tract with 'X' (window 64, threshold 1.0).

**Expected Tool Call:**
```json
{
  "tool": "mask_low_complexity",
  "arguments": { "sequence": "AAAA…(100×)", "windowSize": 64, "threshold": 1.0, "maskChar": "X" }
}
```

**Response:**
```json
{ "masked": "XXXX…(100×)" }
```
The poly-A DUST score (31.0) exceeds the threshold, so every position is masked.

### Example 2: High complexity preserved

**User Prompt:**
> Mask a varied 78-bp sequence at threshold 10.0.

**Expected Tool Call:**
```json
{
  "tool": "mask_low_complexity",
  "arguments": { "sequence": "ATGCTAGCATGCA…(78 bp)", "windowSize": 64, "threshold": 10.0 }
}
```

**Response:**
```json
{ "masked": "ATGCTAGCATGCA…(78 bp, unchanged)" }
```

## Performance

- **Time Complexity:** O(n · windowSize).
- **Space Complexity:** O(n).

## See Also

- [find_low_complexity_regions](find_low_complexity_regions.md)
- [dust_score](dust_score.md)
