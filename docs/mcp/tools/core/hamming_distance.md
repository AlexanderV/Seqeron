# hamming_distance

Calculate Hamming distance between two sequences of equal length.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `hamming_distance` |
| **Method ID** | `ApproximateMatcher.HammingDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the Hamming distance between two strings of equal length. The Hamming distance is the number of positions at which the corresponding characters differ. This metric is useful for comparing sequences of equal length, such as fixed-length DNA barcodes or aligned sequences.

## Core Documentation Reference

- Source: [ApproximateMatcher.cs#L163](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/ApproximateMatcher.cs#L163)
- Algorithm: Character-by-character comparison (case-insensitive)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | The first sequence (min length: 1) |
| `sequence2` | string | Yes | The second sequence (must be same length as sequence1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | integer | Number of positions with different characters (>= 0) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence1 cannot be null or empty |
| 1003 | Sequence2 cannot be null or empty |
| 1004 | Sequences must have equal length for Hamming distance |

## Examples

### Example 1: Identical sequences

**User Prompt:**
> What's the Hamming distance between ATGC and ATGC?

**Expected Tool Call:**
```json
{
  "tool": "hamming_distance",
  "arguments": {
    "sequence1": "ATGC",
    "sequence2": "ATGC"
  }
}
```

**Response:**
```json
{
  "distance": 0
}
```

### Example 2: One mismatch

**User Prompt:**
> Compare ATGC and ATGG

**Expected Tool Call:**
```json
{
  "tool": "hamming_distance",
  "arguments": {
    "sequence1": "ATGC",
    "sequence2": "ATGG"
  }
}
```

**Response:**
```json
{
  "distance": 1
}
```

## Performance

- **Time Complexity:** O(n) where n is the sequence length
- **Space Complexity:** O(1)
- **Note:** Comparison is case-insensitive

## See Also

- [edit_distance](edit_distance.md) - For sequences of different lengths
- [calculate_similarity](calculate_similarity.md) - K-mer based similarity
