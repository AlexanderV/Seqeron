# count_approximate_occurrences

Count approximate occurrences of a pattern in a sequence, allowing mismatches.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `count_approximate_occurrences` |
| **Method ID** | `ApproximateMatcher.CountApproximateOccurrences` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Counts the number of positions where a pattern approximately matches a sequence, allowing up to a specified number of mismatches (substitutions). This is useful for finding motifs with variations, SNP analysis, and fuzzy pattern matching in genomic data.

## Core Documentation Reference

- Source: [ApproximateMatcher.cs#L283](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/ApproximateMatcher.cs#L283)
- Algorithm: Sliding window with mismatch counting

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to search in (min length: 1) |
| `pattern` | string | Yes | The pattern to find (min length: 1) |
| `maxMismatches` | integer | Yes | Maximum number of allowed mismatches (>= 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `count` | integer | Number of approximate occurrences found (>= 0) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Pattern cannot be null or empty |
| 1005 | MaxMismatches cannot be negative |

## Examples

### Example 1: Exact matches only

**User Prompt:**
> How many exact "ATG" matches are in "ATGATGATG"?

**Expected Tool Call:**
```json
{
  "tool": "count_approximate_occurrences",
  "arguments": {
    "sequence": "ATGATGATG",
    "pattern": "ATG",
    "maxMismatches": 0
  }
}
```

**Response:**
```json
{
  "count": 3
}
```

### Example 2: With mismatches allowed

**User Prompt:**
> Find "ATG" with up to 1 mismatch in "ATGCTGATG"

**Expected Tool Call:**
```json
{
  "tool": "count_approximate_occurrences",
  "arguments": {
    "sequence": "ATGCTGATG",
    "pattern": "ATG",
    "maxMismatches": 1
  }
}
```

**Response:**
```json
{
  "count": 3
}
```

## Performance

- **Time Complexity:** O(n * m) where n is sequence length and m is pattern length
- **Space Complexity:** O(1)
- **Note:** Comparison is case-insensitive

## See Also

- [hamming_distance](hamming_distance.md) - Distance between equal-length sequences
- [edit_distance](edit_distance.md) - Full edit distance calculation
- [suffix_tree_count](suffix_tree_count.md) - Exact pattern counting
