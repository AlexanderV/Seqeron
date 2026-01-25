# edit_distance

Calculate edit distance (Levenshtein distance) between two sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `edit_distance` |
| **Method ID** | `ApproximateMatcher.EditDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the edit distance (Levenshtein distance) between two strings. The edit distance is the minimum number of single-character edits (insertions, deletions, or substitutions) required to transform one string into the other. Unlike Hamming distance, sequences can have different lengths.

## Core Documentation Reference

- Source: [ApproximateMatcher.cs#L186](../../../../SuffixTree.Genomics/ApproximateMatcher.cs#L186)
- Algorithm: Wagner-Fischer dynamic programming (space-optimized)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | The first sequence (min length: 1) |
| `sequence2` | string | Yes | The second sequence (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | integer | Minimum number of edits needed (>= 0) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence1 cannot be null or empty |
| 1003 | Sequence2 cannot be null or empty |

## Examples

### Example 1: Identical sequences

**User Prompt:**
> What's the edit distance between ATGC and ATGC?

**Expected Tool Call:**
```json
{
  "tool": "edit_distance",
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

### Example 2: Classic example

**User Prompt:**
> How many edits to transform "kitten" to "sitting"?

**Expected Tool Call:**
```json
{
  "tool": "edit_distance",
  "arguments": {
    "sequence1": "kitten",
    "sequence2": "sitting"
  }
}
```

**Response:**
```json
{
  "distance": 3
}
```

## Performance

- **Time Complexity:** O(m * n) where m, n are sequence lengths
- **Space Complexity:** O(min(m, n)) using two-row optimization
- **Note:** Comparison is case-insensitive

## See Also

- [hamming_distance](hamming_distance.md) - For equal-length sequences only
- [count_approximate_occurrences](count_approximate_occurrences.md) - Find patterns with mismatches
