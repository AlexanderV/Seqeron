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

Optional `insertionCost`, `deletionCost`, `substitutionCost` (default 1 each; unit costs keep the original behaviour) return the weighted Levenshtein distance with the semantics of rapidfuzz `Levenshtein.distance(s1, s2, weights=(insertion, deletion, substitution))`: transforming `sequence1` into `sequence2`, an insertion adds a character of `sequence2` and a deletion removes a character of `sequence1` (so swapping the sequences swaps those two costs). Delegates to `ApproximateMatcher.EditDistance(s1, s2, insertionCost, deletionCost, substitutionCost)`.

## Core Documentation Reference

- Source: [ApproximateMatcher.cs#L186](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L186)
- Algorithm: Myers (1999) bit-parallel engine (unit costs); weighted Wagner–Fischer dynamic programming (non-unit costs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | The first sequence (min length: 1) |
| `sequence2` | string | Yes | The second sequence (min length: 1) |
| `insertionCost` | integer | No | Cost of inserting a character of `sequence2` (>= 0; default 1) |
| `deletionCost` | integer | No | Cost of deleting a character of `sequence1` (>= 0; default 1) |
| `substitutionCost` | integer | No | Cost of a substitution (>= 0; default 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | integer | Minimum number of edits needed (>= 0) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence1 cannot be null or empty |
| 1003 | Sequence2 cannot be null or empty |
| 1007 | insertionCost must be >= 0 |
| 1008 | deletionCost must be >= 0 |
| 1009 | substitutionCost must be >= 0 |

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

### Example 3: Weighted costs

**User Prompt:**
> Edit distance from kitten to sitting if a substitution costs 2?

**Expected Tool Call:**
```json
{
  "tool": "edit_distance",
  "arguments": {
    "sequence1": "kitten",
    "sequence2": "sitting",
    "insertionCost": 1,
    "deletionCost": 1,
    "substitutionCost": 2
  }
}
```

**Response:**
```json
{
  "distance": 5
}
```

(rapidfuzz 3.14.6 `Levenshtein.distance("kitten", "sitting", weights=(1, 1, 2))` = 5.)

## Performance

- **Time Complexity:** O(⌈min(m, n)/64⌉ · max(m, n)) for unit costs (Myers bit-parallel); O(m · n) for weighted costs
- **Space Complexity:** O(min(m, n))
- **Note:** Comparison is case-sensitive (ordinal)

## See Also

- [hamming_distance](hamming_distance.md) - For equal-length sequences only
- [count_approximate_occurrences](count_approximate_occurrences.md) - Find patterns with mismatches
