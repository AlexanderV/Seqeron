# find_common_regions

Find common substrings between two DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_common_regions` |
| **Method ID** | `GenomicAnalyzer.FindCommonRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each start position in `sequence2`, finds the **single longest** contiguous substring of
length ≥ `minLength` that also occurs in `sequence1`, located via a suffix tree. Distinct such
substrings are reported once, with the first occurrence position in `sequence1` and the start
position in `sequence2`. Shorter prefixes sharing a start position with a longer match are not
reported. `minLength` below 1 is treated as 1. Both inputs must be valid DNA.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L306](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L306)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First DNA sequence (min length 1) |
| `sequence2` | string | Yes | Second DNA sequence (min length 1) |
| `minLength` | integer | Yes | Minimum common-region length |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Regions: `sequence, positionInFirst, positionInSecond, length` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: Identical sequences

**Input:** `{ "sequence1": "ACGT", "sequence2": "ACGT", "minLength": 2 }`

**Response:**
```json
{ "items": [
  { "sequence": "ACGT", "positionInFirst": 0, "positionInSecond": 0, "length": 4 },
  { "sequence": "CGT",  "positionInFirst": 1, "positionInSecond": 1, "length": 3 },
  { "sequence": "GT",   "positionInFirst": 2, "positionInSecond": 2, "length": 2 }
] }
```

### Example 2: No common region ≥ minLength

**Input:** `{ "sequence1": "AAAA", "sequence2": "CCCC", "minLength": 2 }`
→ **Response:** `{ "items": [] }`

## Performance

- **Time Complexity:** O(n + m·log m). **Space Complexity:** O(n) for the suffix tree.

## See Also

- [find_repeats](find_repeats.md)
- [generate_dot_plot](generate_dot_plot.md)
