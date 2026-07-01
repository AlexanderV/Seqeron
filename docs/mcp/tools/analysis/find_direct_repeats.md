# find_direct_repeats

Find direct repeats (identical copies separated by a spacer).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_direct_repeats` |
| **Method ID** | `RepeatFinder.FindDirectRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds direct repeats: identical substrings of length in `[minLength, maxLength]` that appear
twice with at least `minSpacing` bases between the two copies. Each result reports the two
0-based start positions, the repeat sequence, its length, and the spacing
(`secondPosition − firstPosition − length`). Matching is case-insensitive via a suffix tree.
`minLength` must be ≥ 2 and `maxLength` ≥ `minLength`.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L796](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L796)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minLength` | integer | No | Minimum repeat length (default 5, ≥ 2) |
| `maxLength` | integer | No | Maximum repeat length (default 50) |
| `minSpacing` | integer | No | Minimum spacing between copies (default 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Repeats: `firstPosition, secondPosition, repeatSequence, length, spacing` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | minLength must be ≥ 2 |
| 1003 | maxLength must be ≥ minLength |

## Examples

### Example 1: One spaced repeat

**Input:** `{ "sequence": "ATGCGGATGC", "minLength": 4, "maxLength": 4, "minSpacing": 1 }`
→ "ATGC" at 0 and 6, spacing 2 →
`{ "items": [ { "firstPosition": 0, "secondPosition": 6, "repeatSequence": "ATGC", "length": 4, "spacing": 2 } ] }`

### Example 2: No repeat

**Input:** `{ "sequence": "ACGTACGT", "minLength": 6 }`
→ **Response:** `{ "items": [] }`

## Performance

- **Time Complexity:** O(n · (maxLength − minLength)). **Space Complexity:** O(n) suffix tree.

## See Also

- [find_repeats](find_repeats.md)
- [find_tandem_repeats](find_tandem_repeats.md)
- [find_inverted_repeats](find_inverted_repeats.md)
