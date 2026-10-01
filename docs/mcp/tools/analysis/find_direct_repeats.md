# find_direct_repeats

Find exact direct repeats as maximal repeated pairs (MUMmer `repeat-match -f`); only A/C/G/T match; repeats longer than `maxLength` are dropped, not split.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_direct_repeats` |
| **Method ID** | `RepeatFinder.FindDirectRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds exact direct repeats reported as **maximal repeated pairs** (Gusfield 1997 §7.12; the
forward-strand output of MUMmer `repeat-match -f`): two copies at 0-based positions `i < j` that
cannot be extended to the left (`i = 0` or `S[i−1] ≠ S[j−1]`) or to the right (the length is the full
common-prefix length). Each repeat is therefore reported once at its full extent — its nested
sub-windows are not. Results are then filtered to `minLength ≤ length ≤ maxLength` (a maximal repeat
longer than `maxLength` is not reported, not truncated) and `spacing = secondPosition − firstPosition − length ≥ minSpacing`
(negative `minSpacing` admits overlapping copies). Matching is case-insensitive; only A/C/G/T match
(N/IUPAC never match — MUMmer `-n` convention). Output sorted by (firstPosition, secondPosition).
`minLength` must be ≥ 2 and `maxLength` ≥ `minLength`.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L3678](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L3678)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minLength` | integer | No | Minimum repeat length (default 5, ≥ 2) |
| `maxLength` | integer | No | Maximum repeat length (default 50) |
| `minSpacing` | integer | No | Minimum spacing between copies (default 1; 0 = abutting, negative = overlap allowed) |

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

### Example 2: Periodic copies — only the maximal pair

**Input:** `{ "sequence": "ACGTATTACGTATTACGTA", "minLength": 4, "maxLength": 50, "minSpacing": 1 }`
→ copies 0/7 and 7/14 are sub-windows of one overlapping maximal repeat (0, 7, length 12, spacing −5,
filtered by `minSpacing`); the only spaced maximal pair is "ACGTA" at 0 and 14 (repeat-match -f: `1 15 5`) →
`{ "items": [ { "firstPosition": 0, "secondPosition": 14, "repeatSequence": "ACGTA", "length": 5, "spacing": 9 } ] }`

### Example 3: No repeat

**Input:** `{ "sequence": "ACGTACGT", "minLength": 6 }`
→ **Response:** `{ "items": [] }`

## Performance

- **Time Complexity:** O(n log² n + z) (suffix array + LCP, bottom-up lcp-interval traversal; z = maximal pairs with length in range). **Space Complexity:** O(n).

## See Also

- [find_repeats](find_repeats.md)
- [find_tandem_repeats](find_tandem_repeats.md)
- [find_inverted_repeats](find_inverted_repeats.md)
