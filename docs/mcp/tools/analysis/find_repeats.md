# find_repeats

All repeated substrings of length ≥ `minLength` in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_repeats` |
| **Method ID** | `GenomicAnalyzer.FindRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **every distinct substring** of length ≥ `minLength` that occurs at least twice
in a DNA sequence, reporting each with all its 0-based occurrence positions. Uses the
sequence's suffix tree / longest-common-prefix of adjacent sorted suffixes (CMU
15-451 Lecture #10 §2.1). Overlapping occurrences are counted (e.g. `AA` in `AAAA`
occurs at 0, 1, 2).

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L49](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L49)
- Algorithm: `docs/algorithms/Repeat_Analysis/Repeat_Detection.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minLength` | integer | Yes | Minimum repeat length (≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Repeats: `{ sequence, positions[], length, count }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | Minimum repeat length must be at least 1 |

## Examples

### Example 1: Homopolymer repeats (AAAA, minLength 2)

**User Prompt:**
> Find repeats of length ≥ 2 in "AAAA".

**Expected Tool Call:**
```json
{
  "tool": "find_repeats",
  "arguments": { "sequence": "AAAA", "minLength": 2 }
}
```

**Response:**
```json
{ "items": [ { "sequence": "AA", "positions": [0, 1, 2], "length": 2, "count": 3 }, { "sequence": "AAA", "positions": [0, 1], "length": 3, "count": 2 } ] }
```

### Example 2: Tandem ATG (ATGATG, minLength 3)

**User Prompt:**
> Repeats of length ≥ 3 in "ATGATG".

**Expected Tool Call:**
```json
{
  "tool": "find_repeats",
  "arguments": { "sequence": "ATGATG", "minLength": 3 }
}
```

**Response:**
```json
{ "items": [ { "sequence": "ATG", "positions": [0, 3], "length": 3, "count": 2 } ] }
```
Only ATG recurs; TGA and GAT occur once each.

## Performance

- **Time Complexity:** O(n²) worst case (adjacent-suffix prefix enumeration).
- **Space Complexity:** O(n) for the suffix tree.

## See Also

- [find_tandem_repeats](find_tandem_repeats.md) — consecutive repeat units
- [find_direct_repeats](find_direct_repeats.md) — spaced identical copies
- [find_motif](find_motif.md) — locate a known motif
