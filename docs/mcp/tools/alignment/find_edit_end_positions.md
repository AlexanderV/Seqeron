# find_edit_end_positions

Sellers k-differences search: end positions of approximate matches with their minimum edit distance.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `find_edit_end_positions` |
| **Method ID** | `ApproximateMatcher.FindEditEndPositions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Sellers (1980) approximate string matching (the k-differences problem; Navarro 2001 §5.1): reports every 0-based end position j of `sequence` at which some substring ending at j is within `maxEdits` Levenshtein edits of `pattern`, together with that minimum distance min_i ed(pattern, sequence[i..j]). Case-insensitive. Engine: Myers (1999) bit-parallel column scan (⌈m/64⌉ words, O(⌈m/64⌉·n)), identical to the Wagner–Fischer DP with a free text start; cross-checked against edlib HW mode. Output ordered by increasing end position.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L209](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L209)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to search in. |
| `pattern` | string | Yes | Pattern to find. |
| `maxEdits` | integer | Yes | Maximum allowed edit distance (>= 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | `{endPosition, distance}` — 0-based end position and the minimum edit distance of a substring ending there |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1005 | Pattern cannot be null or empty |
| 1006 | maxEdits must be >= 0 |

## Examples

### Example 1: Navarro survey/surgery

**User Prompt:**
> Where can 'survey' end in 'surgery' with at most 2 edits?

**Tool Call:**
```json
{
  "tool": "find_edit_end_positions",
  "arguments": {
    "sequence": "surgery",
    "pattern": "survey",
    "maxEdits": 2
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "endPosition": 4,
      "distance": 2
    },
    {
      "endPosition": 5,
      "distance": 2
    },
    {
      "endPosition": 6,
      "distance": 2
    }
  ]
}
```

### Example 2: Exact and 1-edit ends

**Tool Call:**
```json
{
  "tool": "find_edit_end_positions",
  "arguments": {
    "sequence": "GATTACAGATTTACA",
    "pattern": "TTAC",
    "maxEdits": 1
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "endPosition": 4,
      "distance": 1
    },
    {
      "endPosition": 5,
      "distance": 0
    },
    {
      "endPosition": 6,
      "distance": 1
    },
    {
      "endPosition": 12,
      "distance": 1
    },
    {
      "endPosition": 13,
      "distance": 0
    },
    {
      "endPosition": 14,
      "distance": 1
    }
  ]
}
```

## Worked Example

`TTAC` occurs exactly ending at 5 and 13 (distance 0); the neighbouring ends 4/6 and 12/14 are reachable with one deletion/insertion.

## See Also

- [find_with_edits](find_with_edits.md) — All approximate match windows with their alignments
- [edit_alignment](edit_alignment.md) — Optimal alignment of two strings

## References

- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L209](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L209)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
