# find_best_match

Return the single best (minimum Hamming distance) fixed-length window of pattern inside sequence; ties resolve to the leftmost window and an exact (distance 0) match short-circuits the scan.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `find_best_match` |
| **Method ID** | `ApproximateMatcher.FindBestMatch` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Return the single best (minimum Hamming distance) fixed-length window of pattern inside sequence; ties resolve to the leftmost window and an exact (distance 0) match short-circuits the scan.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L248](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L248)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to search in. |
| `pattern` | string | Yes | Pattern to find. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `match` | object|null | Best match (position, matchedSequence, distance, mismatchPositions, mismatchType), or null if pattern is longer than sequence. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1005 | Pattern cannot be null or empty |

## Examples

### Example 1: Exact match

**Tool Call:**
```json
{
  "tool": "find_best_match",
  "arguments": {
    "sequence": "ACGTACGT",
    "pattern": "ACG"
  }
}
```

**Response:**
```json
{
  "match": {
    "position": 0,
    "matchedSequence": "ACG",
    "distance": 0,
    "mismatchPositions": [],
    "mismatchType": "Substitution"
  }
}
```

### Example 2: Best imperfect

**Tool Call:**
```json
{
  "tool": "find_best_match",
  "arguments": {
    "sequence": "AAAA",
    "pattern": "TT"
  }
}
```

**Response:**
```json
{
  "match": {
    "position": 0,
    "matchedSequence": "AA",
    "distance": 2,
    "mismatchPositions": [
      0,
      1
    ],
    "mismatchType": "Substitution"
  }
}
```

## Worked Example

`ACG` occurs exactly at position 0, so the scan short-circuits with distance 0.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L248](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L248)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
