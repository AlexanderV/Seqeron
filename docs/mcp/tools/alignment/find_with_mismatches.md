# find_with_mismatches

Find all occurrences of pattern inside sequence allowing up to maxMismatches substitutions (Hamming-style, fixed-length window). Reports position, matched substring, mismatch count, and 0-based mismatch positions.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `find_with_mismatches` |
| **Method ID** | `ApproximateMatcher.FindWithMismatches` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Find all occurrences of pattern inside sequence allowing up to maxMismatches substitutions (Hamming-style, fixed-length window). Reports position, matched substring, mismatch count, and 0-based mismatch positions.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L26](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L26)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to search in. |
| `pattern` | string | Yes | Pattern to find. |
| `maxMismatches` | integer | Yes | Maximum number of allowed mismatches (>= 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | Matches, each with position, matchedSequence, distance, mismatchPositions, mismatchType. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1005 | Pattern cannot be null or empty |
| 1006 | maxMismatches must be >= 0 |

## Examples

### Example 1: Exact occurrences

**Tool Call:**
```json
{
  "tool": "find_with_mismatches",
  "arguments": {
    "sequence": "ACGTACGT",
    "pattern": "ACG",
    "maxMismatches": 0
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "position": 0,
      "matchedSequence": "ACG",
      "distance": 0,
      "mismatchPositions": [],
      "mismatchType": "Substitution"
    },
    {
      "position": 4,
      "matchedSequence": "ACG",
      "distance": 0,
      "mismatchPositions": [],
      "mismatchType": "Substitution"
    }
  ]
}
```

### Example 2: Allow 1 mismatch

**Tool Call:**
```json
{
  "tool": "find_with_mismatches",
  "arguments": {
    "sequence": "AAAA",
    "pattern": "AA",
    "maxMismatches": 1
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "position": 0,
      "matchedSequence": "AA",
      "distance": 0,
      "mismatchPositions": [],
      "mismatchType": "Substitution"
    },
    {
      "position": 1,
      "matchedSequence": "AA",
      "distance": 0,
      "mismatchPositions": [],
      "mismatchType": "Substitution"
    },
    {
      "position": 2,
      "matchedSequence": "AA",
      "distance": 0,
      "mismatchPositions": [],
      "mismatchType": "Substitution"
    }
  ]
}
```

## Worked Example

`ACG` occurs exactly at positions 0 and 4 of `ACGTACGT`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L26](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L26)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
