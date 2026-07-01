# find_with_edits

Find all approximate matches of pattern in sequence allowing up to maxEdits Levenshtein edits (insertions, deletions, substitutions), scanning variable-length windows of length pattern.Length +/- maxEdits.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `find_with_edits` |
| **Method ID** | `ApproximateMatcher.FindWithEdits` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Find all approximate matches of pattern in sequence allowing up to maxEdits Levenshtein edits (insertions, deletions, substitutions), scanning variable-length windows of length pattern.Length +/- maxEdits.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L118](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L118)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to search in. |
| `pattern` | string | Yes | Pattern to find. |
| `maxEdits` | integer | Yes | Maximum allowed edit distance (>= 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | Matches, each with position, matchedSequence, distance, mismatchPositions, mismatchType. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1005 | Pattern cannot be null or empty |
| 1006 | maxEdits must be >= 0 |

## Examples

### Example 1: Edit distance 1

**Tool Call:**
```json
{
  "tool": "find_with_edits",
  "arguments": {
    "sequence": "ACGT",
    "pattern": "AGT",
    "maxEdits": 1
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "position": 0,
      "matchedSequence": "ACGT",
      "distance": 1,
      "mismatchPositions": [],
      "mismatchType": "Edit"
    },
    {
      "position": 1,
      "matchedSequence": "CGT",
      "distance": 1,
      "mismatchPositions": [],
      "mismatchType": "Substitution"
    },
    {
      "position": 2,
      "matchedSequence": "GT",
      "distance": 1,
      "mismatchPositions": [],
      "mismatchType": "Edit"
    }
  ]
}
```

### Example 2: Exact only

**Tool Call:**
```json
{
  "tool": "find_with_edits",
  "arguments": {
    "sequence": "ACGTACGT",
    "pattern": "ACGT",
    "maxEdits": 0
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "position": 0,
      "matchedSequence": "ACGT",
      "distance": 0,
      "mismatchPositions": [],
      "mismatchType": "Substitution"
    },
    {
      "position": 4,
      "matchedSequence": "ACGT",
      "distance": 0,
      "mismatchPositions": [],
      "mismatchType": "Substitution"
    }
  ]
}
```

## Worked Example

Windows of length 2..4 are scored against `AGT`; three windows have edit distance 1 to the pattern.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L118](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L118)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
