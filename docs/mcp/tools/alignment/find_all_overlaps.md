# find_all_overlaps

Compute all suffix-of-i / prefix-of-j overlaps between ordered read pairs above minOverlap length and minIdentity ratio; each result is an overlap-graph edge i -> j.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `find_all_overlaps` |
| **Method ID** | `SequenceAssembler.FindAllOverlaps` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Compute all suffix-of-i / prefix-of-j overlaps between ordered read pairs above minOverlap length and minIdentity ratio; each result is an overlap-graph edge i -> j.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L142](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L142)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reads` | array<string> | Yes | Input sequence reads. |
| `minOverlap` | integer | No | Minimum overlap length (bp). (default: 20) |
| `minIdentity` | number | No | Minimum identity ratio (0.0-1.0). (default: 0.9) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | List of overlaps, each with readIndex1, readIndex2, overlapLength, position1, position2. |

## Examples

### Example 1: One directed overlap

**Tool Call:**
```json
{
  "tool": "find_all_overlaps",
  "arguments": {
    "reads": [
      "AAAAACCCCC",
      "CCCCCGGGGG"
    ],
    "minOverlap": 5,
    "minIdentity": 0.9
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "readIndex1": 0,
      "readIndex2": 1,
      "overlapLength": 5,
      "position1": 5,
      "position2": 0
    }
  ]
}
```

### Example 2: No overlaps

**Tool Call:**
```json
{
  "tool": "find_all_overlaps",
  "arguments": {
    "reads": [
      "AAAAA",
      "TTTTT"
    ],
    "minOverlap": 5,
    "minIdentity": 0.9
  }
}
```

**Response:**
```json
{
  "items": []
}
```

## Worked Example

The 5-bp suffix `CCCCC` of read 0 equals the prefix of read 1, yielding one edge 0 -> 1 of length 5 starting at position 5 in read 0.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L142](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L142)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
