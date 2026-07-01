# find_overlap

Detect the longest suffix-of-sequence1 / prefix-of-sequence2 overlap satisfying minOverlap and minIdentity thresholds. Returns null if no qualifying overlap exists.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `find_overlap` |
| **Method ID** | `SequenceAssembler.FindOverlap` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Detect the longest suffix-of-sequence1 / prefix-of-sequence2 overlap satisfying minOverlap and minIdentity thresholds. Returns null if no qualifying overlap exists.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L221](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L221)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (suffix candidate). |
| `sequence2` | string | Yes | Second sequence (prefix candidate). |
| `minOverlap` | integer | No | Minimum overlap length (bp). (default: 20) |
| `minIdentity` | number | No | Minimum identity ratio (0.0-1.0). (default: 0.9) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `overlap` | object|null | Overlap descriptor (length, position1, position2) or null if none qualifies. |

## Examples

### Example 1: Overlap found

**Tool Call:**
```json
{
  "tool": "find_overlap",
  "arguments": {
    "sequence1": "AAAAACCCCC",
    "sequence2": "CCCCCGGGGG",
    "minOverlap": 5,
    "minIdentity": 0.9
  }
}
```

**Response:**
```json
{
  "overlap": {
    "length": 5,
    "position1": 5,
    "position2": 0
  }
}
```

### Example 2: No overlap

**Tool Call:**
```json
{
  "tool": "find_overlap",
  "arguments": {
    "sequence1": "AAAAA",
    "sequence2": "TTTTT",
    "minOverlap": 5,
    "minIdentity": 0.9
  }
}
```

**Response:**
```json
{
  "overlap": null
}
```

## Worked Example

The longest qualifying suffix/prefix match is `CCCCC` (length 5), starting at position 5 in sequence1 and position 0 in sequence2.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L221](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L221)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
