# calculate_coverage

Map each read to its best ungapped position on the reference (>= minOverlap matching bases) and return per-base coverage depth as an integer array of length reference.Length.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `calculate_coverage` |
| **Method ID** | `SequenceAssembler.CalculateCoverage` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Map each read to its best ungapped position on the reference (>= minOverlap matching bases) and return per-base coverage depth as an integer array of length reference.Length.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L772](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L772)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference sequence. |
| `reads` | array<string> | Yes | Reads to map onto reference. |
| `minOverlap` | integer | No | Minimum matching bases required to map a read. (default: 20) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `coverage` | array<integer> | Per-base coverage depth; coverage[i] is the number of placed reads spanning reference position i. |

## Examples

### Example 1: Two reads

**Tool Call:**
```json
{
  "tool": "calculate_coverage",
  "arguments": {
    "reference": "ACGTACGTAC",
    "reads": [
      "ACGTA",
      "GTACG"
    ],
    "minOverlap": 5
  }
}
```

**Response:**
```json
{
  "coverage": [
    1,
    1,
    2,
    2,
    2,
    1,
    1,
    0,
    0,
    0
  ]
}
```

### Example 2: No reads

**Tool Call:**
```json
{
  "tool": "calculate_coverage",
  "arguments": {
    "reference": "ACGTACGTAC",
    "reads": [],
    "minOverlap": 5
  }
}
```

**Response:**
```json
{
  "coverage": [
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  ]
}
```

## Worked Example

Read `ACGTA` places at position 0 (covers [0,5)); read `GTACG` places at position 2 (covers [2,7)). Depths sum per position.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L772](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L772)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
