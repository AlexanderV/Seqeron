# assembly_stats

Compute assembly statistics (N50, longest contig, total length, read accounting) for a precomputed list of contigs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `assembly_stats` |
| **Method ID** | `SequenceAssembler.CalculateStats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Compute assembly statistics (N50, longest contig, total length, read accounting) for a precomputed list of contigs.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L574](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L574)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `contigs` | array<string> | Yes | Assembled contig sequences. |
| `totalReads` | integer | Yes | Total number of input reads (used to fill totalReads). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `contigs` | array<string> | The input contigs (echoed). |
| `totalReads` | integer | Total number of input reads. |
| `assembledReads` | integer | Number of reads incorporated (simplified: all reads). |
| `n50` | number | N50 contig length. |
| `longestContig` | integer | Longest contig length (bp). |
| `totalLength` | integer | Total contig length (bp). |

## Examples

### Example 1: Three contigs

**Tool Call:**
```json
{
  "tool": "assembly_stats",
  "arguments": {
    "contigs": [
      "AAAAAAAAAA",
      "AAAAAAAAAAAAAAAAAAAA",
      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    ],
    "totalReads": 5
  }
}
```

**Response:**
```json
{
  "contigs": [
    "AAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  ],
  "totalReads": 5,
  "assembledReads": 5,
  "n50": 30,
  "longestContig": 30,
  "totalLength": 60
}
```

### Example 2: No contigs

**Tool Call:**
```json
{
  "tool": "assembly_stats",
  "arguments": {
    "contigs": [],
    "totalReads": 3
  }
}
```

**Response:**
```json
{
  "contigs": [],
  "totalReads": 3,
  "assembledReads": 0,
  "n50": 0,
  "longestContig": 0,
  "totalLength": 0
}
```

## Worked Example

Lengths 10, 20, 30 (total 60, half 30). Sorted descending, cumulative length reaches the half at the 30-bp contig, so N50 = 30.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L574](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L574)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
