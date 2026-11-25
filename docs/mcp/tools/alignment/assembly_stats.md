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
| `assembledReads` | integer | No | Reads that lie in the contigs, in `[0, totalReads]`; default `totalReads` (0 when there are no contigs). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `contigs` | array<string> | The input contigs (echoed). |
| `totalReads` | integer | Total number of input reads. |
| `assembledReads` | integer | Reads that lie in the contigs (Velvet `usedReads`); the supplied `assembledReads`, else `totalReads` (0 without contigs). |
| `n50` | number | N50 contig length (QUAST `N50`: shortest contig of the fewest longest contigs covering ≥ 50 % of the total; 0 without contigs). |
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

Lengths 10, 20, 30 (total 60). Sorted descending, the 30-bp contig alone covers 30·100 ≥ 60·50, so N50 = 30.
Odd total: lengths 4, 3, 2 (total 9): 4·100 < 9·50, 7·100 ≥ 9·50 → N50 = 3 (QUAST `N50.py`; the
former `total/2 = 4` integer cut-off reported 4). The N50 comes from the single canonical Nx core
`SequenceAssembler.CalculateNx`, shared with `GenomeAssemblyAnalyzer.CalculateNx`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L574](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L574)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
