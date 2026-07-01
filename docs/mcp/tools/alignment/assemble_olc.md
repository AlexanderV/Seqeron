# assemble_olc

Assemble reads with the overlap-layout-consensus approach: detect pairwise suffix-prefix overlaps above thresholds, build the overlap graph, greedily lay out contigs, and report N50 / longest / total length.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `assemble_olc` |
| **Method ID** | `SequenceAssembler.AssembleOLC` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assemble reads with the overlap-layout-consensus approach: detect pairwise suffix-prefix overlaps above thresholds, build the overlap graph, greedily lay out contigs, and report N50 / longest / total length.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L61](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L61)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reads` | array<string> | Yes | Input sequence reads. |
| `minOverlap` | integer | No | Minimum overlap length to accept (bp). (default: 20) |
| `minIdentity` | number | No | Minimum identity ratio for an overlap to be accepted (0.0-1.0). (default: 0.9) |
| `kmerSize` | integer | No | K-mer size (unused by OLC path; shared parameter). (default: 31) |
| `minContigLength` | integer | No | Minimum contig length to keep (bp). (default: 100) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `contigs` | array<string> | Assembled contig sequences. |
| `totalReads` | integer | Number of input reads. |
| `assembledReads` | integer | Number of reads incorporated. |
| `n50` | number | N50 contig length. |
| `longestContig` | integer | Longest contig length (bp). |
| `totalLength` | integer | Total assembled length (bp). |

## Errors

| Code | Message |
|------|---------|
| 4001 | Assembly failed |

## Examples

### Example 1: Two overlapping reads

**Tool Call:**
```json
{
  "tool": "assemble_olc",
  "arguments": {
    "reads": [
      "AAAAACCCCC",
      "CCCCCGGGGG"
    ],
    "minOverlap": 5,
    "minIdentity": 0.9,
    "kmerSize": 4,
    "minContigLength": 5
  }
}
```

**Response:**
```json
{
  "contigs": [
    "AAAAACCCCCGGGGG"
  ],
  "totalReads": 2,
  "assembledReads": 2,
  "n50": 15,
  "longestContig": 15,
  "totalLength": 15
}
```

### Example 2: Empty read set

**Tool Call:**
```json
{
  "tool": "assemble_olc",
  "arguments": {
    "reads": [],
    "minOverlap": 5,
    "minIdentity": 0.9,
    "kmerSize": 4,
    "minContigLength": 5
  }
}
```

**Response:**
```json
{
  "contigs": [],
  "totalReads": 0,
  "assembledReads": 0,
  "n50": 0,
  "longestContig": 0,
  "totalLength": 0
}
```

## Worked Example

Reads `AAAAACCCCC` and `CCCCCGGGGG` share a 5-bp `CCCCC` suffix/prefix overlap, so they merge into the single 15-bp contig `AAAAACCCCCGGGGG`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L61](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L61)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
