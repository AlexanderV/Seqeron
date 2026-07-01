# assemble_de_bruijn

Assemble reads using a de Bruijn graph: shred reads into k-mers, build the (k-1)-mer node graph, trace Eulerian walks per component into contigs, then report N50 / longest / total length.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `assemble_de_bruijn` |
| **Method ID** | `SequenceAssembler.AssembleDeBruijn` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assemble reads using a de Bruijn graph: shred reads into k-mers, build the (k-1)-mer node graph, trace Eulerian walks per component into contigs, then report N50 / longest / total length.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L104](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L104)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reads` | array<string> | Yes | Input sequence reads. |
| `minOverlap` | integer | No | Minimum overlap length (shared parameter; not the operative knob here). (default: 20) |
| `minIdentity` | number | No | Minimum identity (shared parameter). (default: 0.9) |
| `kmerSize` | integer | No | K-mer size used to build the de Bruijn graph. (default: 31) |
| `minContigLength` | integer | No | Minimum contig length to keep (bp). (default: 100) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `contigs` | array<string> | Assembled contig sequences. |
| `totalReads` | integer | Number of input reads. |
| `assembledReads` | integer | Number of reads incorporated (simplified: all reads). |
| `n50` | number | N50 contig length. |
| `longestContig` | integer | Longest contig length (bp). |
| `totalLength` | integer | Total assembled length (bp). |

## Errors

| Code | Message |
|------|---------|
| 4001 | Assembly failed |

## Examples

### Example 1: Two-read Eulerian walk

**Tool Call:**
```json
{
  "tool": "assemble_de_bruijn",
  "arguments": {
    "reads": [
      "AAABBB",
      "AABBBC"
    ],
    "kmerSize": 3,
    "minContigLength": 3
  }
}
```

**Response:**
```json
{
  "contigs": [
    "AAABBBBBBC"
  ],
  "totalReads": 2,
  "assembledReads": 2,
  "n50": 10,
  "longestContig": 10,
  "totalLength": 10
}
```

### Example 2: Empty read set

**Tool Call:**
```json
{
  "tool": "assemble_de_bruijn",
  "arguments": {
    "reads": [],
    "kmerSize": 3,
    "minContigLength": 3
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

With k=3 the graph nodes are 2-mers; the Eulerian walk over reads `AAABBB` and `AABBBC` spells the 10-bp superstring `AAABBBBBBC`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L104](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L104)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
