# assemble_de_bruijn

Assemble reads using a de Bruijn graph: shred reads into their distinct k-mers (case-insensitive; coverage adds no edges), build the (k-1)-mer node graph, output every maximal non-branching path (unitig) and isolated cycle as a contig, then report N50 / longest / total length. Optionally add reverse-complement k-mers (double-stranded reads).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `assemble_de_bruijn` |
| **Method ID** | `SequenceAssembler.AssembleDeBruijn` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assemble reads using a de Bruijn graph: shred reads into their distinct k-mers (case-insensitive; coverage adds no edges), build the (k-1)-mer node graph, output every maximal non-branching path (unitig) and isolated cycle as a contig, then report N50 / longest / total length. Optionally add reverse-complement k-mers (double-stranded reads).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L101](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L101)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reads` | array<string> | Yes | Input sequence reads. |
| `minOverlap` | integer | No | Minimum overlap length (shared parameter; not the operative knob here). (default: 20) |
| `minIdentity` | number | No | Minimum identity (shared parameter). (default: 0.9) |
| `kmerSize` | integer | No | K-mer size used to build the de Bruijn graph (>= 2; odd recommended with includeReverseComplement). (default: 31) |
| `minContigLength` | integer | No | Minimum contig length to keep (bp). (default: 100) |
| `includeReverseComplement` | boolean | No | Treat reads as double-stranded DNA: add each k-mer's reverse complement (Velvet). (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `contigs` | array<string> | Assembled contig sequences. |
| `totalReads` | integer | Number of input reads. |
| `assembledReads` | integer | Reads with at least one k-mer (either strand with `includeReverseComplement`) in a reported contig (Velvet `usedReads`). |
| `n50` | number | N50 contig length. |
| `longestContig` | integer | Longest contig length (bp). |
| `totalLength` | integer | Total assembled length (bp). |

## Errors

| Code | Message |
|------|---------|
| 4001 | Assembly failed |

## Examples

### Example 1: Overlapping reads

**Tool Call:**
```json
{
  "tool": "assemble_de_bruijn",
  "arguments": {
    "reads": [
      "ATGGCGT",
      "GGCGTGCA"
    ],
    "kmerSize": 4,
    "minContigLength": 4
  }
}
```

**Response:**
```json
{
  "contigs": [
    "ATGGCGTGCA"
  ],
  "totalReads": 2,
  "assembledReads": 2,
  "n50": 10,
  "longestContig": 10,
  "totalLength": 10
}
```

### Example 1b: Branching graph yields unitigs

`AAABBB`, `AABBBC` at k = 3: the nodes AA and BB carry self-loops (not 1-in-1-out), so the
contigs are `AAA`, `AABB`, `BBB`, `BBC` (totalLength 13, longestContig 4, n50 3).

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

With k = 4 the reads `ATGGCGT` and `GGCGTGCA` share the 4-mers `GGCG` and `GCGT`, which are
counted once; the distinct 4-mers form one non-branching path of 3-mer nodes that spells
`ATGGCGTGCA` (Compeau & Pevzner BA3K/BA3M; Velvet). With `includeReverseComplement`, reads
`GATTACAGGT` and `TAAGGACCTG` (k = 5) join into `GATTACAGGTCCTTA`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L101](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L101)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
