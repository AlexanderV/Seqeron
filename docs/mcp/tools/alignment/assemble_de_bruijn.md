# assemble_de_bruijn

Assemble reads using a de Bruijn graph: shred reads into their distinct k-mers (case-insensitive; coverage adds no edges; Velvet read rule: every non-ACGT symbol is read as A), build the (k-1)-mer node graph, output every maximal non-branching path (unitig) and isolated cycle as a contig, then report N50 / longest / total length. Double-stranded by default (Velvet: each k-mer is recorded with its reverse complement; each contig/reverse-complement pair is reported once); includeReverseComplement = false gives the single-strand textbook model.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `assemble_de_bruijn` |
| **Method ID** | `SequenceAssembler.AssembleDeBruijn` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assemble reads using a de Bruijn graph: shred reads into their distinct k-mers (case-insensitive; coverage adds no edges; Velvet read rule: every non-ACGT symbol is read as A), build the (k-1)-mer node graph, output every maximal non-branching path (unitig) and isolated cycle as a contig, then report N50 / longest / total length. Double-stranded by default (Velvet: each k-mer is recorded with its reverse complement; each contig/reverse-complement pair is reported once); includeReverseComplement = false gives the single-strand textbook model.

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
| `includeReverseComplement` | boolean | No | Treat reads as double-stranded DNA: add each k-mer's reverse complement (Velvet; default true). False = forward strand only (Rosalind BA3M single-strand model). (default: true) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `contigs` | array<string> | Assembled contig sequences. |
| `totalReads` | integer | Number of input reads. |
| `assembledReads` | integer | Reads (non-ACGT read as A) with at least one k-mer (either strand with `includeReverseComplement`) in a reported contig (Velvet `usedReads`). |
| `n50` | number | N50 contig length. |
| `longestContig` | integer | Longest contig length (bp). |
| `totalLength` | integer | Total assembled length (bp). |

## Errors

| Code | Message |
|------|---------|
| 4001 | Assembly failed |

## Examples

### Example 1: Overlapping reads (single strand)

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
    "minContigLength": 4,
    "includeReverseComplement": false
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

### Example 1b: Branching graph yields unitigs (single strand)

`AAACCC`, `AACCCG` at k = 3 with `includeReverseComplement: false`: the nodes AA and CC carry self-loops
(not 1-in-1-out), so the contigs are `AAA`, `AACC`, `CCC`, `CCG` (totalLength 13, longestContig 4, n50 3).

### Example 1c: Double-stranded default

`GATTACAGGT`, `TAAGGACCTG` at k = 5 (defaults): the second read is the reverse complement of the genome's
3' half, so the strands join into `GATTACAGGTCCTTA` (the pair {contig, reverse complement} is reported once, as
the ordinally smaller string); with `includeReverseComplement: false` the result is the two input reads.

### Example 1d: Non-ACGT symbols

`ACGTNGGTCA` at k = 4 with `includeReverseComplement: false` gives `ACGTAGGTCA`: as in Velvet's read loader
(`velvetifySequence` + `writeNucleotide`) every symbol other than A/C/G/T is read as A.

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
`ATGGCGTGCA` (Compeau & Pevzner BA3K/BA3M; Velvet; single strand). With the default
double-stranded graph, reads `GATTACAGGT` and `TAAGGACCTG` (k = 5) join into `GATTACAGGTCCTTA`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L101](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L101)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
