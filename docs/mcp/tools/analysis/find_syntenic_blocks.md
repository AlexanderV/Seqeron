# find_syntenic_blocks

Collinear runs of orthologous genes (syntenic blocks).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_syntenic_blocks` |
| **Method ID** | `ComparativeGenomics.FindSyntenicBlocks` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **syntenic blocks** — collinear runs of orthologous genes — between two genomes,
given an ortholog map (`genome1 gene id → genome2 gene id`). Anchors that are adjacent
(within `maxGap`) and consistently ordered are chained; a block must contain at least
`minBlockSize` anchors. Each block reports its coordinates in both genomes, gene count,
whether it is inverted, and its identity.

## Core Documentation Reference

- Source: [ComparativeGenomics.cs#L129](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs#L129)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genome1Genes` | array of Gene | Yes | Genes of genome 1 |
| `genome2Genes` | array of Gene | Yes | Genes of genome 2 |
| `orthologMap` | object | Yes | `genome1 gene id → genome2 gene id` |
| `minBlockSize` | integer | No | Minimum block size (default 3) |
| `maxGap` | integer | No | Maximum gap (default 5) |

A Gene is `{ id, genomeId, start, end, strand, sequence? }`.

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Blocks `{ genome1Id, start1, end1, genome2Id, start2, end2, isInverted, geneCount, identity }` |

## Errors

None (empty inputs yield an empty result).

## Examples

### Example 1: Five collinear anchors

**User Prompt:**
> Find syntenic blocks for 5 genes g0..g4 mapping 1:1 to h0..h4 in the same order.

**Expected Tool Call:**
```json
{
  "tool": "find_syntenic_blocks",
  "arguments": { "genome1Genes": "g0..g4 at 0,100,…", "genome2Genes": "h0..h4 at 0,100,…", "orthologMap": { "g0": "h0", "g1": "h1", "g2": "h2", "g3": "h3", "g4": "h4" } }
}
```

**Response:**
```json
{ "items": [ { "geneCount": 5, "isInverted": false } ] }
```
Five adjacent collinear anchors form exactly one forward block.

## Performance

- **Time Complexity:** O(n log n) anchor chaining.
- **Space Complexity:** O(number of anchors).

## See Also

- [find_orthologs](find_orthologs.md)
- [detect_rearrangements](detect_rearrangements.md)
- [reversal_distance](reversal_distance.md)
