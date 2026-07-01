# find_reciprocal_best_hits

Reciprocal best hits (RBH) between two genomes.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_reciprocal_best_hits` |
| **Method ID** | `ComparativeGenomics.FindReciprocalBestHits` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **reciprocal best hits** (RBH) between two genomes — a stricter ortholog
criterion than one-directional best hits: a gene pair is returned only when each gene
is the other's best hit (above `minIdentity`). Each gene's `sequence` must be
populated.

## Core Documentation Reference

- Source: [ComparativeGenomics.cs#L466](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs#L466)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genome1Genes` | array of Gene | Yes | Genes of genome 1 (with sequence) |
| `genome2Genes` | array of Gene | Yes | Genes of genome 2 (with sequence) |
| `minIdentity` | number | No | Minimum identity (default 0.3) |

A Gene is `{ id, genomeId, start, end, strand, sequence }`.

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ gene1Id, gene2Id, identity, coverage, alignmentLength }` |

## Errors

None (empty gene lists yield an empty result).

## Examples

### Example 1: Two reciprocal best hits

**User Prompt:**
> Find reciprocal best hits between two genomes with pairwise-identical genes.

**Expected Tool Call:**
```json
{
  "tool": "find_reciprocal_best_hits",
  "arguments": { "genome1Genes": [ { "id": "a1", "genomeId": "G1", "start": 0, "end": 14, "strand": "+", "sequence": "ACGTACGTACGTAC" }, { "id": "a2", "genomeId": "G1", "start": 0, "end": 16, "strand": "+", "sequence": "TTTTGGGGCCCCAAAA" } ], "genome2Genes": [ { "id": "b1", "genomeId": "G2", "start": 0, "end": 14, "strand": "+", "sequence": "ACGTACGTACGTAC" }, { "id": "b2", "genomeId": "G2", "start": 0, "end": 16, "strand": "+", "sequence": "TTTTGGGGCCCCAAAA" } ] }
}
```

**Response:**
```json
{ "items": [ { "gene1Id": "a1", "gene2Id": "b1" }, { "gene1Id": "a2", "gene2Id": "b2" } ] }
```

## Performance

- **Time Complexity:** O(|G1|·|G2|).
- **Space Complexity:** O(|G1|·|G2|).

## See Also

- [find_orthologs](find_orthologs.md)
- [find_syntenic_blocks](find_syntenic_blocks.md)
