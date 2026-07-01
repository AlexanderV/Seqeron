# find_orthologs

Best-hit ortholog pairs between two genomes.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_orthologs` |
| **Method ID** | `ComparativeGenomics.FindOrthologs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **ortholog pairs** between two genomes by 5-mer Jaccard similarity of the gene
sequences, keeping mutual best hits above `minIdentity` and `minCoverage`. Each gene's
`sequence` must be populated. Each pair reports the two gene ids, the identity,
coverage and alignment length.

## Core Documentation Reference

- Source: [ComparativeGenomics.cs#L335](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs#L335)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genome1Genes` | array of Gene | Yes | Genes of genome 1 (with sequence) |
| `genome2Genes` | array of Gene | Yes | Genes of genome 2 (with sequence) |
| `minIdentity` | number | No | Minimum identity (default 0.3) |
| `minCoverage` | number | No | Minimum coverage (default 0.5) |

A Gene is `{ id, genomeId, start, end, strand, sequence }`.

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ gene1Id, gene2Id, identity, coverage, alignmentLength }` |

## Errors

None (empty gene lists yield an empty result).

## Examples

### Example 1: Two mutual best hits

**User Prompt:**
> Find orthologs between two genomes whose genes have identical sequences pairwise.

**Expected Tool Call:**
```json
{
  "tool": "find_orthologs",
  "arguments": { "genome1Genes": [ { "id": "a1", "genomeId": "G1", "start": 0, "end": 14, "strand": "+", "sequence": "ACGTACGTACGTAC" }, { "id": "a2", "genomeId": "G1", "start": 0, "end": 16, "strand": "+", "sequence": "TTTTGGGGCCCCAAAA" } ], "genome2Genes": [ { "id": "b1", "genomeId": "G2", "start": 0, "end": 14, "strand": "+", "sequence": "ACGTACGTACGTAC" }, { "id": "b2", "genomeId": "G2", "start": 0, "end": 16, "strand": "+", "sequence": "TTTTGGGGCCCCAAAA" } ] }
}
```

**Response:**
```json
{ "items": [ { "gene1Id": "a1", "gene2Id": "b1", "identity": 1.0, "coverage": 1.0 }, { "gene1Id": "a2", "gene2Id": "b2", "identity": 1.0, "coverage": 1.0 } ] }
```
Identical sequences give Jaccard identity and coverage 1.0.

## Performance

- **Time Complexity:** O(|G1|·|G2|) similarity comparisons.
- **Space Complexity:** O(|G1|·|G2|).

## See Also

- [find_reciprocal_best_hits](find_reciprocal_best_hits.md)
- [find_syntenic_blocks](find_syntenic_blocks.md)
- [compare_genomes](compare_genomes.md)
