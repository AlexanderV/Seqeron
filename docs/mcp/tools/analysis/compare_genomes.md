# compare_genomes

End-to-end comparative-genomics pipeline for two gene sets.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `compare_genomes` |
| **Method ID** | `ComparativeGenomics.CompareGenomes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Runs the full comparative pipeline over two gene sets: reciprocal-best-hit (RBH)
ortholog detection, syntenic-block identification, rearrangement detection, and
summary statistics. Returns the syntenic blocks, ortholog pairs, rearrangement
events, the overall syntenic-gene fraction, the number of conserved (core) genes,
and each genome's specific (dispensable) gene count. The core/dispensable partition
follows Tettelin et al. (2005); orthology is by RBH (Moreno-Hagelsieb & Latimer 2008).

## Core Documentation Reference

- Source: [ComparativeGenomics.cs#L766](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs#L766)
- Evidence: `docs/Evidence/COMPGEN-COMPARE-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genome1Genes` | Gene[] | Yes | Genes of genome 1 (≥ 1) |
| `genome2Genes` | Gene[] | Yes | Genes of genome 2 (≥ 1) |
| `minOrthologIdentity` | number | No | Minimum ortholog identity (default 0.3) |
| `minSyntenicBlockSize` | integer | No | Minimum syntenic block size (default 3) |

Each **Gene** is `{ id, genomeId, start, end, strand, sequence? }`.

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `syntenicBlocks` | array | Syntenic blocks between the two genomes |
| `orthologs` | array | RBH ortholog gene pairs |
| `rearrangements` | array | Detected rearrangement events |
| `overallSynteny` | number | Fraction of syntenic genes in [0, 1] |
| `conservedGenes` | integer | Number of core genes (RBH ortholog pairs) |
| `genomeSpecificGenes1` | integer | Dispensable genes unique to genome 1 |
| `genomeSpecificGenes2` | integer | Dispensable genes unique to genome 2 |

## Errors

| Code | Message |
|------|---------|
| 1001 | Genome must contain at least one gene |

## Examples

### Example 1: One shared + one unique each → core/dispensable partition

**User Prompt:**
> Compare two genomes that share one gene and each have one unique gene.

**Expected Tool Call (sequences abbreviated):**
```json
{
  "tool": "compare_genomes",
  "arguments": {
    "genome1Genes": [
      { "id": "a0", "genomeId": "G1", "start": 0,   "end": 60,  "strand": "+", "sequence": "<shared>" },
      { "id": "a1", "genomeId": "G1", "start": 100, "end": 160, "strand": "+", "sequence": "<unique1>" }
    ],
    "genome2Genes": [
      { "id": "c0", "genomeId": "G2", "start": 0,   "end": 60,  "strand": "+", "sequence": "<shared>" },
      { "id": "c1", "genomeId": "G2", "start": 100, "end": 160, "strand": "+", "sequence": "<unique2>" }
    ]
  }
}
```

**Response (partition fields):**
```json
{ "conservedGenes": 1, "genomeSpecificGenes1": 1, "genomeSpecificGenes2": 1 }
```
The shared gene is core; each unique gene is dispensable to its genome (Tettelin 2005).

### Example 2: Disjoint content → no core

Two genomes with mutually dissimilar genes give
`{ "conservedGenes": 0, "genomeSpecificGenes1": 2, "genomeSpecificGenes2": 2 }`.

## Performance

- **Time Complexity:** O(n·m) pairwise gene comparison for RBH.

## See Also

- [find_orthologs](find_orthologs.md)
- [find_syntenic_blocks](find_syntenic_blocks.md)
- [calculate_ani](calculate_ani.md)
