# cluster_genes

Cluster genes from multiple genomes into ortholog groups.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `cluster_genes` |
| **Method ID** | `PanGenomeAnalyzer.ClusterGenes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Groups genes into ortholog clusters with the CD-HIT greedy incremental algorithm
(Li & Godzik 2006): sequences are processed longest → shortest, the longest becomes the
representative of the first cluster, and each remaining gene joins the first existing
representative whose global sequence identity meets `identityThreshold`, otherwise it starts a
new cluster. Global identity = identical residues in the ungapped alignment divided by the
length of the shorter sequence (CD-HIT `-G 1`). Each cluster reports its member gene/genome
ids, mean pairwise identity, and the representative (longest) member as consensus.

## Core Documentation Reference

- Source: [PanGenomeAnalyzer.cs#L210](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs#L210)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genomes` | GenomeInput[] | Yes | Genomes with ordered gene lists. |
| `identityThreshold` | number | No | Global identity cutoff in [0,1] (CD-HIT -c, default 0.9). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].clusterId` | string | Cluster id (`cluster_1`, …). |
| `items[].geneIds` | string[] | Member gene ids (representative first). |
| `items[].genomeIds` | string[] | Distinct genomes contributing members. |
| `items[].genomeCount` | integer | Number of distinct genomes. |
| `items[].averageIdentity` | number | Mean pairwise global identity. |
| `items[].consensusSequence` | string | Representative (longest) member sequence. |

## Errors

None. Empty input returns an empty list.

## Example

**User Prompt:**
> Cluster genes `ATGCATGC` and `ATGCATGG` (87.5 % identical) at threshold 0.8.

**Response:**
```json
{ "items": [ { "genomeCount": 2, "averageIdentity": 0.875 } ] }
```

They meet the 0.8 cutoff, so they form a single cluster with mean identity 7/8 = 0.875.

## References

- [PanGenomeAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs) — `ClusterGenes`
- Li W. & Godzik A. (2006) Bioinformatics 22:1658 (CD-HIT); Page A.J. et al. (2015) Bioinformatics 31:3691 (Roary).

## See Also

- [gene_presence_absence_matrix](gene_presence_absence_matrix.md)
- [construct_pangenome](construct_pangenome.md)
