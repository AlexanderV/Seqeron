# accessory_genes

Summarise accessory gene clusters of a pan-genome.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `accessory_genes` |
| **Method ID** | `PanGenomeAnalyzer.AnalyzeAccessoryGenes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Accessory (shell/cloud) genes are the clusters present in **more than one but not all**
genomes (`1 < genomeCount < totalGenomes`). For each such cluster the tool reports the ids
of the genomes that carry it and its frequency `genomeCount / totalGenomes`
(Tettelin et al. 2005; Page et al. 2015, Roary).

Core clusters (present in all genomes) and unique clusters (present in exactly one genome)
are excluded.

## Core Documentation Reference

- Source: [PanGenomeAnalyzer.cs#L855](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs#L855)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `clusters` | GeneCluster[] | Yes | All gene clusters (e.g. from `cluster_genes`). |
| `totalGenomes` | integer | Yes | Total number of genomes in the analysis. |

Each `GeneCluster` carries `clusterId`, `geneIds`, `genomeIds`, `genomeCount`,
`averageIdentity`, and `consensusSequence`.

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | One entry per accessory cluster. |
| `items[].clusterId` | string | Accessory cluster id. |
| `items[].genomesWithGene` | string[] | Ids of genomes carrying this cluster. |
| `items[].frequency` | number | `genomeCount / totalGenomes`. |

## Errors

None (empty or fully core/unique inputs simply return an empty list).

## Example

**User Prompt:**
> Which gene clusters are accessory across my 3 genomes?

Given clusters with occupancies 3/3, 2/3, and 1/3:

**Response:**
```json
{
  "items": [
    { "clusterId": "c2", "genomesWithGene": ["genome1", "genome2"], "frequency": 0.6666666666666666 }
  ]
}
```

Only the 2-of-3 cluster is accessory; the 3/3 cluster is core and the 1/3 cluster is unique.

## References

- [PanGenomeAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs) — `AnalyzeAccessoryGenes`
- Tettelin H. et al. (2005) PNAS 102:13950; Page A.J. et al. (2015) Bioinformatics 31:3691 (Roary).

## See Also

- [core_gene_clusters](core_gene_clusters.md) — the core partition
- [find_genome_specific_genes](find_genome_specific_genes.md) — the unique (singleton) partition
- [construct_pangenome](construct_pangenome.md) — full core/accessory/unique partition
