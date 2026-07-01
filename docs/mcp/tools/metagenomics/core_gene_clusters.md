# core_gene_clusters

Filter gene clusters down to the core set.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `core_gene_clusters` |
| **Method ID** | `PanGenomeAnalyzer.GetCoreGeneClusters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the clusters whose occupancy fraction `genomeCount / totalGenomes` is at least
`threshold` — the Roary "present in at least N % of samples" core definition (Page et al. 2015).
This is a **fractional** test, not `floor(threshold × totalGenomes)`: at threshold 0.99 over 3
genomes only a 3/3 cluster qualifies, and over 100 genomes a 99/100 cluster is core but 98/100
is not.

## Core Documentation Reference

- Source: [PanGenomeAnalyzer.cs#L810](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs#L810)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `clusters` | GeneCluster[] | Yes | All gene clusters. |
| `totalGenomes` | integer | Yes | Total number of genomes. |
| `threshold` | number | No | Core fraction threshold (default 0.99). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | GeneCluster[] | The clusters that meet the core threshold. |

## Errors

None. Empty input returns an empty list.

## Example

**User Prompt:**
> Which clusters are core across 3 genomes at 100 % occupancy?

Clusters with occupancy 3/3, 2/3, 1/3:

**Response:**
```json
{ "items": [ { "clusterId": "c1", "genomeCount": 3 } ] }
```

Only the 3/3 cluster meets `occupancy / 3 ≥ 1.0`.

## References

- [PanGenomeAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs) — `GetCoreGeneClusters`
- Page A.J. et al. (2015) Bioinformatics 31:3691 (Roary).

## See Also

- [accessory_genes](accessory_genes.md)
- [core_genome_alignment](core_genome_alignment.md)
- [select_phylogenetic_markers](select_phylogenetic_markers.md)
