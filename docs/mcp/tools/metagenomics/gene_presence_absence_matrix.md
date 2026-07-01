# gene_presence_absence_matrix

Build a per-genome gene presence/absence matrix.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `gene_presence_absence_matrix` |
| **Method ID** | `PanGenomeAnalyzer.CreatePresenceAbsenceMatrix` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Produces one row per genome. For each supplied cluster (the columns), the cluster is marked
**present** in a genome when the genome owns any of the cluster's gene ids. Each row carries the
flattened `(clusterId → present)` entries plus the total column count and the number present.

## Core Documentation Reference

- Source: [PanGenomeAnalyzer.cs#L352](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs#L352)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genomes` | GenomeInput[] | Yes | Genomes with ordered gene lists. |
| `clusters` | GeneCluster[] | Yes | Clusters forming the columns (e.g. from `cluster_genes`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].genomeId` | string | Genome id (row). |
| `items[].genePresence[]` | array | `{ clusterId, present }` per cluster. |
| `items[].totalGenes` | integer | Number of clusters (columns). |
| `items[].presentGenes` | integer | Clusters present in this genome. |

## Errors

None.

## Example

**User Prompt:**
> Build the presence/absence matrix for g1 {a,b} and g2 {c} against clusters c1 {a,c} and c2 {b}.

**Response:**
```json
{
  "items": [
    { "genomeId": "g1", "totalGenes": 2, "presentGenes": 2 },
    { "genomeId": "g2", "totalGenes": 2, "presentGenes": 1 }
  ]
}
```

g1 carries both clusters; g2 carries only c1 (via gene c).

## References

- [PanGenomeAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs) — `CreatePresenceAbsenceMatrix`
- Page A.J. et al. (2015) Bioinformatics 31:3691 (Roary).

## See Also

- [cluster_genes](cluster_genes.md)
- [fit_heaps_law](fit_heaps_law.md)
