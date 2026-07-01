# core_genome_alignment

Build a per-genome core-genome alignment block.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `core_genome_alignment` |
| **Method ID** | `PanGenomeAnalyzer.CreateCoreGenomeAlignment` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For the requested genome, concatenates — in cluster order — the sequence of the first member
gene that the genome contributes to each supplied core cluster. This produces the genome's row
of a concatenated core-genome alignment (as used for core-genome phylogenies). If the genome id
is not present, an empty string is returned.

## Core Documentation Reference

- Source: [PanGenomeAnalyzer.cs#L824](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs#L824)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genomes` | GenomeInput[] | Yes | Genomes with ordered gene lists. |
| `coreClusters` | GeneCluster[] | Yes | Core clusters (e.g. from `core_gene_clusters`). |
| `genomeId` | string | Yes | Genome to extract the alignment block for. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `result` | string | Concatenated core-cluster sequences (`""` if the genome is absent). |

## Errors

None. An unknown genome id yields an empty string.

## Example

**User Prompt:**
> Build g1's core-genome alignment block from clusters c1 (gene1 = ATGC) and c2 (gene2 = GCTA).

**Response:**
```json
{ "result": "ATGCGCTA" }
```

## References

- [PanGenomeAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs) — `CreateCoreGenomeAlignment`
- Page A.J. et al. (2015) Bioinformatics 31:3691 (Roary).

## See Also

- [core_gene_clusters](core_gene_clusters.md)
- [select_phylogenetic_markers](select_phylogenetic_markers.md)
