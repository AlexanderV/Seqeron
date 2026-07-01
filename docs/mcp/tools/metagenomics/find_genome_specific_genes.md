# find_genome_specific_genes

List the cluster ids unique to each genome.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `find_genome_specific_genes` |
| **Method ID** | `PanGenomeAnalyzer.FindGenomeSpecificGenes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns, per genome, the ids of the singleton clusters (`genomeCount == 1`) owned only by that
genome — the strain-specific (unique) portion of the pan-genome. Genomes with no unique clusters
are omitted from the output.

## Core Documentation Reference

- Source: [PanGenomeAnalyzer.cs#L870](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs#L870)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genomes` | GenomeInput[] | Yes | Genomes with ordered gene lists. |
| `clusters` | GeneCluster[] | Yes | Gene clusters covering the genomes. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].genomeId` | string | Genome id. |
| `items[].uniqueGeneIds` | string[] | Singleton cluster ids owned only by this genome. |

## Errors

None. No clusters returns an empty list.

## Example

**User Prompt:**
> Which clusters are strain-specific in g1 and g2?

g1 owns cluster c1, g2 owns c2 (both singletons):

**Response:**
```json
{ "items": [ { "genomeId": "g1", "uniqueGeneIds": ["c1"] }, { "genomeId": "g2", "uniqueGeneIds": ["c2"] } ] }
```

## References

- [PanGenomeAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs) — `FindGenomeSpecificGenes`
- Tettelin H. et al. (2005) PNAS 102:13950.

## See Also

- [accessory_genes](accessory_genes.md)
- [construct_pangenome](construct_pangenome.md)
