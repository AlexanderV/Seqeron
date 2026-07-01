# construct_pangenome

Construct a pan-genome from a set of genomes.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `construct_pangenome` |
| **Method ID** | `PanGenomeAnalyzer.ConstructPanGenome` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Clusters the genes of all genomes into ortholog groups, then partitions the clusters by
occupancy:

- **core** — present in at least `coreFraction` of the genomes (Roary fractional rule, Page 2015)
- **unique** — present in exactly one genome
- **accessory** — everything in between

It also reports **genome fluidity** (Kislyuk 2011) and an **open vs closed** classification
from the Heaps'-law decay exponent of the new-gene curve (Tettelin 2008; open when α < 1).

## Core Documentation Reference

- Source: [PanGenomeAnalyzer.cs#L103](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs#L103)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genomes` | GenomeInput[] | Yes | Genomes with ordered gene lists. |
| `identityThreshold` | number | No | Ortholog-clustering identity (default 0.9). |
| `coreFraction` | number | No | Core occupancy fraction (default 0.99). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `coreGenes` / `accessoryGenes` / `uniqueGenes` | string[] | Cluster ids per partition. |
| `genomeToGenes[]` | array | Per-genome gene id lists. |
| `statistics` | object | Counts, `coreFraction`, `genomeFluidity`, `type` (`Open`/`Closed`). |

## Errors

None. Empty/null input yields zeroed statistics and a `Closed` type.

## Example

**User Prompt:**
> Build the pan-genome of 3 genomes sharing one core gene, one in-two-of-three gene, and three
> strain-specific genes (coreFraction = 1.0).

**Response (statistics):**
```json
{ "totalGenomes": 3, "totalGenes": 5, "coreGeneCount": 1, "accessoryGeneCount": 1, "uniqueGeneCount": 3, "coreFraction": 0.2 }
```

1 core + 1 accessory + 3 unique = 5 clusters; coreFraction = 1/5 = 0.2.

## References

- [PanGenomeAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs) — `ConstructPanGenome`
- Tettelin H. et al. (2005) PNAS 102:13950; Tettelin H. et al. (2008) Curr. Opin. Microbiol. 11:472; Kislyuk A.O. et al. (2011) BMC Genomics 12:32; Page A.J. et al. (2015) Bioinformatics 31:3691.

## See Also

- [cluster_genes](cluster_genes.md)
- [core_gene_clusters](core_gene_clusters.md)
- [fit_heaps_law](fit_heaps_law.md)
