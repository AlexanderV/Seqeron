# select_phylogenetic_markers

Select single-copy core clusters as phylogenetic markers.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `select_phylogenetic_markers` |
| **Method ID** | `PanGenomeAnalyzer.SelectPhylogeneticMarkers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Filters the candidate core clusters to **single-copy core** clusters — present in all
`totalGenomes` genomes with exactly one gene per genome (panX "all strains represented exactly
once"; Roary paralog filtering) — that contain at least one **parsimony-informative site**
(a column with ≥ 2 states each supported by ≥ 2 sequences, Zvelebil & Baum 2008). Qualifying
markers are returned ordered by descending parsimony-informative-site count (ties broken by
ordinal cluster id), capped at `maxMarkers`.

## Core Documentation Reference

- Source: [PanGenomeAnalyzer.cs#L997](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs#L997)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genomes` | GenomeInput[] | Yes | Genomes (used to recover member sequences). |
| `coreClusters` | GeneCluster[] | Yes | Candidate core clusters. |
| `totalGenomes` | integer | Yes | Total number of genomes. |
| `maxMarkers` | integer | No | Maximum markers to return (default 100). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | GeneCluster[] | Selected marker clusters, most informative first. |

## Errors

None. Empty candidates returns an empty list.

## Example

**User Prompt:**
> From two single-copy core clusters, pick the phylogenetic markers.

`aHi` has 2 informative sites, `zLo` has 1:

**Response:**
```json
{ "items": [ { "clusterId": "aHi" }, { "clusterId": "zLo" } ] }
```

With `maxMarkers = 1` only `aHi` (the most informative) is returned.

## References

- [PanGenomeAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs) — `SelectPhylogeneticMarkers`
- Ding W. et al. (2018) NAR 46:e5 (panX); Page A.J. et al. (2015) Bioinformatics 31:3691 (Roary); Zvelebil & Baum (2008).

## See Also

- [core_gene_clusters](core_gene_clusters.md)
- [core_genome_alignment](core_genome_alignment.md)
