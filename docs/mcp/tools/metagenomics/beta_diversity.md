# beta_diversity

Compute beta-diversity (between-sample) distances for two samples.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `beta_diversity` |
| **Method ID** | `MetagenomicsAnalyzer.CalculateBetaDiversity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Given two per-species abundance vectors, computes:

- **Bray-Curtis** dissimilarity `1 - 2·Σ min(aᵢ, bᵢ) / Σ (aᵢ + bᵢ)` (Bray & Curtis 1957)
- **Jaccard distance** `1 - |shared| / |union|` on presence/absence (Jaccard 1901)
- **Shared / unique species** counts (a species is "present" when its abundance is `> 0`)

UniFrac is reported as `0` because no phylogenetic tree is supplied.

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L572](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L572)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sample1Name` | string | Yes | Identifier of sample 1. |
| `sample1` | AbundanceItem[] | Yes | Sample 1 abundance vector. |
| `sample2Name` | string | Yes | Identifier of sample 2. |
| `sample2` | AbundanceItem[] | Yes | Sample 2 abundance vector. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sample1` / `sample2` | string | Echoed sample names. |
| `brayCurtis` | number | Bray-Curtis dissimilarity. |
| `jaccardDistance` | number | Jaccard distance on presence/absence. |
| `uniFracDistance` | number | Always 0 (no tree). |
| `sharedSpecies` | integer | Species present in both. |
| `uniqueToSample1` | integer | Species present only in sample 1. |
| `uniqueToSample2` | integer | Species present only in sample 2. |

## Errors

None. Empty vectors yield zero distances.

## Example

**User Prompt:**
> Compare the fish communities of my two tanks.

Tank 1 = Goldfish 6, Guppy 7, Rainbow 4; Tank 2 = Goldfish 10, Guppy 0, Rainbow 6:

**Response:**
```json
{
  "sample1": "Tank 1",
  "sample2": "Tank 2",
  "brayCurtis": 0.3939393939393939,
  "jaccardDistance": 0.3333333333333333,
  "uniFracDistance": 0,
  "sharedSpecies": 2,
  "uniqueToSample1": 1,
  "uniqueToSample2": 0
}
```

Bray-Curtis = 1 − 2·10/33 = 13/33; Jaccard = 1 − 2/3 = 1/3 (Guppy absent in Tank 2).

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `CalculateBetaDiversity`
- Bray J.R. & Curtis J.T. (1957) Ecol. Monogr. 27:325; Jaccard P. (1901); Whittaker R.H. (1960).

## See Also

- [alpha_diversity](alpha_diversity.md) — within-sample diversity
- [differential_abundance](differential_abundance.md) — per-taxon condition comparison
