# bootstrap_support

Non-parametric bootstrap support for the clades of a phylogenetic tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `bootstrap_support` |
| **Method ID** | `PhylogeneticAnalyzer.Bootstrap` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Estimates clade support via Felsenstein's non-parametric bootstrap (Felsenstein 1985). A reference tree is built from the input alignment; then the alignment columns are resampled with replacement to a pseudo-alignment of the same length, a replicate tree is rebuilt, and each non-trivial clade of the reference tree is scored by the proportion of replicates that recover it. Support values lie in `[0, 1]` (multiply by 100 for the published percentage). The internal RNG is seeded with a fixed constant (42), so results are deterministic for identical inputs. Cost is `O(replicates · n² · L)`.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L1184](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L1184)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | object (map string→string) | Yes | Aligned sequences keyed by taxon name (≥2, equal length) |
| `replicates` | integer | No | Number of bootstrap replicates (≥1, default 100) |
| `distanceMethod` | string | No | `PDistance` \| `JukesCantor` \| `Kimura2Parameter` \| `Hamming` (default `JukesCantor`) |
| `treeMethod` | string | No | `UPGMA` \| `NeighborJoining` (default `UPGMA`) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `support` | array | One entry per non-trivial reference clade |
| `support[].clade` | string | Sorted, `\|`-joined leaf names of the clade |
| `support[].support` | number | Proportion of replicates recovering the clade, in `[0, 1]` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | At least 2 sequences required |
| 1003 | At least 1 replicate required |
| 1004 | Unknown distance method |
| 1005 | Unknown tree method |

## Example

**User Prompt:**
> Bootstrap the support for a tree of A=B (all A) and C=D (all G) with 100 UPGMA replicates.

**Tool Call:**
```json
{
  "tool": "bootstrap_support",
  "arguments": {
    "sequences": { "A": "AAAAAAAAAA", "B": "AAAAAAAAAA", "C": "GGGGGGGGGG", "D": "GGGGGGGGGG" },
    "replicates": 100,
    "distanceMethod": "JukesCantor",
    "treeMethod": "UPGMA"
  }
}
```

**Response:**
```json
{
  "support": [
    { "clade": "A|B", "support": 1.0 },
    { "clade": "C|D", "support": 1.0 }
  ]
}
```

Because the two groups are invariant, column resampling never changes the distance matrix, so both clades are recovered in every replicate (support = 1.0). With `treeMethod = NeighborJoining` the tree is an unrooted trifurcation `((A,B),C,D)`, whose only non-trivial rooted clade is `{A,B}`.

## References

- Felsenstein J (1985). *Confidence limits on phylogenies: an approach using the bootstrap.* Evolution 39(4):783–791.
- [PhylogeneticAnalyzer.Bootstrap](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L1184)

## See Also

- [build_phylogenetic_tree](build_phylogenetic_tree.md)
- [robinson_foulds_distance](robinson_foulds_distance.md)
