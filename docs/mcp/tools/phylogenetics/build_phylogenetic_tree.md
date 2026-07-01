# build_phylogenetic_tree

Build a phylogenetic tree from a set of named, pre-aligned sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `build_phylogenetic_tree` |
| **Method ID** | `PhylogeneticAnalyzer.BuildTree` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes a pairwise distance matrix using the chosen substitution model (`PDistance` | `JukesCantor` | `Kimura2Parameter` | `Hamming`) and constructs a tree with **UPGMA** (rooted, ultrametric, strictly bifurcating) or **Neighbor-Joining** (unrooted; the central node is a trifurcation of the last three OTUs, Saitou & Nei 1987). Returns the Newick string, the taxa, the distance matrix, and the method name. All input sequences must be the same length (aligned).

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L136](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L136)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | object (map string→string) | Yes | Aligned sequences keyed by taxon name (≥2, equal length) |
| `distanceMethod` | string | No | `PDistance` \| `JukesCantor` \| `Kimura2Parameter` \| `Hamming` (default `JukesCantor`) |
| `treeMethod` | string | No | `UPGMA` \| `NeighborJoining` (default `UPGMA`) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `newick` | string | Tree in Newick format (ends with `;`) |
| `taxa` | array of string | Taxon names |
| `distanceMatrix` | number[][] | Symmetric pairwise distance matrix; diagonal is zero |
| `method` | string | `UPGMA` or `NeighborJoining` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | At least 2 sequences required |
| 1003 | All sequences must have the same length (aligned) |
| 1004 | Unknown distance method |
| 1005 | Unknown tree method |

## Example

**User Prompt:**
> Build a UPGMA tree for A and B (identical) plus C.

**Tool Call:**
```json
{
  "tool": "build_phylogenetic_tree",
  "arguments": {
    "sequences": { "A": "ACGTACGT", "B": "ACGTACGT", "C": "TTTTTTTT" },
    "distanceMethod": "JukesCantor",
    "treeMethod": "UPGMA"
  }
}
```

**Response (excerpt):**
```json
{
  "taxa": ["A", "B", "C"],
  "method": "UPGMA",
  "distanceMatrix": [[0, 0, "Infinity"], [0, 0, "Infinity"], ["Infinity", "Infinity", 0]]
}
```

`distanceMatrix[0][1] = 0` because A and B are identical, so UPGMA joins them first. A vs C differs at 6 of 8 comparable sites (p = 0.75), which saturates the Jukes-Cantor correction to `+Infinity`.

## References

- Sokal R, Michener C (1958). UPGMA. — Saitou N, Nei M (1987). *The neighbor-joining method.* Mol Biol Evol 4(4):406–425.
- [PhylogeneticAnalyzer.BuildTree](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L136)

## See Also

- [build_tree_from_matrix](build_tree_from_matrix.md)
- [distance_matrix](distance_matrix.md)
- [to_newick](to_newick.md)
