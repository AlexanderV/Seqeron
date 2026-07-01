# build_tree_from_matrix

Build a phylogenetic tree directly from a precomputed symmetric distance matrix.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `build_tree_from_matrix` |
| **Method ID** | `PhylogeneticAnalyzer.BuildTreeFromMatrix` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Constructs a tree (UPGMA or Neighbor-Joining) directly from a symmetric square distance matrix, bypassing sequence-to-distance computation. Useful for verifying tree-construction algorithms against reference matrices (e.g. the Wikipedia worked examples). The matrix must be square with size equal to the number of taxa, and the taxa order must match the matrix rows/columns.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L174](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L174)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `taxa` | array of string | Yes | Taxon names in matrix order (≥2) |
| `distanceMatrix` | number[][] | Yes | Symmetric square matrix; size = `taxa.length` |
| `treeMethod` | string | No | `UPGMA` \| `NeighborJoining` (default `UPGMA`) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `newick` | string | Tree in Newick format |
| `taxa` | array of string | Taxon names |
| `distanceMatrix` | number[][] | The input matrix, echoed back |
| `method` | string | `UPGMA` or `NeighborJoining` |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least 2 taxa required |
| 1002 | Distance matrix is required |
| 1003 | Distance matrix dimensions must match the number of taxa |
| 1004 | Unknown tree method |

## Example

**User Prompt:**
> Build a UPGMA tree from the Wikipedia 5S rRNA distance matrix.

**Tool Call:**
```json
{
  "tool": "build_tree_from_matrix",
  "arguments": {
    "taxa": ["a", "b", "c", "d", "e"],
    "distanceMatrix": [
      [0, 17, 21, 31, 23],
      [17, 0, 30, 34, 21],
      [21, 30, 0, 28, 39],
      [31, 34, 28, 0, 43],
      [23, 21, 39, 43, 0]
    ],
    "treeMethod": "UPGMA"
  }
}
```

**Response (excerpt):**
```json
{
  "newick": "(((a:8.5000,b:8.5000):2.5000,e:11.0000):5.5000,(c:14.0000,d:14.0000):2.5000);",
  "method": "UPGMA"
}
```

This matches the hand-derived ultrametric tree of the Wikipedia UPGMA working example (all tips 16.5 from the root). With `NeighborJoining` on the additive NJ example the result is a trifurcation: `(((a:2.0000,b:3.0000):3.0000,c:4.0000):2.0000,d:2.0000,e:1.0000);`.

## References

- Wikipedia: UPGMA — Neighbor joining (worked examples).
- [PhylogeneticAnalyzer.BuildTreeFromMatrix](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L174)

## See Also

- [build_phylogenetic_tree](build_phylogenetic_tree.md)
- [distance_matrix](distance_matrix.md)
