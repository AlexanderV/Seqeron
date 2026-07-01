# patristic_distance

Sum of branch lengths along the unique path between two taxa in a rooted tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `patristic_distance` |
| **Method ID** | `PhylogeneticAnalyzer.PatristicDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the patristic (tree-path) distance between two taxa: the sum of the branch lengths on the unique path connecting them, routed through their MRCA. Returns `NaN` when the MRCA cannot be found (one or both taxa missing). Identical taxa return `0`. Multifurcations (polytomies) are traversed over all children.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L1125](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L1125)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `newick` | string | Yes | Tree in Newick format (min length 1) |
| `taxon1` | string | Yes | First taxon name |
| `taxon2` | string | Yes | Second taxon name |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `distance` | number | Patristic distance; `NaN` if a taxon is missing; `0` for identical taxa |

## Errors

| Code | Message |
|------|---------|
| 2001 | Newick string is empty |

## Example

**User Prompt:**
> Patristic distance between A and B in `(A:1,B:2)`.

**Tool Call:**
```json
{ "tool": "patristic_distance", "arguments": { "newick": "(A:1,B:2);", "taxon1": "A", "taxon2": "B" } }
```

**Response:**
```json
{ "distance": 3 }
```

The path A→root→B sums the two leaf branch lengths (1 + 2 = 3). A missing taxon returns `NaN`; `taxon1 == taxon2` returns `0`.

## References

- [PhylogeneticAnalyzer.PatristicDistance](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L1125)

## See Also

- [mrca](mrca.md)
- [tree_length](tree_length.md)
