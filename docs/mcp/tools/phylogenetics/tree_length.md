# tree_length

Sum of all branch lengths in a tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `tree_length` |
| **Method ID** | `PhylogeneticAnalyzer.CalculateTreeLength` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses a Newick tree and returns its total length: the sum of the branch lengths of every node in the tree (including the root node's own branch length). This is the quantity minimized by the minimum-evolution criterion. Negative branch lengths (possible from Neighbor-Joining output) are summed as-is. Undefined branch lengths default to 0.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L830](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L830)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `newick` | string | Yes | Tree in Newick format (min length 1) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `length` | number | Sum of all branch lengths (root edge included) |

## Errors

| Code | Message |
|------|---------|
| 2001 | Newick string is empty |

## Example

**User Prompt:**
> Total branch length of `((A:1,B:1):1,(C:1,D:1):1)`?

**Tool Call:**
```json
{ "tool": "tree_length", "arguments": { "newick": "((A:1,B:1):1,(C:1,D:1):1);" } }
```

**Response:**
```json
{ "length": 6 }
```

Four leaf edges (1 each) plus two internal edges (1 each) sum to 6. A caterpillar `(A:1,(B:1,(C:1,D:1):0.5):0.5)` sums to 5.

## References

- DendroPy `Tree.length()`; Biopython `Tree.total_branch_length`; Rzhetsky & Nei (1992).
- [PhylogeneticAnalyzer.CalculateTreeLength](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L830)

## See Also

- [tree_depth](tree_depth.md)
- [tree_leaves](tree_leaves.md)
- [patristic_distance](patristic_distance.md)
