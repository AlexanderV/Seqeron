# tree_depth

Tree height: the maximum number of edges from the root to any leaf.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `tree_depth` |
| **Method ID** | `PhylogeneticAnalyzer.GetTreeDepth` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses a Newick tree and returns its height: the number of edges on the longest downward path from the root to a leaf. A single-leaf tree has depth `0`. (An empty/null tree would be `-1` by convention, but the tool always parses a non-empty Newick string.)

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L861](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L861)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `newick` | string | Yes | Tree in Newick format (min length 1) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `depth` | integer | Tree height in edges; a single-leaf tree has depth 0 |

## Errors

| Code | Message |
|------|---------|
| 2001 | Newick string is empty |

## Example

**User Prompt:**
> What is the depth of `((A,B),(C,D))`?

**Tool Call:**
```json
{ "tool": "tree_depth", "arguments": { "newick": "((A,B),(C,D));" } }
```

**Response:**
```json
{ "depth": 2 }
```

The longest root→leaf path is root → internal → leaf = 2 edges. A caterpillar `(A,(B,(C,D)))` has depth 3; a two-leaf tree `(A,B)` has depth 1.

## References

- Tree (graph theory) / Tree (abstract data type): height definition.
- [PhylogeneticAnalyzer.GetTreeDepth](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L861)

## See Also

- [tree_length](tree_length.md)
- [tree_leaves](tree_leaves.md)
- [parse_newick](parse_newick.md)
