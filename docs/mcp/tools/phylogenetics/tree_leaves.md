# tree_leaves

Enumerate the leaf (taxon) nodes of a tree with their branch lengths.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `tree_leaves` |
| **Method ID** | `PhylogeneticAnalyzer.GetLeaves` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses a Newick tree and returns each leaf (terminal/taxon) node in left-to-right pre-order, with its name and its branch length to the parent. A leaf is a node with no children.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L801](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L801)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `newick` | string | Yes | Tree in Newick format (min length 1) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `leaves` | array | Leaf nodes in pre-order |
| `leaves[].name` | string | Leaf (taxon) name |
| `leaves[].branchLength` | number | Branch length to the parent |

## Errors

| Code | Message |
|------|---------|
| 2001 | Newick string is empty |

## Example

**User Prompt:**
> List the leaves and branch lengths of `((A:0.1,B:0.2,C:0.3):0.0,(D:0.4,E:0.5):0.6)`.

**Tool Call:**
```json
{ "tool": "tree_leaves", "arguments": { "newick": "((A:0.1,B:0.2,C:0.3):0.0,(D:0.4,E:0.5):0.6);" } }
```

**Response:**
```json
{
  "leaves": [
    { "name": "A", "branchLength": 0.1 },
    { "name": "B", "branchLength": 0.2 },
    { "name": "C", "branchLength": 0.3 },
    { "name": "D", "branchLength": 0.4 },
    { "name": "E", "branchLength": 0.5 }
  ]
}
```

Leaves are returned in left-to-right pre-order with their per-leaf branch lengths.

## References

- Biopython `Tree.get_terminals`; Tree (graph theory): leaf definition.
- [PhylogeneticAnalyzer.GetLeaves](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L801)

## See Also

- [tree_depth](tree_depth.md)
- [tree_length](tree_length.md)
- [parse_newick](parse_newick.md)
