# parse_newick

Parse a Newick-format tree string and report a structural summary.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `parse_newick` |
| **Method ID** | `PhylogeneticAnalyzer.ParseNewick` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses a Newick tree and returns: a canonical re-serialization (round-trip), the leaf/taxon list (in pre-order), the leaf count, the tree depth (height in edges), and the total branch length. A trailing `;` is stripped and an optional root branch length is supported. Malformed input (unbalanced parentheses, trailing garbage) is rejected with a format error; an empty/whitespace string is rejected as an argument error.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L657](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L657)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `newick` | string | Yes | Tree in Newick format (min length 1) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `newick` | string | Canonical re-serialization of the parsed tree |
| `taxa` | array of string | Leaf names in pre-order |
| `leafCount` | integer | Number of leaves |
| `depth` | integer | Tree height in edges (single leaf = 0) |
| `totalLength` | number | Sum of all branch lengths |

## Errors

| Code | Message |
|------|---------|
| 2001 | Newick string is empty |
| 2002 | Malformed Newick string (unbalanced parentheses or trailing garbage) |

## Example

**User Prompt:**
> Summarize the tree `((A:1,B:1):1,(C:1,D:1):1)`.

**Tool Call:**
```json
{ "tool": "parse_newick", "arguments": { "newick": "((A:1,B:1):1,(C:1,D:1):1);" } }
```

**Response (excerpt):**
```json
{ "taxa": ["A", "B", "C", "D"], "leafCount": 4, "depth": 2, "totalLength": 6 }
```

Four leaf edges (1 each) plus two internal edges (1 each) sum to a total length of 6; the deepest leaf is 2 edges from the root.

## References

- Olsen (1990); Wikipedia: Newick format.
- [PhylogeneticAnalyzer.ParseNewick](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L657)

## See Also

- [to_newick](to_newick.md)
- [tree_leaves](tree_leaves.md)
- [tree_depth](tree_depth.md)
- [tree_length](tree_length.md)
