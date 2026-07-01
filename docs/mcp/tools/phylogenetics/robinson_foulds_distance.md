# robinson_foulds_distance

Robinson-Foulds (symmetric clade-difference) distance between two rooted trees.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `robinson_foulds_distance` |
| **Method ID** | `PhylogeneticAnalyzer.RobinsonFouldsDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the Robinson-Foulds distance between two **rooted** trees over the same taxon set as the size of the symmetric difference of their clade (cluster) sets. Only non-trivial clades are compared; trivial clades (a single leaf, or the full taxon set) are excluded. The result is a non-negative even integer. Both trees are supplied as Newick strings and parsed internally.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L879](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L879)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `tree1` | string | Yes | First tree in Newick format |
| `tree2` | string | Yes | Second tree in Newick format |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `distance` | integer | Symmetric clade difference (non-negative even integer) |

## Errors

| Code | Message |
|------|---------|
| 2001 | Newick string is empty |

## Example

**User Prompt:**
> How different are `((A,B),(C,D))` and `((A,C),(B,D))`?

**Tool Call:**
```json
{ "tool": "robinson_foulds_distance", "arguments": { "tree1": "((A,B),(C,D));", "tree2": "((A,C),(B,D));" } }
```

**Response:**
```json
{ "distance": 4 }
```

The two trees share no non-trivial clade (`{A,B},{C,D}` vs `{A,C},{B,D}`), giving the maximum RF of `2(n-2) = 4` for `n = 4`. Identical trees return `0`; collapsing one internal edge changes RF by 1.

## References

- Robinson DF, Foulds LR (1981). *Comparison of phylogenetic trees.* Math. Biosci. 53:131–147.
- [PhylogeneticAnalyzer.RobinsonFouldsDistance](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L879)

## See Also

- [bootstrap_support](bootstrap_support.md)
- [parse_newick](parse_newick.md)
