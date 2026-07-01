# to_newick

Serialize a phylogenetic tree to canonical Newick format.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `to_newick` |
| **Method ID** | `PhylogeneticAnalyzer.ToNewick` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Accepts a Newick string, parses it internally, and re-emits canonical Newick. Branch lengths are formatted with `F4` (four decimals, invariant culture) when `includeBranchLengths` is true, and omitted otherwise. Internal node names are emitted only if they are valid unquoted Newick labels (Olsen 1990); auto-generated names containing Newick metacharacters are dropped. Round-trips with `parse_newick` for trees whose labels are valid unquoted labels.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L577](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L577)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `newick` | string | Yes | Input tree in Newick format (min length 1) |
| `includeBranchLengths` | boolean | No | Whether to include branch lengths (default `true`) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `newick` | string | Canonical Newick serialization |

## Errors

| Code | Message |
|------|---------|
| 2001 | Newick string is empty |
| 2002 | Malformed Newick string |

## Example

**User Prompt:**
> Normalize `(A:0.1,B:0.2,C:0.3)` to canonical Newick.

**Tool Call:**
```json
{ "tool": "to_newick", "arguments": { "newick": "(A:0.1,B:0.2,C:0.3);", "includeBranchLengths": true } }
```

**Response:**
```json
{ "newick": "(A:0.1000,B:0.2000,C:0.3000);" }
```

Branch lengths are re-emitted with four decimals. With `includeBranchLengths = false` the same input yields `(A,B,C);`.

## References

- Olsen (1990); Wikipedia: Newick format.
- [PhylogeneticAnalyzer.ToNewick](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L577)

## See Also

- [parse_newick](parse_newick.md)
- [build_phylogenetic_tree](build_phylogenetic_tree.md)
