# mrca

Most Recent Common Ancestor (MRCA) of two taxa in a rooted tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `mrca` |
| **Method ID** | `PhylogeneticAnalyzer.FindMRCA` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses a Newick tree and returns the deepest node whose subtree contains both taxa. Returns `found = false` (with empty `name`/`taxa`) when either taxon is missing from the tree. Self-MRCA (`taxon1 == taxon2`) is the leaf itself when that taxon exists.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L1076](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L1076)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `newick` | string | Yes | Tree in Newick format (min length 1) |
| `taxon1` | string | Yes | First taxon name |
| `taxon2` | string | Yes | Second taxon name |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `found` | boolean | True when both taxa are present |
| `name` | string | Name of the MRCA node (empty for unnamed internal nodes) |
| `taxa` | array of string | Taxa contained in the MRCA subtree |

## Errors

| Code | Message |
|------|---------|
| 2001 | Newick string is empty |

## Example

**User Prompt:**
> In `((A,B),(C,D))`, what is the MRCA of A and C?

**Tool Call:**
```json
{ "tool": "mrca", "arguments": { "newick": "((A,B),(C,D));", "taxon1": "A", "taxon2": "C" } }
```

**Response:**
```json
{ "found": true, "name": "", "taxa": ["A", "B", "C", "D"] }
```

A and C sit in different clades, so their MRCA is the root (all four taxa). Siblings A and B instead return `taxa = ["A", "B"]`; a missing taxon returns `found = false`.

## References

- [PhylogeneticAnalyzer.FindMRCA](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L1076)

## See Also

- [patristic_distance](patristic_distance.md)
- [parse_newick](parse_newick.md)
