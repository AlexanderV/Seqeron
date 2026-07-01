# distance_matrix

Compute the symmetric pairwise distance matrix for a list of aligned sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `distance_matrix` |
| **Method ID** | `PhylogeneticAnalyzer.CalculateDistanceMatrix` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the symmetric `n × n` pairwise distance matrix for aligned sequences under the chosen substitution model (`PDistance` | `JukesCantor` | `Kimura2Parameter` | `Hamming`). The diagonal is zero and `d(i,j) = d(j,i)`. Gaps (`-`) and non-standard/ambiguous bases are skipped per column; case is ignored. `JukesCantor` and `Kimura2Parameter` may return `+Infinity` at saturation.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L199](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L199)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedSequences` | array of string | Yes | Aligned sequences, equal length (≥1) |
| `method` | string | No | `PDistance` \| `JukesCantor` \| `Kimura2Parameter` \| `Hamming` (default `JukesCantor`) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `matrix` | number[][] | Symmetric `n × n` distance matrix; diagonal is zero |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one aligned sequence is required |
| 1002 | Sequences must have the same length |
| 1003 | Unknown distance method |

## Example

**User Prompt:**
> What are the Hamming distances between AAAA, AAAC and CCCC?

**Tool Call:**
```json
{
  "tool": "distance_matrix",
  "arguments": { "alignedSequences": ["AAAA", "AAAC", "CCCC"], "method": "Hamming" }
}
```

**Response:**
```json
{ "matrix": [[0, 1, 4], [1, 0, 3], [4, 3, 0]] }
```

`AAAA↔AAAC` differ at one site (1), `AAAA↔CCCC` at four (4), `AAAC↔CCCC` at three (3).

## References

- Jukes & Cantor (1969); Kimura (1980).
- [PhylogeneticAnalyzer.CalculateDistanceMatrix](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L199)

## See Also

- [pairwise_distance](pairwise_distance.md)
- [build_tree_from_matrix](build_tree_from_matrix.md)
