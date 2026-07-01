# pairwise_distance

Evolutionary distance between two aligned sequences under a chosen substitution model.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Phylogenetics |
| **Tool Name** | `pairwise_distance` |
| **Method ID** | `PhylogeneticAnalyzer.CalculatePairwiseDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the distance between two equal-length aligned sequences:

- **Hamming** — raw count of differing sites.
- **PDistance** — proportion of differing sites (`differences / comparable sites`).
- **JukesCantor** (JC69) — `d = -3/4 · ln(1 - 4p/3)`; `+Infinity` when `p ≥ 0.75` (saturation).
- **Kimura2Parameter** (K80) — `d = -1/2 · ln((1 - 2S - V) · √(1 - 2V))`, distinguishing transitions (S) from transversions (V); `+Infinity` at saturation.

Gaps (`-`) and non-standard/ambiguous bases are skipped; case is ignored.

## Core Documentation Reference

- Source: [PhylogeneticAnalyzer.cs#L223](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L223)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seq1` | string | Yes | First aligned sequence |
| `seq2` | string | Yes | Second aligned sequence (same length as `seq1`) |
| `method` | string | No | `PDistance` \| `JukesCantor` \| `Kimura2Parameter` \| `Hamming` (default `JukesCantor`) |

## Returns

| Field | Type | Description |
|-------|------|-------------|
| `distance` | number | The evolutionary distance (may be `+Infinity` at JC69/K2P saturation) |

## Errors

| Code | Message |
|------|---------|
| 1001 | First sequence cannot be null |
| 1002 | Second sequence cannot be null |
| 1003 | Sequences must have the same length |
| 1004 | Unknown distance method |

## Example

**User Prompt:**
> Jukes-Cantor distance between `ACGTACGT` and `TCGTACGT`?

**Tool Call:**
```json
{ "tool": "pairwise_distance", "arguments": { "seq1": "ACGTACGT", "seq2": "TCGTACGT", "method": "JukesCantor" } }
```

**Response:**
```json
{ "distance": 0.13674 }
```

One difference in 8 sites gives `p = 1/8`; the JC69 correction `-3/4·ln(5/6) ≈ 0.13674`.

## References

- Jukes & Cantor (1969); Kimura (1980).
- [PhylogeneticAnalyzer.CalculatePairwiseDistance](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs#L223)

## See Also

- [distance_matrix](distance_matrix.md)
