# find_tandem_repeats

Identify tandem repeats.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_tandem_repeats` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindTandemRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans for tandem arrays with a repeat unit of `minUnitLength`–`maxUnitLength` bp occurring at least
`minCopies` times, reporting the unit, copy number and purity (fraction of matching bases) per
occurrence.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L753](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L753)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Sequences `{ id, sequence }`. |
| `minUnitLength` | integer | No | Minimum unit length (default 2, > 0). |
| `maxUnitLength` | integer | No | Maximum unit length (default 50, ≥ min). |
| `minCopies` | integer | No | Minimum copies (default 3, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, start, end, unit, copies, purity }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | Minimum unit length must be positive |
| 1003 | Maximum unit length must be >= minimum unit length |
| 1004 | Minimum copies must be positive |

## Example

`"CAG"` × 30 (90 bp) → `{ unit: "CAG", copies: 30, start: 0, end: 89, purity: 1.0 }`.

## References

- [GenomeAssemblyAnalyzer.FindTandemRepeats](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L753)
