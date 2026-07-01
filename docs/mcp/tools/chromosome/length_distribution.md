# length_distribution

Bucket sequence lengths into bins.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `length_distribution` |
| **Method ID** | `GenomeAssemblyAnalyzer.CalculateLengthDistribution` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assigns each length to the first `<bin` bucket it falls under (default upper bounds
`100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000`), or to the `>=<maxBin>` overflow
bucket when it exceeds the largest bin.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L1202](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1202)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `lengths` | integer[] | Yes | Sequence lengths. |
| `bins` | integer[] | No | Bin upper bounds; empty/omitted → defaults. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distribution` | object | Map of bucket label (`"<100"`, …, `">=1000000"`) → count. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Lengths cannot be null |

## Example

`[50, 150, 600]` → `{ "<100": 1, "<500": 1, "<1000": 1, … }`.

## References

- [GenomeAssemblyAnalyzer.CalculateLengthDistribution](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1202)
