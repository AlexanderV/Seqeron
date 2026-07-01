# nx_curve

Compute Nx/Lx for multiple thresholds.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `nx_curve` |
| **Method ID** | `GenomeAssemblyAnalyzer.CalculateNxCurve` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes `Nx`/`Lx` (via `CalculateNx`) for each threshold, sorted ascending. When no thresholds are
supplied, the default deciles `10, 20, …, 90` are used.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L302](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L302)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `lengths` | integer[] | Yes | Sequence lengths (any order). |
| `thresholds` | integer[] | No | Thresholds; empty/omitted → 10..90 step 10. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ threshold, nx, lx, cumulativeLength }` per threshold, sorted ascending. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Lengths cannot be null |

## Example

`[100,90,80,70,60,50,40,30,20,10]` with default thresholds returns 9 entries (10..90); the N50 entry
has `nx = 70, lx = 4`.

## References

- [GenomeAssemblyAnalyzer.CalculateNxCurve](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L302)
- [nx_statistics](nx_statistics.md)
