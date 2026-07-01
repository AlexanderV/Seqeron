# find_heterochromatin_regions

Identify heterochromatin regions by k-mer repeat content.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_heterochromatin_regions` |
| **Method ID** | `ChromosomeAnalyzer.FindHeterochromatinRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Slides a window (step = windowSize/2) and marks windows whose k-mer repeat content ≥ `minRepeatContent`.
Contiguous marked windows are merged into a region and classified by its midpoint position:
`Telomeric` near either end (< 5% or > 95%), `Centromeric` near the middle (45–55%), else
`Constitutive`.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L1290](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1290)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Chromosome sequence. |
| `windowSize` | integer | No | Window size in bp (default 100000, > 0). |
| `minRepeatContent` | number | No | Minimum repeat content fraction (default 0.5). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ start, end, type }` (type ∈ Telomeric/Centromeric/Constitutive). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be positive |

## Example

300 kb of uniform `A` → a single `Constitutive` region spanning `[0, 299999]`.

## References

- [ChromosomeAnalyzer.FindHeterochromatinRegions](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1290)
