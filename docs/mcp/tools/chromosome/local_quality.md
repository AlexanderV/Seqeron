# local_quality

Compute per-window local quality metrics.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `local_quality` |
| **Method ID** | `GenomeAssemblyAnalyzer.CalculateLocalQuality` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Slides a window (step = windowSize/2) over each sequence and reports the GC fraction, N count and
linguistic complexity (distinct 4-mers / possible 4-mers) per window.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L1029](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1029)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Sequences `{ id, sequence }`. |
| `windowSize` | integer | No | Window size (default 1000, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, position, windowSize, gcContent, nCount, complexity }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | Window size must be positive |

## Example

`"GCGCGCGCGC"` (10 bp) with `windowSize = 1000` → one window `{ position: 0, windowSize: 10,
gcContent: 1.0, nCount: 0 }`.

## References

- [GenomeAssemblyAnalyzer.CalculateLocalQuality](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1029)
- [find_suspicious_regions](find_suspicious_regions.md)
