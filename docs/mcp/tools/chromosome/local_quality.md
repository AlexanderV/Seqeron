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

Slides a window over each sequence (starts 0, step, 2·step, …; step defaults to windowSize/2; the last window is
truncated at the sequence end, so every base is covered) and reports per window: the GC fraction
(canonical windowed GC, (G+C)/(A+C+G+T) — N/IUPAC codes are not in the denominator, the QUAST rule), the N count
and the canonical linguistic complexity (Orlov & Potapov 2004 / Rosalind LING, m = 10). `position` is the 0-based
window start and `windowSize` the actual (possibly truncated) length.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Sequences `{ id, sequence }`. |
| `windowSize` | integer | No | Window size (default 1000, > 0). |
| `stepSize` | integer | No | Distance between window starts (default windowSize/2, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, position, windowSize, gcContent, nCount, complexity }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | Window size must be positive |
| 1002 | Step size must be positive |

## Example

`"GCGCGCGCGC"` (10 bp) with `windowSize = 1000` → one window `{ position: 0, windowSize: 10,
gcContent: 1.0, nCount: 0 }`.

## References

- [GenomeAssemblyAnalyzer.CalculateLocalQuality](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs)
- [find_suspicious_regions](find_suspicious_regions.md)
