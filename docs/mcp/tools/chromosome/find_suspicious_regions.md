# find_suspicious_regions

Flag potentially misassembled regions.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_suspicious_regions` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindSuspiciousRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Slides a window over each sequence and flags windows whose GC deviates from the global GC by more than
`gcDeviation`, or whose linguistic complexity is below `minComplexity` (with low N content). Adjacent
flagged windows merge into a region carrying the combined reason and a severity score.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L1086](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1086)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Sequences `{ id, sequence }`. |
| `gcDeviation` | number | No | Allowed GC deviation from global GC (default 0.15). |
| `minComplexity` | number | No | Minimum linguistic complexity (default 0.3). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, start, end, reason, score }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |

## Example

A 2 kb poly-A run embedded in random sequence is flagged with a reason that includes `Low complexity`.

## References

- [GenomeAssemblyAnalyzer.FindSuspiciousRegions](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1086)
- [local_quality](local_quality.md)
