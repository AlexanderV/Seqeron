# gap_distribution

Summarize a list of gaps.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `gap_distribution` |
| **Method ID** | `GenomeAssemblyAnalyzer.AnalyzeGapDistribution` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Reports the gap count, mean/median/max gap length and per-type counts. The median is the sample
median (`numpy.median`: mean of the two central values for an even count). An empty list returns all zeros.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L388](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L388)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gaps` | array | Yes | Gaps (typically from `find_gaps`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `count` | integer | Number of gaps. |
| `meanLength` | number | Mean gap length. |
| `medianLength` | number | Sample median gap length. |
| `maxLength` | integer | Maximum gap length. |
| `typeCounts` | object | Map of gap type → count. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Gaps cannot be null |

## Example

Gaps of length 5 (`Short`) and 100 (`Long`) → `count 2, mean 52.5, median 52.5, max 100`.

## References

- [GenomeAssemblyAnalyzer.AnalyzeGapDistribution](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L388)
- [find_gaps](find_gaps.md)
