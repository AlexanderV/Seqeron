# find_gaps

Find gaps (N-runs) in assembled sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_gaps` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindGaps` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans each sequence for maximal runs of `N`/`n`. A run of length ≥ `minGapLength` is reported as a
gap with inclusive `[start, end]`, its length, and a length class (`< 10` Short, `< 100` Medium,
`< 1000` Long, else Scaffold).

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L341](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L341)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Assembled sequences `{ id, sequence }`. |
| `minGapLength` | integer | No | Minimum N-run length to report (default 1, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, start, end, length, gapType }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | Minimum gap length must be positive |

## Example

`"AAANNNNNGGG"` → one gap `[3, 7]`, length 5, type `Short`.

## References

- [GenomeAssemblyAnalyzer.FindGaps](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L341)
- [gap_distribution](gap_distribution.md)
