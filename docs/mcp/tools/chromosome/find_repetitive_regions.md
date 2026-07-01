# find_repetitive_regions

Identify repetitive regions using k-mer copy-number frequency.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_repetitive_regions` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindRepetitiveRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Counts k-mer occurrences across the sequences and flags spans where k-mers recur at least `minCopies`
times, merging adjacent hits within `windowSize` into a region with a representative copy count.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L674](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L674)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Sequences `{ id, sequence }`. |
| `kmerSize` | integer | No | K-mer size (default 15, > 0). |
| `minCopies` | integer | No | Minimum k-mer copies (default 3, > 0). |
| `windowSize` | integer | No | Merge window (default 100, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, start, end, copies }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | K-mer size must be positive |
| 1003 | Minimum copies must be positive |
| 1004 | Window size must be positive |

## Example

`"ACGTACGT"` × 100 (800 bp) with `kmerSize = 8` → one region `[0, 789]` with `copies = 199`.

## References

- [GenomeAssemblyAnalyzer.FindRepetitiveRegions](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L674)
