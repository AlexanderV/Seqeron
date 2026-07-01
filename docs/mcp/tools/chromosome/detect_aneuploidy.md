# detect_aneuploidy

Detect copy-number states from binned read-depth data.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `detect_aneuploidy` |
| **Method ID** | `ChromosomeAnalyzer.DetectAneuploidy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Groups depth samples by chromosome, bins them by `position / binSize`, and for each bin computes
`logRatio = log2(meanDepth / medianDepth)` and `copyNumber = round(2^logRatio × 2)` (clamped to
`[0, 10]`), with a confidence based on the deviation from the nearest integer copy number. A bin at
exactly the median depth is diploid with confidence 1.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L1537](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1537)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `depthData` | array | Yes | `{ chromosome, position, depth }` samples. |
| `medianDepth` | number | Yes | Genome-wide median depth (> 0). |
| `binSize` | integer | No | Bin size in bp (default 1,000,000, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ chromosome, start, end, copyNumber, logRatio, confidence }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Depth data cannot be null |
| 1002 | Median depth must be positive |
| 1003 | Bin size must be positive |

## Example

Five positions in bin 0 all at depth 30 with `medianDepth = 30` →
`{ start: 0, end: 999999, copyNumber: 2, logRatio: 0.0, confidence: 1.0 }`.

## References

- [ChromosomeAnalyzer.DetectAneuploidy](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1537)
- [identify_whole_chromosome_aneuploidy](identify_whole_chromosome_aneuploidy.md)
