# predict_g_bands

Predict a cytogenetic G-band pattern from sequence GC content.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `predict_g_bands` |
| **Method ID** | `ChromosomeAnalyzer.PredictGBands` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits the sequence into `bandSize` windows and stains each by GC: `< darkBandGcThreshold` → `gpos100`
(dark), `< lightBandGcThreshold` → `gpos50`, else `gneg` (light). Bands are named `{chr}{arm}{n}`,
switching from the `p` arm to the `q` arm at the sequence midpoint. This is a simplified model, not a
substitute for reference cytogenetic ideograms.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L1230](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1230)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chromosomeName` | string | Yes | Chromosome name. |
| `sequence` | string | Yes | Chromosome sequence. |
| `bandSize` | integer | No | Band size in bp (default 5,000,000, > 0). |
| `darkBandGcThreshold` | number | No | GC cutoff for dark band (default 0.37). |
| `lightBandGcThreshold` | number | No | GC cutoff for light band (default 0.45). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ chromosome, start, end, name, stain, gcContent, geneDensity }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Chromosome name cannot be null or empty |
| 1003 | Band size must be positive |

## Example

`"ATATATATATGCGCGCGCGC"` with `bandSize = 10` → band `chr1p1` (`gpos100`, GC 0) and band `chr1q1`
(`gneg`, GC 1).

## References

- [ChromosomeAnalyzer.PredictGBands](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1230)
