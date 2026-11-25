# find_suspicious_regions

Sequence-only screen for anomalous assembly regions (not a misassembly caller).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_suspicious_regions` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindSuspiciousRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

A composition screen — QUAST/GAGE misassemblies need an alignment to a reference, and no published threshold
defines "suspicious" GC or N content, so the cut-offs are parameters and no composite score is reported. Three
criteria, merged into inclusive 0-based regions (union of the flagged windows/intervals; overlapping or adjacent
ones are joined; reasons are listed in the fixed order `High N content; GC deviation; Low complexity`):

- **High N content** — a window (default 500, step windowSize/2) whose N fraction exceeds `maxNFraction`.
- **GC deviation** — a window whose GC fraction differs from the pooled assembly GC, Σ(G+C)/Σ(A+C+G+T), by more than
  `gcDeviation`. A truncated trailing window already covered by the previous window and a window without any
  A/C/G/T base are not tested.
- **Low complexity** — exact intervals of the symmetric DUST algorithm (Morgulis et al. 2006; lh3/sdust, NCBI
  dustmasker; window 64, DUST score > `dustThreshold`, 2.0 = level 20).

Each region reports the measured value of every criterion that fired: `maxGcDeviation`, `maxNFraction`,
`maxDustScore` (0 when the criterion did not fire).

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Sequences `{ id, sequence }`. |
| `gcDeviation` | number | No | Allowed \|window GC − pooled GC\| as a fraction (default 0.15). |
| `dustThreshold` | number | No | DUST score threshold (default 2.0 = dustmasker level 20). |
| `maxNFraction` | number | No | Maximum N fraction per window (default 0.1). |
| `windowSize` | integer | No | GC/N window size (default 500). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, start, end, reason, maxGcDeviation, maxNFraction, maxDustScore }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |

## Example

A 600-bp poly-A sequence → one region `[0, 599]`, reason `Low complexity`, `maxDustScore` 299 (= 598·597/2/597).

## References

- [GenomeAssemblyAnalyzer.FindSuspiciousRegions](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs)
- [local_quality](local_quality.md)
