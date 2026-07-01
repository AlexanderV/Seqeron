# assembly_statistics

Compute comprehensive assembly statistics.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `assembly_statistics` |
| **Method ID** | `GenomeAssemblyAnalyzer.CalculateStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Reports total sequence count and length (with and without N gaps), N50/L50 and N90/L90 contiguity
metrics (Miller, Koren & Sutton 2010), largest/smallest contig, mean and median length, GC content
(over non-N bases), and gap counts/lengths.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L121](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L121)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Assembled sequences `{ id, sequence }`. Empty → all-zero stats. |

## Output Schema

Fields: `totalSequences`, `totalLength`, `totalLengthNoGaps`, `n50`, `l50`, `n90`, `l90`,
`largestContig`, `smallestContig`, `meanLength`, `medianLength`, `gcContent`, `totalGaps`,
`totalGapLength`, `gapPercentage`.

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |

## Example

`s1 = "ACGTACGTGG"` (10 bp, 6 GC), `s2 = "ACGTAC"` (6 bp, 3 GC):

```json
{
  "totalSequences": 2,
  "totalLength": 16,
  "largestContig": 10,
  "smallestContig": 6,
  "meanLength": 8.0,
  "gcContent": 0.5625,
  "n50": 10,
  "l50": 1
}
```

## References

- [GenomeAssemblyAnalyzer.CalculateStatistics](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L121)
