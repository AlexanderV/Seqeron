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
metrics (Miller, Koren & Sutton 2010; QUAST `N50.py`), largest/smallest contig, mean and sample
median length (mean of the two central lengths for an even count), GC content as QUAST `GC_content`
— pooled (G+C)/(A+C+G+T), with N, IUPAC ambiguity codes and gap symbols in neither count — and gap
(maximal N-run) counts/lengths. As in QUAST `correct_seq`, every symbol other than A/C/G/T (IUPAC codes, U, ...)
is an N, i.e. a gap symbol, unless `nonAcgtAsGap = false` (then only `N`/`n` is a gap).

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L121](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L121)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Assembled sequences `{ id, sequence }`. Empty → all-zero stats. |
| `nonAcgtAsGap` | boolean | No | QUAST `correct_seq`: treat every non-ACGT symbol (IUPAC codes, U, ...) as N, i.e. a gap symbol (default true); false = only N/n. |

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
