# local_align_zscore

Shuffle Z-score of the best local alignment score of two sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `local_align_zscore` |
| **Method ID** | `AlignmentSignificance.ShuffleZScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The shuffle Z-score of the best local alignment score of two sequences (the idea of EMBOSS prss and FASTA ssearch): sequence2 is shuffled `shuffles` times (a seeded Fisher-Yates shuffle, reproducible), every shuffle is aligned locally to sequence1 with the same linear-gap scoring as local_align, and z = (score - mean) / standard deviation of the shuffled scores. A large z (above about 5-8) means the score is far above what the composition of the sequences gives by chance. Use local_align_affine / local_align / water_align with significance=true for the analytic bit score and E-value of the Karlin-Altschul statistics.

## Core Documentation Reference

- Source: [AlignmentSignificance.Gapped.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/AlignmentSignificance.Gapped.cs)
- Algorithm: [Alignment_Significance_Statistics.md](../../../algorithms/Alignment/Alignment_Significance_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence (the one that is shuffled). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapExtend` | number | No | Gap score charged per gap position (<= 0; the linear model of local_align). (default: -1) |
| `shuffles` | integer | No | Number of shuffles (>= 2; default 1000). (default: 1000) |
| `seed` | integer | No | Seed of the shuffles (default 1). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | integer | The best local score of the pair |
| `mean` | number | Mean of the best local scores of the shuffles of sequence2 |
| `standardDeviation` | number | Sample standard deviation of the shuffled scores |
| `zScore` | number | (score - mean) / standardDeviation (0 when every shuffled score is equal) |
| `shuffles` | integer | Number of shuffles |
| `seed` | integer | Seed of the shuffles |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required (null); an empty string is accepted |
| 1002 | shuffles must be at least 2 / gapExtend is not an integer |

## Examples

### Example 1: Shared 40-mer

**Input:** `{ "sequence1": "TTTTACGTTGCAAGGCTTAACCGGTTAACGTACGGTACCATGCAGGGG", "sequence2": "CCCCAAAAACGTTGCAAGGCTTAACCGGTTAACGTACGGTACCATGCATTTTCCCC", "shuffles": 200, "seed": 7 }`

**Response:** a z-score above 5 (the value is reproducible for the seed).

## See Also

- [local_align](local_align.md)
- [local_align_affine](local_align_affine.md)
