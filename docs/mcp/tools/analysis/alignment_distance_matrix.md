# alignment_distance_matrix

Pairwise distance matrix of a multiple alignment as EMBOSS 6.6.0 distmat: matches = columns where both rows carry the same non-gap symbol, g = columns where either row has a gap; 'uncorrected' D = 1 - matches/(L - g + g*gapWeight); 'jukes-cantor' -b*ln(1 - D/b) with b = 3/4 (nucleotide) or 19/20 (protein). Fractions, not percent; entries with a zero denominator or a saturated logarithm are null.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `alignment_distance_matrix` |
| **Method ID** | `MultipleAlignmentStatistics.CalculateDistanceMatrix` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Pairwise distance matrix of a multiple alignment as EMBOSS 6.6.0 distmat: matches = columns where both rows carry the same non-gap symbol, g = columns where either row has a gap; 'uncorrected' D = 1 - matches/(L - g + g*gapWeight); 'jukes-cantor' -b*ln(1 - D/b) with b = 3/4 (nucleotide) or 19/20 (protein). Fractions, not percent; entries with a zero denominator or a saturated logarithm are null.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MultipleAlignmentStatistics.cs#L220](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MultipleAlignmentStatistics.cs#L220)
- Algorithm: [Alignment_Statistics.md §5.6](../../../algorithms/Alignment/Alignment_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | At least two aligned sequences of equal length. |
| `method` | string | No | 'uncorrected' (default) or 'jukes-cantor'. (default: uncorrected) |
| `residueType` | string | No | 'auto' (default), 'nucleotide' or 'protein' (selects the Jukes-Cantor constant b). (default: auto) |
| `gapWeight` | number | No | distmat -gapweight (default 0). (default: 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `method` | string | Method used. |
| `matrix` | array<array<number>> | Symmetric distance matrix (fractions, zero diagonal; null = undefined or saturated). |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least two sequences are required |
| 1002 | method must be 'uncorrected' or 'jukes-cantor' |

## Examples

### Example 1: EMBOSS infoalign / distmat reference values

**Tool Call:**
```json
{
  "tool": "alignment_distance_matrix",
  "arguments": {
    "sequences": [
      "ACGT-ACGTAC",
      "ACGTTACG-AC",
      "AGGT-ACGTAC",
      "A-GTTTCGTAA"
    ]
  }
}
```

**Result (checked against EMBOSS 6.6.0):** matrix[0][2] = 0.1, matrix[1][2] = 0.1111 (EMBOSS distmat 6.6.0: 10.00, 11.11 percent).
