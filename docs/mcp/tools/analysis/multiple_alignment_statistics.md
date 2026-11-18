# multiple_alignment_statistics

Statistics of a gapped multiple alignment: per-sequence comparison with a reference row or the EMBOSS cons consensus (EMBOSS infoalign: seqLength, alignLength, gap runs, gap length, identities, similarities, differences, % change; verified against infoalign 6.6.0), per-column residue counts / gap fraction / most common residue / ClustalW conservation symbol ('*' identical, ':' strong, '.' weak group; ClustalW 2.1; DNA columns only '*'), gap percent, fully conserved and gap-free columns, and the pairwise percent-identity matrix (distmat uncorrected, identity over columns where neither row has a gap; gapWeight 1 = Biopython DistanceCalculator('identity')). Rows are upper-cased, '.'/'~' read as '-'.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `multiple_alignment_statistics` |
| **Method ID** | `MultipleAlignmentStatistics.Calculate` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Statistics of a gapped multiple alignment: per-sequence comparison with a reference row or the EMBOSS cons consensus (EMBOSS infoalign: seqLength, alignLength, gap runs, gap length, identities, similarities, differences, % change; verified against infoalign 6.6.0), per-column residue counts / gap fraction / most common residue / ClustalW conservation symbol ('*' identical, ':' strong, '.' weak group; ClustalW 2.1; DNA columns only '*'), gap percent, fully conserved and gap-free columns, and the pairwise percent-identity matrix (distmat uncorrected, identity over columns where neither row has a gap; gapWeight 1 = Biopython DistanceCalculator('identity')). Rows are upper-cased, '.'/'~' read as '-'.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MultipleAlignmentStatistics.cs#L110](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MultipleAlignmentStatistics.cs#L110)
- Algorithm: [Alignment_Statistics.md §5.6](../../../algorithms/Alignment/Alignment_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | At least two aligned sequences of equal length (DNA or protein, gaps '-'). |
| `residueType` | string | No | 'auto' (type of the first row as EMBOSS), 'nucleotide' (EDNAFULL) or 'protein' (EBLOSUM62; also enables ':'/'.' conservation). (default: auto) |
| `referenceSequence` | integer | No | infoalign -refseq: 0 = the EMBOSS cons consensus (default), k = the k-th sequence (1-based). (default: 0) |
| `gapWeight` | number | No | Weight of gap columns in the identity matrix denominator L - g + g*gapWeight (distmat -gapweight; default 0). (default: 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceCount` | integer | Rows. |
| `alignmentLength` | integer | Columns. |
| `fullyConservedColumns` | integer | Columns with ClustalW symbol '*'. |
| `gapFreeColumns` | integer | Columns without a gap. |
| `gapPercent` | number | 100 * gap symbols / (rows * columns). |
| `meanPairwiseIdentity` | number | Mean off-diagonal percent identity (null when undefined). |
| `conservationLine` | string | ClustalW conservation line ('*', ':', '.', ' '). |
| `pairwiseIdentity` | array<array<number>> | Percent-identity matrix (null = undefined). |
| `sequences` | array<object> | Per-row infoalign statistics: index, sequenceLength, alignmentLength, gapRuns, gapLength, identities, similarities, differences, percentChange. |
| `columns` | array<object> | Per-column: index, counts, gapCount, gapFraction, mostCommonResidue, mostCommonCount, conservation. |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least two sequences are required |
| 1002 | residueType must be 'nucleotide', 'protein' or 'auto' |
| — | referenceSequence outside 0..n (ArgumentOutOfRangeException); unequal row lengths (ArgumentException) |

## Examples

### Example 1: EMBOSS infoalign / distmat reference values

**Tool Call:**
```json
{
  "tool": "multiple_alignment_statistics",
  "arguments": {
    "sequences": [
      "ACGT-ACGTAC",
      "ACGTTACG-AC",
      "AGGT-ACGTAC",
      "A-GTTTCGTAA"
    ],
    "residueType": "nucleotide",
    "referenceSequence": 1
  }
}
```

**Result (checked against EMBOSS 6.6.0):** per-sequence row 2 (infoalign -refseq 1): sequenceLength 10, alignmentLength 11, gapRuns 1, gapLength 1, identities 9, similarities 0, differences 1, percentChange 18.181818; pairwiseIdentity[0][2] = 90 (distmat 10.00).
