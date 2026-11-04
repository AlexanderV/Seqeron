# multiple_align

Center-star multiple sequence alignment (Gusfield 1993) with exact pairwise alignments: the center is the sequence with the highest total optimal Needleman-Wunsch score to the others, every other sequence is aligned to it and the alignments are merged column-wise (once a gap, always a gap). Linear gap model (gapOpen ignored); for two sequences the result is an optimal pairwise global alignment. Returns rows, consensus and column sum-of-pairs score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align` |
| **Method ID** | `SequenceAligner.MultipleAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Center-star multiple sequence alignment (Gusfield 1993) with exact pairwise alignments: the center is the sequence with the highest total optimal Needleman-Wunsch score to the others, every other sequence is aligned to it and the alignments are merged column-wise (once a gap, always a gap). Linear gap model (gapOpen ignored); for two sequences the result is an optimal pairwise global alignment. Returns rows, consensus and column sum-of-pairs score.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1203](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1203)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences to align (at least one). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty (ignored: linear-gap model; kept for schema compatibility). (default: -2) |
| `gapExtend` | integer | No | Gap-extend penalty. (default: -1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA. |
| `consensus` | string | Majority-vote consensus of the aligned rows. |
| `totalScore` | integer | Sum-of-pairs score of the MSA. |

## Errors

| Code | Message |
|------|---------|
| 1010 | Invalid DNA in sequences |

## Examples

### Example 1: Three identical

**Tool Call:**
```json
{
  "tool": "multiple_align",
  "arguments": {
    "sequences": [
      "ACGTACGT",
      "ACGTACGT",
      "ACGTACGT"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTACGT",
    "ACGTACGT",
    "ACGTACGT"
  ],
  "consensus": "ACGTACGT",
  "totalScore": 24
}
```

### Example 2: Empty set

**Tool Call:**
```json
{
  "tool": "multiple_align",
  "arguments": {
    "sequences": []
  }
}
```

**Response:**
```json
{
  "alignedSequences": [],
  "consensus": "",
  "totalScore": 0
}
```

## Worked Example

Three identical 8-mers align gaplessly; SP = C(3,2)=3 pairs x 8 matched columns x (+1) = 24. For two sequences the result is the optimal Needleman-Wunsch alignment (GGAAGCCAATAATCGTGCACGAAGGTGACTGCAGGAGAGTACGAAGGT vs AGACAGCCAAAATCTTCAGAGTTCGACTGGAGGTGACTACGAAGAT scores 24, Biopython 1.88).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1203](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1203)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
