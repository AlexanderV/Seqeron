# multiple_align

Center-star multiple sequence alignment (Gusfield 1993) with exact pairwise alignments: the center is the sequence with the highest total optimal pairwise global score to the others, every other sequence is aligned to it with the canonical global aligner (Needleman–Wunsch, or Gotoh affine with gapModel='affine') and the pairwise alignments are merged column-wise (once a gap, always a gap). Default linear gap model (every gap position costs gapExtend; gapOpen is ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score. For two sequences the result is an optimal pairwise global alignment under either model. Returns the aligned rows, majority consensus and sum-of-pairs score (gap-gap columns score 0).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align` |
| **Method ID** | `SequenceAligner.MultipleAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Center-star multiple sequence alignment (Gusfield 1993) with exact pairwise alignments: the center is the sequence with the highest total optimal pairwise global score to the others, every other sequence is aligned to it with the canonical global aligner (Needleman–Wunsch, or Gotoh affine with gapModel='affine') and the pairwise alignments are merged column-wise (once a gap, always a gap). Default linear gap model (every gap position costs gapExtend; gapOpen is ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score. For two sequences the result is an optimal pairwise global alignment under either model. Returns the aligned rows, majority consensus and sum-of-pairs score (gap-gap columns score 0).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1758](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1758)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences to align (at least one). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty. Used only with gapModel='affine' (a gap of length k then costs gapOpen + k * gapExtend, Gotoh); ignored by the default linear model. (default: -2) |
| `gapExtend` | integer | No | Gap penalty per gap position (linear model: every gap position; affine model: the per-position part of gapOpen + k * gapExtend). (default: -1) |
| `gapModel` | string | No | 'linear' (default): every gap position costs gapExtend and gapOpen is ignored; totalScore is the column sum of pairs (gap-gap columns 0). 'affine': gap runs cost gapOpen + k * gapExtend; the pairwise / profile alignments use affine (Gotoh) dynamic programming and totalScore is the sum of pairs with natural affine gap costs (Altschul 1989: each pair of rows is projected out of the alignment, gap-gap columns skipped, every gap run of the projection pays gapOpen once). For two sequences the affine result is an optimal pairwise alignment; for three or more sequences every method is a heuristic. (default: linear) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA. |
| `consensus` | string | Majority-vote consensus of the aligned rows. |
| `totalScore` | integer | Sum-of-pairs score of the MSA: linear model = column score (match/mismatch from the matrix, residue-gap = gapExtend, gap-gap = 0); gapModel='affine' = natural affine gap costs (every gap run of each pair projection pays gapOpen once). |

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

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1758](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1758)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
