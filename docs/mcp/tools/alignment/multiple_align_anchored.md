# multiple_align_anchored

Approximate, faster center-star MSA: picks the center by 4-mer count cosine similarity, builds a suffix tree on it and aligns every other sequence with exact-match anchors (MUMmer-style seed, weighted colinear chain, Needleman–Wunsch — or Gotoh affine with gapModel='affine' — on the gaps between anchors). The pairwise alignments are NOT guaranteed optimal, so the sum-of-pairs score can be lower than multiple_align's. Default linear gap model (gapOpen ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_anchored` |
| **Method ID** | `SequenceAligner.MultipleAlignAnchored` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Approximate, faster center-star MSA: picks the center by 4-mer count cosine similarity, builds a suffix tree on it and aligns every other sequence with exact-match anchors (MUMmer-style seed, weighted colinear chain, Needleman–Wunsch — or Gotoh affine with gapModel='affine' — on the gaps between anchors). The pairwise alignments are NOT guaranteed optimal, so the sum-of-pairs score can be lower than multiple_align's. Default linear gap model (gapOpen ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1821](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1821)

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
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a nucleotide beats a gap on ties). |
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
  "tool": "multiple_align_anchored",
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

## Worked Example

Three identical 8-mers align gaplessly; SP = 3 pairs x 8 columns x (+1) = 24.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1821](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1821)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
