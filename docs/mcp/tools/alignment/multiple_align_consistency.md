# multiple_align_consistency

Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): primary library of global (Needleman-Wunsch) and local (Smith-Waterman) pairwise alignments weighted by percent identity, triplet extension (direct weight + sum of min(W1, W2) over third sequences), progressive alignment along a UPGMA guide tree with the extended weights as position-specific scores and zero gap penalties. The scoring parameters only seed the pairwise library; totalScore is the ordinary column sum-of-pairs score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_consistency` |
| **Method ID** | `SequenceAligner.MultipleAlignConsistency` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): primary library of global (Needleman-Wunsch) and local (Smith-Waterman) pairwise alignments weighted by percent identity, triplet extension (direct weight + sum of min(W1, W2) over third sequences), progressive alignment along a UPGMA guide tree with the extended weights as position-specific scores and zero gap penalties. The scoring parameters only seed the pairwise library; totalScore is the ordinary column sum-of-pairs score.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2425](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2425)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences to align (at least one). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty (ignored: linear-gap model). (default: -2) |
| `gapExtend` | integer | No | Gap penalty charged per gap position. (default: -1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA. |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a nucleotide beats a gap on ties). |
| `totalScore` | integer | Column sum-of-pairs score of the MSA (match/mismatch from the matrix, residue-gap = gapExtend, gap-gap = 0). |

## Errors

| Code | Message |
|------|---------|
| 1010 | Invalid DNA in sequences |

## Examples

### Example 1: Three short sequences

**Tool Call:**
```json
{
  "tool": "multiple_align_consistency",
  "arguments": {
    "sequences": [
      "CTT",
      "ACT",
      "CAA"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "-CTT",
    "AC-T",
    "-CAA"
  ],
  "consensus": "-CAT",
  "totalScore": -3
}
```

## Worked Example

Residue pairs supported by other sequences through the triplet extension (e.g. S1.C-S2.C through S0.C: min(200, 133) = 133) outweigh unsupported pairs; no substitution matrix or gap penalty is used in the progressive DP.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2425](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2425)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
