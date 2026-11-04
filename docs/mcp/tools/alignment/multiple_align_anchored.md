# multiple_align_anchored

Approximate, faster center-star MSA: 4-mer count cosine center, suffix-tree exact-match anchors (MUMmer-style seed, weighted colinear chain, Needleman-Wunsch between anchors). Pairwise alignments are not guaranteed optimal, so the sum-of-pairs score can be lower than multiple_align's. Linear gap model (gapOpen ignored).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_anchored` |
| **Method ID** | `SequenceAligner.MultipleAlignAnchored` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Approximate, faster center-star MSA: 4-mer count cosine center, suffix-tree exact-match anchors (MUMmer-style seed, weighted colinear chain, Needleman-Wunsch between anchors). Pairwise alignments are not guaranteed optimal, so the sum-of-pairs score can be lower than multiple_align's. Linear gap model (gapOpen ignored).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1242](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1242)

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

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1242](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1242)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
