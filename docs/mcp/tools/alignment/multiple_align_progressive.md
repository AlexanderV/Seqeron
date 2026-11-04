# multiple_align_progressive

Guide-tree progressive MSA (Feng & Doolittle 1987 / CLUSTAL): ClustalW pairwise distances from all Needleman-Wunsch alignments (1 - identities / shorter length), a UPGMA guide tree, then profile-profile Needleman-Wunsch with sum-of-pairs column scores (once a gap, always a gap). Linear gap model (gapOpen ignored).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_progressive` |
| **Method ID** | `SequenceAligner.MultipleAlignProgressive` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Guide-tree progressive MSA (Feng & Doolittle 1987 / CLUSTAL): ClustalW pairwise distances from all Needleman-Wunsch alignments (1 - identities / shorter length), a UPGMA guide tree, then profile-profile Needleman-Wunsch with sum-of-pairs column scores (once a gap, always a gap). Linear gap model (gapOpen ignored).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1676](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1676)

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

### Example 1: One deletion

**Tool Call:**
```json
{
  "tool": "multiple_align_progressive",
  "arguments": {
    "sequences": [
      "ACGT",
      "ACGT",
      "AGT"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGT",
    "ACGT",
    "A-GT"
  ],
  "consensus": "ACGT",
  "totalScore": 8
}
```

## Worked Example

ACGT, ACGT merge first; AGT is aligned to that profile with the gap opposite C: SP = +3 (A) - 1 (C,C,-) + 3 (G) + 3 (T) = 8.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1676](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1676)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
