# multiple_align_iterative

Progressive MSA followed by MUSCLE-style tree-dependent restricted-partitioning refinement (Edgar 2004, stage 3): each guide-tree edge splits the alignment into two profiles that are re-aligned and kept only if the sum-of-pairs score strictly improves; repeated until convergence or maxIterations passes. Never scores below multiple_align_progressive. Linear gap model (gapOpen ignored).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_iterative` |
| **Method ID** | `SequenceAligner.MultipleAlignIterative` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Progressive MSA followed by MUSCLE-style tree-dependent restricted-partitioning refinement (Edgar 2004, stage 3): each guide-tree edge splits the alignment into two profiles that are re-aligned and kept only if the sum-of-pairs score strictly improves; repeated until convergence or maxIterations passes. Never scores below multiple_align_progressive. Linear gap model (gapOpen ignored).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2111](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2111)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences to align (at least one). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty (ignored: linear-gap model). (default: -2) |
| `gapExtend` | integer | No | Gap penalty charged per gap position. (default: -1) |
| `maxIterations` | integer | No | Upper bound on full refinement passes (>= 1). (default: 16) |

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
| 1002 | maxIterations must be >= 1 |

## Examples

### Example 1: Gap relocation

**Tool Call:**
```json
{
  "tool": "multiple_align_iterative",
  "arguments": {
    "sequences": [
      "TCA",
      "TCT",
      "GCAC",
      "TAAC"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "TCA-",
    "TCT-",
    "GCAC",
    "TAAC"
  ],
  "consensus": "TCAC",
  "totalScore": -3
}
```

## Worked Example

The progressive seed T--CA / T--CT / GCAC- / TAAC- scores -7; refinement moves the gap block to the last column: TCA- / TCT- / GCAC / TAAC, SP -3.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2111](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2111)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
