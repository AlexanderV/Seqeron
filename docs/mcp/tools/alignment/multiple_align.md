# multiple_align

Anchor-based progressive multiple sequence alignment: pick a center sequence by 4-mer cosine similarity, build a suffix tree on it, reconcile per-sequence anchored alignments into an MSA with consensus and sum-of-pairs score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align` |
| **Method ID** | `SequenceAligner.MultipleAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Anchor-based progressive multiple sequence alignment: pick a center sequence by 4-mer cosine similarity, build a suffix tree on it, reconcile per-sequence anchored alignments into an MSA with consensus and sum-of-pairs score.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L702](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L702)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences to align (at least one). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty. (default: -2) |
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

Three identical 8-mers align gaplessly; sum-of-pairs score = C(3,2)=3 pairs * 8 matched columns * (+1) = 24.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L702](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L702)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
