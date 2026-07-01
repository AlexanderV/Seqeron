# local_align

Local pairwise alignment (Smith-Waterman): find the best-scoring substring alignment with a zero floor, ignoring poorly matching flanks.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `local_align` |
| **Method ID** | `SequenceAligner.LocalAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Local pairwise alignment (Smith-Waterman): find the best-scoring substring alignment with a zero floor, ignoring poorly matching flanks.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L301](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L301)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty. (default: -2) |
| `gapExtend` | integer | No | Gap-extend penalty. (default: -1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Aligned local substring of sequence1. |
| `alignedSequence2` | string | Aligned local substring of sequence2. |
| `score` | integer | Optimal local alignment score. |
| `alignmentType` | string | Alignment type ("Local"). |
| `startPosition1` | integer | 0-based start of the local region in sequence1. |
| `startPosition2` | integer | 0-based start of the local region in sequence2. |
| `endPosition1` | integer | 0-based end of the local region in sequence1. |
| `endPosition2` | integer | 0-based end of the local region in sequence2. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Shared core

**Tool Call:**
```json
{
  "tool": "local_align",
  "arguments": {
    "sequence1": "TTTACGTTTT",
    "sequence2": "GGGACGTGGG"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "ACGT",
  "alignedSequence2": "ACGT",
  "score": 4,
  "alignmentType": "Local",
  "startPosition1": 3,
  "startPosition2": 3,
  "endPosition1": 6,
  "endPosition2": 6
}
```

### Example 2: Identical

**Tool Call:**
```json
{
  "tool": "local_align",
  "arguments": {
    "sequence1": "ACGT",
    "sequence2": "ACGT"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "ACGT",
  "alignedSequence2": "ACGT",
  "score": 4,
  "alignmentType": "Local",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 3,
  "endPosition2": 3
}
```

## Worked Example

The shared core `ACGT` scores 4 and is reported without the differing `TTT.../GGG...` flanks.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L301](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L301)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
