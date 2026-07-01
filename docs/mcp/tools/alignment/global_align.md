# global_align

Global pairwise alignment (Needleman-Wunsch): align the full length of both sequences end-to-end under a configurable match/mismatch/gap scoring.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `global_align` |
| **Method ID** | `SequenceAligner.GlobalAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Global pairwise alignment (Needleman-Wunsch): align the full length of both sequences end-to-end under a configurable match/mismatch/gap scoring.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L72](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L72)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (DNA, uppercased internally). |
| `sequence2` | string | Yes | Second sequence. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty. (default: -2) |
| `gapExtend` | integer | No | Gap-extend penalty (linear gap per position). (default: -1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Aligned sequence1 (with '-' gaps). |
| `alignedSequence2` | string | Aligned sequence2 (with '-' gaps). |
| `score` | integer | Optimal global alignment score. |
| `alignmentType` | string | Alignment type ("Global"). |
| `startPosition1` | integer | 0-based start in sequence1. |
| `startPosition2` | integer | 0-based start in sequence2. |
| `endPosition1` | integer | 0-based end in sequence1. |
| `endPosition2` | integer | 0-based end in sequence2. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Identical

**Tool Call:**
```json
{
  "tool": "global_align",
  "arguments": {
    "sequence1": "GATTACA",
    "sequence2": "GATTACA"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "GATTACA",
  "alignedSequence2": "GATTACA",
  "score": 7,
  "alignmentType": "Global",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 6,
  "endPosition2": 6
}
```

### Example 2: One deletion

**Tool Call:**
```json
{
  "tool": "global_align",
  "arguments": {
    "sequence1": "GATTACA",
    "sequence2": "GATACA"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "GATTACA",
  "alignedSequence2": "GA-TACA",
  "score": 5,
  "alignmentType": "Global",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 6,
  "endPosition2": 5
}
```

## Worked Example

With default scoring (match +1), two identical 7-mers align with score 7. `GATTACA` vs `GATACA` needs one gap: score = 6 - 1 = 5.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L72](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L72)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
