# semi_global_align

Semi-global (fitting/glocal) pairwise alignment with free end gaps in sequence2; useful for fitting a shorter query into a longer reference.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `semi_global_align` |
| **Method ID** | `SequenceAligner.SemiGlobalAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Semi-global (fitting/glocal) pairwise alignment with free end gaps in sequence2; useful for fitting a shorter query into a longer reference.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L412](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L412)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | Query sequence (typically shorter). |
| `sequence2` | string | Yes | Reference sequence (typically longer). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty. (default: -2) |
| `gapExtend` | integer | No | Gap-extend penalty. (default: -1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Aligned query (leading/trailing '-' for free end gaps). |
| `alignedSequence2` | string | Aligned reference (reproduces the full reference). |
| `score` | integer | Optimal fitting score. |
| `alignmentType` | string | Alignment type ("SemiGlobal"). |
| `startPosition1` | integer | 0-based start in sequence1. |
| `startPosition2` | integer | 0-based start in sequence2. |
| `endPosition1` | integer | 0-based end in sequence1. |
| `endPosition2` | integer | 0-based end in sequence2. |

## Errors

| Code | Message |
|------|---------|
| 1010 | Invalid DNA in sequence |

## Examples

### Example 1: Fit query in reference

**Tool Call:**
```json
{
  "tool": "semi_global_align",
  "arguments": {
    "sequence1": "ACGT",
    "sequence2": "TTACGTTT"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "--ACGT--",
  "alignedSequence2": "TTACGTTT",
  "score": 4,
  "alignmentType": "SemiGlobal",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 3,
  "endPosition2": 7
}
```

## Worked Example

Free end gaps let the query `ACGT` match its exact occurrence in `TTACGTTT` (score 4); the reference flanks are carried as end gaps.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L412](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L412)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
