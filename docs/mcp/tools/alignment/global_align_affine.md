# global_align_affine

Global pairwise alignment with affine gap costs (Gotoh 1982): a gap of length k scores gapOpen + k·gapExtend (BLAST+ / Biopython convention), border gaps are charged an opening each (Flouri et al. 2015). Equals Biopython PairwiseAligner(mode='global') and parasail nw_trace.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `global_align_affine` |
| **Method ID** | `SequenceAligner.GlobalAlignAffine` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Global pairwise alignment with affine gap costs (Gotoh 1982): a gap of length k scores gapOpen + k·gapExtend (BLAST+ / Biopython convention), border gaps are charged an opening each (Flouri et al. 2015). Equals Biopython PairwiseAligner(mode='global') and parasail nw_trace.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L387](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L387)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (DNA, uppercased internally). |
| `sequence2` | string | Yes | Second sequence. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open score added once per gap run (negative = penalty). (default: -2) |
| `gapExtend` | integer | No | Gap score per gap position, in addition to gapOpen for the first one. (default: -1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Aligned sequence1 (with '-' gaps). |
| `alignedSequence2` | string | Aligned sequence2 (with '-' gaps). |
| `score` | integer | Alignment score. |
| `alignmentType` | string | Alignment type (Global, Local or SemiGlobal). |
| `startPosition1` | integer | 0-based start in sequence1. |
| `startPosition2` | integer | 0-based start in sequence2. |
| `endPosition1` | integer | 0-based end in sequence1. |
| `endPosition2` | integer | 0-based end in sequence2. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Long gap, blastn scoring

**Tool Call:**
```json
{
  "tool": "global_align_affine",
  "arguments": {
    "sequence1": "TTTTTTTTTTACGT",
    "sequence2": "ACGT",
    "match": 2,
    "mismatch": -3,
    "gapOpen": -5,
    "gapExtend": -2
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "TTTTTTTTTTACGT",
  "alignedSequence2": "----------ACGT",
  "score": -17,
  "alignmentType": "Global",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 13,
  "endPosition2": 3
}
```

## Worked Example

TTTTTTTTTTACGT vs ACGT with blastn scoring (2/-3/-5/-2): 4 matches (+8) plus one 10-base gap (-5 + 10·-2 = -25) = -17 (Biopython = parasail). The linear `global_align` ignores gapOpen and scores -12.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L387](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L387)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
