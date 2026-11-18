# semi_global_align_affine

Semi-global (fitting) pairwise alignment with affine gap costs: sequence1 is aligned end to end inside sequence2 (free leading/trailing sequence2 residues); a gap of length k scores gapOpen + k·gapExtend. Equals Biopython PairwiseAligner(mode='global', end_insertion_score=0) and parasail sg_dx. startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1); startPosition1/endPosition1 span the whole query.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `semi_global_align_affine` |
| **Method ID** | `SequenceAligner.SemiGlobalAlignAffine` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Semi-global (fitting) pairwise alignment with affine gap costs: sequence1 is aligned end to end inside sequence2 (free leading/trailing sequence2 residues); a gap of length k scores gapOpen + k·gapExtend. Equals Biopython PairwiseAligner(mode='global', end_insertion_score=0) and parasail sg_dx. startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1); startPosition1/endPosition1 span the whole query.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L910](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L910)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | Query sequence (typically shorter). |
| `sequence2` | string | Yes | Reference sequence (typically longer). |
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
| `startPosition1` | integer | 0 (the query is aligned in full). |
| `startPosition2` | integer | First sequence2 index consumed by the alignment (fitted-interval start = number of free leading reference residues). |
| `endPosition1` | integer | Query length - 1. |
| `endPosition2` | integer | Last sequence2 index consumed by the alignment (fitted-interval end); an empty interval has endPosition2 = startPosition2 - 1. |

## Errors

| Code | Message |
|------|---------|
| 1010 | Invalid DNA in sequence1/sequence2 |

## Examples

### Example 1: Fit GATTACA, blastn scoring

**Tool Call:**
```json
{
  "tool": "semi_global_align_affine",
  "arguments": {
    "sequence1": "GATTACA",
    "sequence2": "CCGATTTACAGG",
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
  "alignedSequence1": "--GA-TTACA--",
  "alignedSequence2": "CCGATTTACAGG",
  "score": 7,
  "alignmentType": "SemiGlobal",
  "startPosition1": 0,
  "startPosition2": 2,
  "endPosition1": 6,
  "endPosition2": 9
}
```

## Worked Example

GATTACA fitted into CCGATTTACAGG: 6 matches (+12) and one 1-base gap (-5 - 2 = -7) -> 7 (Biopython = parasail).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L910](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L910)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
