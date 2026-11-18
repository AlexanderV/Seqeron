# local_align

Performs local pairwise alignment using the Smith–Waterman algorithm (finds best-scoring substring alignment, with zero floor). Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use local_align_affine for gapOpen + k·gapExtend). With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `local_align` |
| **Method ID** | `SequenceAligner.LocalAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs local pairwise alignment using the Smith–Waterman algorithm (finds best-scoring substring alignment, with zero floor). Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use local_align_affine for gapOpen + k·gapExtend). With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L654](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L654)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given. (default: -2) |
| `gapExtend` | number | No | Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: BLOSUM62, PAM250 (protein) or EDNAFULL (IUPAC DNA). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |

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

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L654](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L654)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
