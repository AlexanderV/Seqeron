# global_align_affine

Global pairwise alignment with affine gap costs (Gotoh 1982): a gap of length k scores gapOpen + k·gapExtend (the BLAST+ / Biopython open_gap_score = gapOpen + gapExtend convention), so BLAST scoring (match 2, mismatch -3, gapOpen -5, gapExtend -2) can be applied. Border gaps are charged an opening each (Flouri et al. 2015). Equals Biopython PairwiseAligner(mode='global') and parasail nw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `global_align_affine` |
| **Method ID** | `SequenceAligner.GlobalAlignAffine` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Global pairwise alignment with affine gap costs (Gotoh 1982): a gap of length k scores gapOpen + k·gapExtend (the BLAST+ / Biopython open_gap_score = gapOpen + gapExtend convention), so BLAST scoring (match 2, mismatch -3, gapOpen -5, gapExtend -2) can be applied. Border gaps are charged an opening each (Flouri et al. 2015). Equals Biopython PairwiseAligner(mode='global') and parasail nw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L387](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L387)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (DNA, uppercased internally). |
| `sequence2` | string | Yes | Second sequence. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E. (default: -2) |
| `gapExtend` | number | No | Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: BLOSUM62, PAM250 (protein) or EDNAFULL (IUPAC DNA). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |

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
| `exactScore` | number | Real-valued score (equals `score` for integer scoring; with matrix / fractional gaps `score` is this value rounded half away from zero). |

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
