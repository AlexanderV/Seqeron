# local_align_affine

Local pairwise alignment with affine gap costs (Smith–Waterman 1981 with Gotoh 1982 gaps): a gap of length k scores gapOpen + k·gapExtend; zero floor. Equals Biopython PairwiseAligner(mode='local') and parasail sw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `local_align_affine` |
| **Method ID** | `SequenceAligner.LocalAlignAffine` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Local pairwise alignment with affine gap costs (Smith–Waterman 1981 with Gotoh 1982 gaps): a gap of length k scores gapOpen + k·gapExtend; zero floor. Equals Biopython PairwiseAligner(mode='local') and parasail sw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L807](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L807)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
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

### Example 1: Local window, blastn scoring

**Tool Call:**
```json
{
  "tool": "local_align_affine",
  "arguments": {
    "sequence1": "GGGACGTTTACGTCCC",
    "sequence2": "TTACGTACGTAA",
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
  "alignedSequence1": "TTACGT",
  "alignedSequence2": "TTACGT",
  "score": 12,
  "alignmentType": "Local",
  "startPosition1": 7,
  "startPosition2": 0,
  "endPosition1": 12,
  "endPosition2": 5
}
```

## Worked Example

The best local alignment is TTACGT (6 matches x 2 = 12) at sequence1[7..12] / sequence2[0..5] (Biopython = parasail).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L807](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L807)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
