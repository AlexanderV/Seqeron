# semi_global_align_affine

Semi-global (fitting) pairwise alignment with affine gap costs: sequence1 is aligned end to end inside sequence2 (free leading/trailing sequence2 residues); a gap of length k scores gapOpen + k·gapExtend. Equals Biopython PairwiseAligner(mode='global', end_insertion_score=0) and parasail sg_dx. startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1); startPosition1/endPosition1 span the whole query. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Matrix mode accepts protein sequences (no DNA validation).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `semi_global_align_affine` |
| **Method ID** | `SequenceAligner.SemiGlobalAlignAffine` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Semi-global (fitting) pairwise alignment with affine gap costs: sequence1 is aligned end to end inside sequence2 (free leading/trailing sequence2 residues); a gap of length k scores gapOpen + k·gapExtend. Equals Biopython PairwiseAligner(mode='global', end_insertion_score=0) and parasail sg_dx. startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1); startPosition1/endPosition1 span the whole query. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Matrix mode accepts protein sequences (no DNA validation).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L910](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L910)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | Query sequence (typically shorter). |
| `sequence2` | string | Yes | Reference sequence (typically longer). |
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
