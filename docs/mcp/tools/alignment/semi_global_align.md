# semi_global_align

Performs semi-global (fitting / glocal) pairwise alignment with free end gaps in sequence2; useful for fitting a shorter query into a longer reference. Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use semi_global_align_affine for gapOpen + k·gapExtend). startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1); startPosition1/endPosition1 span the whole query. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Matrix mode accepts protein sequences (no DNA validation).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `semi_global_align` |
| **Method ID** | `SequenceAligner.SemiGlobalAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs semi-global (fitting / glocal) pairwise alignment with free end gaps in sequence2; useful for fitting a shorter query into a longer reference. Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use semi_global_align_affine for gapOpen + k·gapExtend). startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1); startPosition1/endPosition1 span the whole query. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Matrix mode accepts protein sequences (no DNA validation).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L873](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L873)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | Query sequence (typically shorter). |
| `sequence2` | string | Yes | Reference sequence (typically longer). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given. (default: -2) |
| `gapExtend` | number | No | Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: BLOSUM62, PAM250 (protein) or EDNAFULL (IUPAC DNA). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Aligned query (leading/trailing '-' for free end gaps). |
| `alignedSequence2` | string | Aligned reference (reproduces the full reference). |
| `score` | integer | Optimal fitting score. |
| `alignmentType` | string | Alignment type ("SemiGlobal"). |
| `startPosition1` | integer | 0 (the query is aligned in full). |
| `startPosition2` | integer | First sequence2 index consumed by the alignment (fitted-interval start = number of free leading reference residues). |
| `endPosition1` | integer | Query length - 1. |
| `endPosition2` | integer | Last sequence2 index consumed by the alignment (fitted-interval end); an empty interval has endPosition2 = startPosition2 - 1. |

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
  "startPosition2": 2,
  "endPosition1": 3,
  "endPosition2": 5
}
```

## Worked Example

Free end gaps let the query `ACGT` match its exact occurrence in `TTACGTTT` (score 4); the reference flanks are carried as end gaps.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L873](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L873)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
