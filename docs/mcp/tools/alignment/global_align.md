# global_align

Performs global pairwise alignment using the Needleman–Wunsch dynamic programming algorithm (aligns full length of both sequences end-to-end). Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use global_align_affine for gapOpen + k·gapExtend). With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `global_align` |
| **Method ID** | `SequenceAligner.GlobalAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs global pairwise alignment using the Needleman–Wunsch dynamic programming algorithm (aligns full length of both sequences end-to-end). Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use global_align_affine for gapOpen + k·gapExtend). With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L63](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L63)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (DNA, uppercased internally). |
| `sequence2` | string | Yes | Second sequence. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given. (default: -2) |
| `gapExtend` | number | No | Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS or NCBI name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500 (protein) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4 (IUPAC DNA). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `format` | string | No | Optional output format; the alignment written in that format is returned in 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name) or EMBOSS -aformat names (pair, simple, srspair, srs, markx0-3, markx10, match, score, sam, tcoffee, msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon). Omit for no formatted text. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |

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
| `exactScore` | number | Real-valued score (equals `score` for integer scoring; with matrix / fractional gaps `score` is this value rounded half away from zero). |
| `formatted` | string | No | The alignment written in the requested format (null when format is omitted). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required (null). An empty string is accepted: the library aligns it to all gaps (the local tools return no alignment: score 0 / empty list) |

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

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L63](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L63)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
