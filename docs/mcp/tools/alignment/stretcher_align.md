# stretcher_align

EMBOSS 6.6 stretcher: Myers-Miller linear-space global alignment with integer scores.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `stretcher_align` |
| **Method ID** | `SequenceAligner.StretcherAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6 stretcher: global alignment by the Myers-Miller linear-space algorithm (Myers & Miller 1988) with affine gaps and INTEGER scores; a gap of length k costs gapOpen + (k-1)*gapExtend (the EMBOSS convention), end gaps are charged like any gap, residues are upper-cased, a residue that is not in the matrix scores 0. The rows are the rows of the installed stretcher 6.6.0 (its midpoint tie rules differ from global_align_affine; the score is the same). VERIFICATION: rows and score identical on 689 / 689 committed random DNA / IUPAC / protein pairs (1-260 residues, gap penalties 0-20 / 0-6, 10 matrices).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Stretcher.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Stretcher.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: `scripts/alignment/generate_stretcher_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (stretcher_oracle.tsv)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (-asequence). |
| `sequence2` | string | Yes | Second sequence (-bsequence). |
| `matrix` | string | No | EMBOSS data matrix (-datafile) or NCBI name; default EDNAFULL (nucleotide) / EBLOSUM62 (protein). |
| `gapOpen` | integer | No | -gapopen (>= 0): cost of a gap of length 1; default 16 (nucleotide) / 12 (protein). |
| `gapExtend` | integer | No | -gapextend (>= 0): cost of each further gap position; default 4 (nucleotide) / 2 (protein). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Row of sequence1. |
| `alignedSequence2` | string | Row of sequence2. |
| `score` | integer | Alignment score. |
| `alignmentType` | string | Global. |
| `startPosition1` | integer | 0. |
| `startPosition2` | integer | 0. |
| `endPosition1` | integer | Length of sequence1. |
| `endPosition2` | integer | Length of sequence2. |
| `exactScore` | number | The score as a real number. |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequence1 / sequence2 is null |
| 1002 | unknown matrix; gapOpen or gapExtend negative |

## Examples

### Example 1: DNA pair (stretcher a.fa b.fa)

**Tool Call:**
```json
{
  "tool": "stretcher_align",
  "arguments": {
    "sequence1": "ACGTACGTAGCTAGCTAGCATCGATCGATCGAT",
    "sequence2": "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "ACGTACGTAGCTAGCTAGCATCGA-TCGATCGAT",
  "alignedSequence2": "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT",
  "score": 140,
  "alignmentType": "Global",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 33,
  "endPosition2": 34,
  "exactScore": 140.0
}
```

## See Also

- [global_align_affine](global_align_affine.md)
- [write_alignment](write_alignment.md)
