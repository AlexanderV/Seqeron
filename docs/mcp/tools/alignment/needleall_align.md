# needleall_align

EMBOSS 6.6 needleall: global alignment of every sequence of one set against every sequence of another.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `needleall_align` |
| **Method ID** | `SequenceAligner.NeedleAll` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6 needleall: the global alignment (Needleman-Wunsch with affine gaps, the program's own scoring: gap of length k costs gapOpen + (k-1)*gapExtend, end gaps free unless endWeight) of every sequence of the set A against every sequence of B, for each B sequence every A sequence, exactly as needle aligns each pair (same rows, ties and single-precision scores). A pair whose score is not above minScore is not reported; the program's error-file line ('Alignment score (S) is less than minimum score(M) for sequences A vs B') is returned in errorText. brief = false adds the Longest_/Shortest_ Identity / Similarity percentages to the sub-header of the formatted output. VERIFICATION: output files (every -aformat) and error files identical to the installed needleall 6.6.0 on 260 / 260 committed random DNA / IUPAC / protein sets (1-3 x 1-4 sequences, every option: matrix, gapOpen, gapExtend, endWeight / endOpen / endExtend, minScore, brief). Use write_alignment (program 'needle') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.NeedleAll.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.NeedleAll.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `setA` | array of string | Yes | The -asequence set (first member of each pair). |
| `setB` | array of string | Yes | The -bsequence sequences (second member of each pair). |
| `namesA` | array of string | No | Optional names of the A sequences (default a1, a2, ...). |
| `namesB` | array of string | No | Optional names of the B sequences (default b1, b2, ...). |
| `matrix` | string | No | Optional EMBOSS data matrix (-datafile). Omit for EDNAFULL when every A sequence is nucleotide, EBLOSUM62 otherwise. |
| `gapOpen` | number | No | Gap opening penalty (-gapopen, 1 - 100). (default: 10) |
| `gapExtend` | number | No | Gap extension penalty (-gapextend, 0 - 10). (default: 0.5) |
| `endWeight` | boolean | No | -endweight: apply the end-gap penalties. (default: false) |
| `endOpen` | number | No | -endopen (1 - 100). (default: 10) |
| `endExtend` | number | No | -endextend (0 - 10). (default: 0.5) |
| `minScore` | number | No | -minscore (-10 - 100): minimum score to report an alignment. (default: 0) |
| `brief` | boolean | No | -brief: false adds the Longest_/Shortest_ Identity / Similarity lines to the sub-header. (default: true) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignments` | array of object | Alignments in the program order (for each B sequence every A sequence): names, rows, offsets, score. |
| `count` | integer | Number of alignments. |
| `errorText` | string | The error-file text (pairs not above minScore). |
| `formatted` | string | Always null (use write_alignment for the files). |

## Errors

| Code | Message |
|------|---------|
| 1001 | a sequence is null or empty |
| 1002 | unknown matrix; a penalty or minScore out of range |

## Examples

### Example 1: One identical pair

**Tool Call:**
```json
{
  "tool": "needleall_align",
  "arguments": {
    "setA": [
      "ACGT"
    ],
    "setB": [
      "ACGT"
    ]
  }
}
```

**Response:**
```json
{
  "alignments": [
    {
      "name1": "a1",
      "name2": "b1",
      "alignedSequence1": "ACGT",
      "alignedSequence2": "ACGT",
      "offset1": 0,
      "offset2": 0,
      "offend1": 0,
      "offend2": 0,
      "score": 20.0
    }
  ],
  "count": 1,
  "errorText": ""
}
```

## See Also

- [write_alignment](write_alignment.md)
- [stretcher_align](stretcher_align.md)
