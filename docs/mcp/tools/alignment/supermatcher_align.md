# supermatcher_align

EMBOSS 6.6 supermatcher: Rabin-Karp word seeds and a banded Smith-Waterman of large sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `supermatcher_align` |
| **Method ID** | `SequenceAligner.Supermatcher` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6 supermatcher: approximate local alignment of large sequences. The words of wordLength residues of the TARGET sequences (-bsequence) are searched in every QUERY sequence (-asequence) with the Rabin-Karp multi-pattern search of embword.c; the longest word match (the first of equal length) fixes the diagonal, and a banded Smith-Waterman of `width` cells (embAlignPathCalcSWFast, single precision) along that diagonal, from the sequence starts to the end of the shorter overhang, gives the alignment, reported only when its score exceeds minScore. A pair without a shared word is reported in errorText ('No wordmatch start points for A vs B. No alignment'). The result is the program's heuristic alignment (the best one inside the band), not the Smith-Waterman optimum; use local_align_affine / local_align_top_n for that. Alignments are returned in the program's order (for each query, every target). VERIFICATION: output files (every -aformat) and error files identical to the installed supermatcher 6.6.0 on 260 / 260 committed random DNA / IUPAC / protein sets (1-3 queries x 1-4 targets, 4-600 residues, every option: matrix, gapOpen, gapExtend, wordLength 3-15, width 1-100, minScore). Use write_alignment (program 'supermatcher') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Supermatcher.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Supermatcher.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `querySequences` | array of string | Yes | The -asequence sequences (first member of each pair). |
| `targetSequences` | array of string | Yes | The -bsequence sequence set (second member of each pair). |
| `queryNames` | array of string | No | Optional names of the query sequences (default query1, query2, ...). |
| `targetNames` | array of string | No | Optional names of the target sequences (default target1, ...). |
| `matrix` | string | No | Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62. |
| `gapOpen` | number | No | Gap opening penalty (-gapopen, 0 - 100): the cost of a gap of length 1. (default: 10) |
| `gapExtend` | number | No | Gap extension penalty (-gapextend, 0 - 10). (default: 0.5) |
| `wordLength` | integer | No | Word length of the seed matching (-wordlen, at least 3). (default: 6) |
| `width` | integer | No | Alignment band width (-width, at least 1). (default: 16) |
| `minScore` | number | No | Minimum alignment score to report an alignment (-minscore). (default: 0) |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers as on the command line, with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 -sprotein1 and the same with 2 for the second sequence (begin / end: 1-based, negative = from the end; -sreverse: reverse-complement a nucleotide sequence, the range then refers to the original; -snucleotide / -sprotein force the type). Measured on EMBOSS 6.6.0: the alignment covers the range, and the formatted alignment (default format 'pair') reports the positions of the ORIGINAL sequence, counting down for a reversed one. Omit for the whole sequences. Here 1 = the query sequences, 2 = the target sequences (every member takes the same qualifiers). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignments` | array of object | Alignments in the program order: name1, name2, alignedSequence1, alignedSequence2, offset1, offset2, offend1, offend2 (residues before / after the aligned part), score. |
| `count` | integer | Number of alignments. |
| `errorText` | string | The error-file text (pairs without a word match). |
| `formatted` | string | Always null (use write_alignment for the files). |

## Errors

| Code | Message |
|------|---------|
| 1001 | a sequence is null or empty |
| 1002 | unknown matrix; wordLength < 3; width < 1; gapOpen / gapExtend / minScore out of range |

## Examples

### Example 1: Identical 28-mer (one word-seeded alignment)

**Tool Call:**
```json
{
  "tool": "supermatcher_align",
  "arguments": {
    "querySequences": [
      "ACGTACGTACGGTTAACCGGTTAAACCC"
    ],
    "targetSequences": [
      "ACGTACGTACGGTTAACCGGTTAAACCC"
    ]
  }
}
```

**Response:**
```json
{
  "alignments": [
    {
      "name1": "query1",
      "name2": "target1",
      "alignedSequence1": "ACGTACGTACGGTTAACCGGTTAAACCC",
      "alignedSequence2": "ACGTACGTACGGTTAACCGGTTAAACCC",
      "offset1": 0,
      "offset2": 0,
      "offend1": 0,
      "offend2": 0,
      "score": 140.0
    }
  ],
  "count": 1,
  "errorText": ""
}
```

## See Also

- [write_alignment](write_alignment.md)
- [stretcher_align](stretcher_align.md)
