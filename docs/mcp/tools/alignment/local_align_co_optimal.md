# local_align_co_optimal

Enumerates ALL distinct local alignments whose score equals the Smith–Waterman optimum (the set Biopython PairwiseAligner(mode='local').align(a, b) iterates; equal on 1 345 oracle pairs), capped at maxAlignments because the number of co-optimal alignments grows exponentially with repeats (tandem repeats, long gap-placement ties). Every alignment ends with an aligned pair, begins right after a zero-score cell, and is not a zero-score extension of a shorter optimal alignment (Biopython rule). Order is deterministic: end cell in row-major order, then diagonal > up > left (affine: M > X > Y) at each branch, so the first alignment equals local_align / local_align_affine. affine = false: linear gaps (gapExtend per position, gapOpen ignored); affine = true: gapOpen + k*gapExtend. With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table and gaps may be fractional (exactScore is the real score). Returns alignments (same fields as local_align); empty list when no pair scores > 0. O(len1*len2) memory: inputs above 25,000,000 cells ((len1+1)*(len2+1)) are rejected with a clear error.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `local_align_co_optimal` |
| **Method ID** | `SequenceAligner.LocalAlignCoOptimal` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Enumerates ALL distinct local alignments whose score equals the Smith–Waterman optimum (the set Biopython PairwiseAligner(mode='local').align(a, b) iterates; equal on 1 345 oracle pairs), capped at maxAlignments because the number of co-optimal alignments grows exponentially with repeats (tandem repeats, long gap-placement ties). Every alignment ends with an aligned pair, begins right after a zero-score cell, and is not a zero-score extension of a shorter optimal alignment (Biopython rule). Order is deterministic: end cell in row-major order, then diagonal > up > left (affine: M > X > Y) at each branch, so the first alignment equals local_align / local_align_affine. affine = false: linear gaps (gapExtend per position, gapOpen ignored); affine = true: gapOpen + k*gapExtend. With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table and gaps may be fractional (exactScore is the real score). Returns alignments (same fields as local_align); empty list when no pair scores > 0. O(len1*len2) memory: inputs above 25,000,000 cells ((len1+1)*(len2+1)) are rejected with a clear error.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.TopN.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TopN.cs)
- Algorithm doc: [Local_Alignment_Smith_Waterman.md](../../../algorithms/Alignment/Local_Alignment_Smith_Waterman.md) (section 5.5)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `maxAlignments` | integer | No | Cap on the number of alignments returned (>= 1). (default: 100) |
| `affine` | boolean | No | false = linear gap model (gapExtend per gap position); true = affine (gapOpen + k*gapExtend). (default: false) |
| `match` | integer | No | Score for matching base (ignored with matrix). (default: 1) |
| `mismatch` | integer | No | Score for mismatching base (ignored with matrix). (default: -1) |
| `gapOpen` | number | No | Gap-open score (affine = true only; ignored by the linear model). May be fractional when matrix is given. (default: -2) |
| `gapExtend` | number | No | Gap score charged per gap position (<= 0). May be fractional when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignments` | array | Local alignments, best first; each has alignedSequence1/2, score, alignmentType, startPosition1/2, endPosition1/2, exactScore (same fields as local_align). |
| `count` | integer | Number of alignments returned. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required (null). An empty string is accepted: the library aligns it to all gaps (the local tools return no alignment: score 0 / empty list) |
| 1002 | maxAlignments must be at least 1 |
| 1003 | Waterman-Eggert declumping and co-optimal enumeration keep the O(m*n) matrices: (len1 + 1) x (len2 + 1) cells exceed the limit of 25,000,000 cells. |

## Examples

### Example 1

**Tool Call:**
```json
{
  "tool": "local_align_co_optimal",
  "arguments": {
    "sequence1": "GAAC",
    "sequence2": "GAAAC",
    "maxAlignments": 100,
    "match": 5,
    "mismatch": -3,
    "gapExtend": -1
  }
}
```

**Response:**
```json
{
  "alignments": [
    {
      "alignedSequence1": "G-AAC",
      "alignedSequence2": "GAAAC",
      "score": 19,
      "alignmentType": "Local",
      "startPosition1": 0,
      "startPosition2": 0,
      "endPosition1": 3,
      "endPosition2": 4,
      "exactScore": 19.0
    },
    {
      "alignedSequence1": "GA-AC",
      "alignedSequence2": "GAAAC",
      "score": 19,
      "alignmentType": "Local",
      "startPosition1": 0,
      "startPosition2": 0,
      "endPosition1": 3,
      "endPosition2": 4,
      "exactScore": 19.0
    },
    {
      "alignedSequence1": "GAA-C",
      "alignedSequence2": "GAAAC",
      "score": 19,
      "alignmentType": "Local",
      "startPosition1": 0,
      "startPosition2": 0,
      "endPosition1": 3,
      "endPosition2": 4,
      "exactScore": 19.0
    }
  ],
  "count": 3
}
```

## Worked Example

The single gap of GAAC/GAAAC can sit in three places (4 matches x 5 - 1 gap = 19; the ungapped alternatives only reach 15); the tool returns all three, the set Biopython PairwiseAligner(mode='local').align(a, b) enumerates.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.TopN.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TopN.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
