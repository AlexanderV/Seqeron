# general_gap_align

Pairwise alignment with an ARBITRARY gap-cost function (Waterman–Smith–Beyer 1976): the score of a gap of length k is gapScores[k-1], so concave (logarithmic, square-root), convex (quadratic), capped, step-wise or any other gap model can be used instead of linear/affine. Every maximal gap run is charged once (a run is never split), the model of Biopython PairwiseAligner with a callable gap_score (1.88, global and local: equal on 960 oracle pairs with 8 gap-function families). mode = global (Needleman–Wunsch), local (Smith–Waterman) or semi_global (free end gaps chosen by endGaps; default Fitting). gapScores are scores (<= 0 = penalty) for k = 1 .. max(length1, length2) (at least that many values; extra values are ignored); a value of -Infinity is not expressible in JSON, use a very large penalty. O(m·n·(m+n)) time: practical to a few hundred residues; a table that is linear or affine (gapScores[k-1] = w1 + (k-1)·e) is recognised and aligned with the O(m·n) linear/affine programs (same score). With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table instead of match/mismatch (protein or IUPAC DNA); exactScore is the real-valued score (score = rounded). Coordinates as global_align / local_align / semi_global_align. Limit: a non-affine gap function needs O(m*n) memory (33 bytes per cell) and O(m*n*(m+n)) time, so (len1+1)*(len2+1) above 4,000,000 cells (about 2,000 x 2,000) is rejected with a clear error; linear and affine functions are delegated to the linear-space programs and have no limit.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `general_gap_align` |
| **Method ID** | `SequenceAligner.GlobalAlignGeneralGap` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Pairwise alignment with an ARBITRARY gap-cost function (Waterman–Smith–Beyer 1976): the score of a gap of length k is gapScores[k-1], so concave (logarithmic, square-root), convex (quadratic), capped, step-wise or any other gap model can be used instead of linear/affine. Every maximal gap run is charged once (a run is never split), the model of Biopython PairwiseAligner with a callable gap_score (1.88, global and local: equal on 960 oracle pairs with 8 gap-function families). mode = global (Needleman–Wunsch), local (Smith–Waterman) or semi_global (free end gaps chosen by endGaps; default Fitting). gapScores are scores (<= 0 = penalty) for k = 1 .. max(length1, length2) (at least that many values; extra values are ignored); a value of -Infinity is not expressible in JSON, use a very large penalty. O(m·n·(m+n)) time: practical to a few hundred residues; a table that is linear or affine (gapScores[k-1] = w1 + (k-1)·e) is recognised and aligned with the O(m·n) linear/affine programs (same score). With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table instead of match/mismatch (protein or IUPAC DNA); exactScore is the real-valued score (score = rounded). Coordinates as global_align / local_align / semi_global_align. Limit: a non-affine gap function needs O(m*n) memory (33 bytes per cell) and O(m*n*(m+n)) time, so (len1+1)*(len2+1) above 4,000,000 cells (about 2,000 x 2,000) is rejected with a clear error; linear and affine functions are delegated to the linear-space programs and have no limit.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.GeneralGap.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.GeneralGap.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (DNA unless matrix is given). |
| `sequence2` | string | Yes | Second sequence (DNA unless matrix is given). |
| `gapScores` | number[] | Yes | Gap scores: gapScores[k-1] is the score of a gap of length k (negative = penalty), for k = 1 .. max(length1, length2) at least. |
| `mode` | string | No | Alignment mode: global, local or semi_global. (default: global) |
| `match` | integer | No | Score for matching base (ignored with matrix). (default: 1) |
| `mismatch` | integer | No | Score for mismatching base (ignored with matrix). (default: -1) |
| `matrix` | string | No | Optional substitution matrix: BLOSUM62, PAM250 (protein) or EDNAFULL (IUPAC DNA). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `endGaps` | string | No | semi_global only: which end gaps are free, see semi_global_align (default 'Fitting'). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Aligned sequence1 ('-' = gap). |
| `alignedSequence2` | string | Aligned sequence2 ('-' = gap). |
| `score` | integer | Optimal score (rounded half away from zero when fractional). |
| `alignmentType` | string | Global, Local or SemiGlobal. |
| `startPosition1` | integer | 0-based start of the aligned segment in sequence1 (-1 for an empty local result). |
| `startPosition2` | integer | 0-based start of the aligned segment in sequence2 (-1 for an empty local result). |
| `endPosition1` | integer | 0-based inclusive end in sequence1. |
| `endPosition2` | integer | 0-based inclusive end in sequence2. |
| `exactScore` | number | Real-valued score (equals score for integer scoring). |

## Errors

| Code | Message |
|------|---------|
| 1010 | Invalid DNA in sequence |

## Examples

### Example 1: Logarithmic (concave) gap cost: one long gap

**Tool Call:**
```json
{
  "tool": "general_gap_align",
  "arguments": {
    "sequence1": "AAAATTTT",
    "sequence2": "AAAACCCCTTTT",
    "gapScores": [
      -2.0,
      -4.0794415417,
      -5.295836866,
      -6.1588830834,
      -6.8283137373,
      -7.3752784077,
      -7.8377304472,
      -8.238324625,
      -8.591673732,
      -8.907755279,
      -9.1936858184,
      -9.4547199494
    ],
    "mode": "global",
    "match": 2,
    "mismatch": -3
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "AAAA----TTTT",
  "alignedSequence2": "AAAACCCCTTTT",
  "score": 10,
  "alignmentType": "Global",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 7,
  "endPosition2": 11,
  "exactScore": 9.841116916640328
}
```

## Worked Example

With w(k) = -(2 + 3 ln k) the four inserted C's are one gap run: 8 matches x 2 = 16, minus w(4) = 2 + 3 ln 4 = 6.1589, so 9.841 (Biopython 1.88 `PairwiseAligner(gap_score=lambda pos, k: -(2 + 3 * log(k)))` gives 9.841116916640328; with the convex w(k) = -k^2 the optimum is 5 and no run of four is used).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.GeneralGap.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.GeneralGap.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
- Waterman MS, Smith TF, Beyer WA (1976) Some biological sequence metrics.
