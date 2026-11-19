# global_align_co_optimal

Enumerates ALL distinct co-optimal GLOBAL (Needleman–Wunsch / Gotoh) alignments - the set Biopython PairwiseAligner(mode='global').align(a, b) iterates; equal as a set on 1 920 oracle pairs covering all 16 end-gap policies, linear and affine gaps, scalar and BLOSUM62 scoring with fractional gaps, plus exhaustive enumeration of every column sequence of tiny inputs - capped at maxAlignments because the number of co-optimal alignments grows exponentially with ties (tandem repeats, long gap-placement ties). An alignment is a gapped column sequence; distinct gapped rows are listed once; an alignment is optimal when its score (substitution columns + gaps; leading/trailing gap columns of a FREE end cost 0) is maximal. Order is deterministic: end cell (last row ascending, then last column ascending when end gaps are free), then diagonal > up > left (affine: M > X > Y) at each branch of the traceback, so the first alignment equals global_align / global_align_affine / semi_global_align(_affine) with the same endGaps. affine = false: linear gaps (every gap position costs gapExtend, gapOpen ignored); affine = true: a gap of length k scores gapOpen + k*gapExtend (a gap in one row directly followed by a gap in the other row is two runs). endGaps (default 'None' = global; any of the 16 end-gap policies, as semi_global_align) selects which end gaps are free; startPosition/endPosition are the consumed intervals of each alignment. With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table and gaps may be fractional (exactScore is the real score). An empty sequence has the single all-gap alignment. O(len1*len2) memory: inputs above 25,000,000 cells ((len1+1)*(len2+1)) are rejected with a clear error (use global_align, which falls back to linear space, when one optimal alignment is enough).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `global_align_co_optimal` |
| **Method ID** | `SequenceAligner.GlobalAlignCoOptimal` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Enumerates ALL distinct co-optimal GLOBAL (Needleman–Wunsch / Gotoh) alignments - the set Biopython PairwiseAligner(mode='global').align(a, b) iterates; equal as a set on 1 920 oracle pairs covering all 16 end-gap policies, linear and affine gaps, scalar and BLOSUM62 scoring with fractional gaps, plus exhaustive enumeration of every column sequence of tiny inputs - capped at maxAlignments because the number of co-optimal alignments grows exponentially with ties (tandem repeats, long gap-placement ties). An alignment is a gapped column sequence; distinct gapped rows are listed once; an alignment is optimal when its score (substitution columns + gaps; leading/trailing gap columns of a FREE end cost 0) is maximal. Order is deterministic: end cell (last row ascending, then last column ascending when end gaps are free), then diagonal > up > left (affine: M > X > Y) at each branch of the traceback, so the first alignment equals global_align / global_align_affine / semi_global_align(_affine) with the same endGaps. affine = false: linear gaps (every gap position costs gapExtend, gapOpen ignored); affine = true: a gap of length k scores gapOpen + k*gapExtend (a gap in one row directly followed by a gap in the other row is two runs). endGaps (default 'None' = global; any of the 16 end-gap policies, as semi_global_align) selects which end gaps are free; startPosition/endPosition are the consumed intervals of each alignment. With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table and gaps may be fractional (exactScore is the real score). An empty sequence has the single all-gap alignment. O(len1*len2) memory: inputs above 25,000,000 cells ((len1+1)*(len2+1)) are rejected with a clear error (use global_align, which falls back to linear space, when one optimal alignment is enough).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.GlobalCoOptimal.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.GlobalCoOptimal.cs)
- Algorithm doc: [Global_Alignment_Needleman_Wunsch.md](../../../algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md) (section 5.6)

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
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS or NCBI name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500 (protein) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4 (IUPAC DNA). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `endGaps` | string | No | Which end gaps are free; default 'None' = global (all end gaps penalised). Which end gaps are free (default 'Fitting' = the classic semi-global: free leading/trailing sequence2 residues). Comma-separated, case-insensitive names of the flags FreeStart1 (leading sequence1 residues free), FreeEnd1 (trailing sequence1), FreeStart2 (leading sequence2), FreeEnd2 (trailing sequence2) or of a preset: None (= global, all end gaps penalised; result tagged Global), Fitting (FreeStart2+FreeEnd2; sequence1 inside sequence2; parasail sg_dx, Biopython end_insertion_score=0), FittingSequence2 (FreeStart1+FreeEnd1; sequence2 inside sequence1; sg_qx), OverlapSuffix1Prefix2 (FreeStart1+FreeEnd2: a suffix of sequence1 overlaps a prefix of sequence2; sg_qb_de), OverlapSuffix2Prefix1 (FreeEnd1+FreeStart2; sg_qe_db), All (all four free; parasail sg, EMBOSS needle's default -endweight N). Example: 'FreeStart1,FreeEnd2'. startPosition1/endPosition1 and startPosition2/endPosition2 are the consumed intervals of each sequence (the free flanks lie outside them; an empty interval has end = start - 1); a free flank is shown against '-' in the aligned rows. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignments` | array | Co-optimal alignments in the deterministic order of the traceback (first = global_align / global_align_affine / semi_global_align with the same endGaps); each has alignedSequence1/2, score, alignmentType, startPosition1/2, endPosition1/2, exactScore. |
| `count` | integer | Number of alignments returned. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required (null). An empty string is accepted: the library aligns it to all gaps (the local tools return no alignment: score 0 / empty list) |
| 1002 | maxAlignments must be at least 1 |
| 1003 | Unknown endGaps |
| 1004 | Co-optimal enumeration keeps the O(m*n) matrices: (len1 + 1) x (len2 + 1) cells exceed the limit of 25,000,000 cells. |

## Examples

### Example 1

**Tool Call:**
```json
{
  "tool": "global_align_co_optimal",
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
      "alignmentType": "Global",
      "startPosition1": 0,
      "startPosition2": 0,
      "endPosition1": 3,
      "endPosition2": 4,
      "exactScore": 19
    },
    {
      "alignedSequence1": "GA-AC",
      "alignedSequence2": "GAAAC",
      "score": 19,
      "alignmentType": "Global",
      "startPosition1": 0,
      "startPosition2": 0,
      "endPosition1": 3,
      "endPosition2": 4,
      "exactScore": 19
    },
    {
      "alignedSequence1": "GAA-C",
      "alignedSequence2": "GAAAC",
      "score": 19,
      "alignmentType": "Global",
      "startPosition1": 0,
      "startPosition2": 0,
      "endPosition1": 3,
      "endPosition2": 4,
      "exactScore": 19
    }
  ],
  "count": 3
}
```

## Worked Example

The single gap of GAAC/GAAAC can sit in three places (4 matches x 5 - 1 gap = 19); the tool returns all three, the set Biopython PairwiseAligner(mode='global').align(a, b) enumerates (identical sets on 1 920 oracle pairs, all 16 end-gap policies).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.GlobalCoOptimal.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.GlobalCoOptimal.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
