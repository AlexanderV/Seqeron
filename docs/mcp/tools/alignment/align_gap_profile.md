# align_gap_profile

Pairwise alignment with per-direction and per-end gap scores (Biopython PairwiseAligner's twelve gap-score properties).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `align_gap_profile` |
| **Method ID** | `SequenceAligner.GlobalAlignAffine(GapScoreProfile)` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Pairwise alignment whose gap scores differ by direction and by position: Biopython 1.88 PairwiseAligner's twelve properties {open,extend}_{internal,left,right}_{insertion,deletion}_score (insertion = a gap in sequence1 = the target, deletion = a gap in sequence2 = the query; the deprecated names target_left_open_gap_score, query_right_extend_gap_score, ... and every aggregate such as open_gap_score, end_extend_gap_score, left_gap_score, insertion_score are accepted too, more specific names override general ones). gapScores = name -> value, Biopython convention (an open score is the score of the FIRST gap position, extend the score of each further one; all <= 0 up to Open + Extend <= 0). A gap run that begins the alignment is a left (start) run, one that finishes it a right (end) run, every other run is internal; unnamed slots use the interior gapOpen / gapExtend of this tool (Seqeron convention: a gap of length k scores gapOpen + k*gapExtend, linear when affine = false: k*gapExtend). A left / right score of 0 is a free overhang (the semi-global family); endGaps (comma-separated flags FreeStart1, FreeEnd1, FreeStart2, FreeEnd2, or the names Fitting, Overlap..., All) makes the matching slots free unless gapScores names them. mode = global (default), semi_global (the result is tagged SemiGlobal when an end is free) or local (a local alignment never starts or ends in a gap: only the internal scores apply, as in Biopython's local mode). affine = false is the linear program (the per-position score is the extend value; open is ignored). With matrix = BLOSUM62 | PAM250 | EDNAFULL ... residues are scored by that table and every gap score may be fractional (exactScore = the real score). Coordinates: a free start / end run is an overhang outside startPosition / endPosition, a priced one is inside. No new dynamic program: the existing linear / Gotoh / Smith-Waterman cores take the gap scores as a function of the matrix position. VERIFICATION: Biopython 1.88 PairwiseAligner score on 1000 random pairs (scalar, EDNAFULL, BLOSUM62; global / local; linear / affine; free, equal and independently priced ends) equal on 998; the two others are a Biopython linear Needleman-Wunsch corner (its first column of the last row is charged the right score nA*right_deletion_extend, its own Gotoh routine the left score). Inputs above 25,000,000 dynamic-programming cells are aligned by the exact linear-space program (Hirschberg / Myers-Miller; the gap score of each sub-problem is resolved by its position in the whole matrix; the same rows as the full matrix, O(m + n) memory).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.GapProfile.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.GapProfile.cs)
- Algorithm doc: [Semi_Global_Alignment.md](../../../algorithms/Alignment/Semi_Global_Alignment.md)
- Oracle: `scripts/alignment/generate_gap_profile_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/gap_profile_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (Biopython target); uppercased internally. |
| `sequence2` | string | Yes | Second sequence (Biopython query). |
| `mode` | string | No | global (default) \| semi_global \| local (case-insensitive). |
| `affine` | boolean | No | true (default) = affine gaps (Gotoh), false = linear (per-position extend score). |
| `match` | integer | No | Score for matching base (ignored when matrix is given). |
| `mismatch` | integer | No | Score for mismatching base (ignored when matrix is given). |
| `gapOpen` | number | No | Interior gap-open score for slots not named in gapScores (Seqeron convention: added once per run; <= 0). May be fractional when matrix is given. |
| `gapExtend` | number | No | Interior gap score per position for slots not named in gapScores (<= 0). May be fractional when matrix is given. |
| `gapScores` | object<string,string> | No | Biopython gap-score properties, name -> value (e.g. {"open_left_insertion_score": 0, "extend_left_insertion_score": 0, "end_deletion_score": -3, "open_gap_score": -5, "extend_gap_score": -1}). |
| `endGaps` | string | No | Optional free end gaps (comma-separated FreeStart1, FreeEnd1, FreeStart2, FreeEnd2 or None, Fitting, FittingSequence2, OverlapSuffix1Prefix2, OverlapSuffix2Prefix1, All); slots named in gapScores keep their score. |
| `matrix` | string | No | Optional substitution matrix by EMBOSS / NCBI name (BLOSUM62, PAM250, EDNAFULL, ...); fractional gap scores are then allowed. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Gapped row of sequence 1 |
| `alignedSequence2` | string | Gapped row of sequence 2 |
| `score` | integer | Alignment score (rounded half away from zero for fractional scores) |
| `alignmentType` | string | Global \| SemiGlobal \| Local |
| `startPosition1 / endPosition1 / startPosition2 / endPosition2` | integer | 0-based inclusive consumed intervals (a free overhang lies outside, a priced end run inside) |
| `exactScore` | number | The real-valued score |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required (null) |
| 1001 | mode must be global, semi_global or local. |
| 1001 | Unknown Biopython gap-score property '...' |
| 1001 | ... must be finite with Extend <= 0 and Open + Extend <= 0 |

## Examples

### Example 1: Free left insertion

**User Prompt:**
> Align AC against GGAC with gap score -5 everywhere except a free gap in sequence 1 at the left end (Biopython left_insertion_score = 0).

**Expected Tool Call:**
```json
{"tool": "align_gap_profile", "arguments": {"sequence1": "AC", "sequence2": "GGAC", "gapScores": {"open_gap_score": -5, "extend_gap_score": -5, "open_left_insertion_score": 0, "extend_left_insertion_score": 0}}}
```

**Response:**
```json
{"alignedSequence1": "--AC", "alignedSequence2": "GGAC", "score": 2, "alignmentType": "SemiGlobal", "startPosition1": 0, "startPosition2": 2, "endPosition1": 1, "endPosition2": 3, "exactScore": 2}
```

Biopython PairwiseAligner(match 1, mismatch -1, gaps -5, left_insertion_score 0).score('AC', 'GGAC') = 2.0; the overhang GG is outside the coordinates.

### Example 2: The reverse pair

**User Prompt:**
> The same scores with the sequences swapped.

**Expected Tool Call:**
```json
{"tool": "align_gap_profile", "arguments": {"sequence1": "GGAC", "sequence2": "AC", "gapScores": {"open_gap_score": -5, "extend_gap_score": -5, "open_left_insertion_score": 0, "extend_left_insertion_score": 0}}}
```

**Response:**
```json
{"alignedSequence1": "GGAC", "alignedSequence2": "--AC", "score": -8, "alignmentType": "SemiGlobal", "startPosition1": 0, "startPosition2": 0, "endPosition1": 3, "endPosition2": 1, "exactScore": -8}
```

A gap in sequence 2 at the left end (deletion) is still charged -5 per position: -10 + 2 = -8 (Biopython 1.88 score = -8.0).

## See Also

- [semi_global_align](semi_global_align.md) - the EndGaps policies (free ends only)
- [global_align_affine](global_align_affine.md)
