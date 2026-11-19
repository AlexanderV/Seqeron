# dual_affine_align

ksw2 / minimap2-style two-piece (dual) affine gap alignment, global or local, optional band.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `dual_affine_align` |
| **Method ID** | `SequenceAligner.GlobalAlignDualAffine / LocalAlignDualAffine` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Pairwise alignment with the TWO-PIECE (dual) affine gap model of ksw2 / minimap2 (ksw2_extd2): a gap of length k scores max(gapOpen1 + k*gapExtend1, gapOpen2 + k*gapExtend2) (a short-gap and a long-gap regime; all four gap numbers are scores <= 0, ksw2's q + k*e cost is gapOpen = -q, gapExtend = -e; defaults = minimap2 -A 2 -B 4 -O 4,24 -E 2,1), Gotoh's O(m n) program with two gap states per direction, mode 'global' (end to end) or 'local' (Smith-Waterman), optional diagonal band (ksw2 w: only the cells |i - j| <= band; a global alignment needs |m - n| <= band). VERIFICATION: 1 600 committed cases (400 global, 400 local, 400 banded global, 400 banded local; scores equal Biopython 1.88 PairwiseAligner with the callable gap score and an independent Waterman-Smith-Beyer program, every alignment re-scores to its score). Z-drop / end-bonus extension alignment is not included (B09.md A2-29b).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.DualAffine.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.DualAffine.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md) (section 4)
- Oracle: `scripts/alignment/generate_dual_affine_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/dual_affine_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `mode` | string | No | 'global' or 'local'. |
| `match` | integer | No | Score of an identical pair. |
| `mismatch` | integer | No | Score of a different pair. |
| `gapOpen1` | integer | No | Short regime: gap open score (<= 0). |
| `gapExtend1` | integer | No | Short regime: extension score per position (<= 0). |
| `gapOpen2` | integer | No | Long regime: gap open score (<= 0). |
| `gapExtend2` | integer | No | Long regime: extension score per position (<= 0). |
| `band` | integer | No | Band half-width (|i - j| <= band); -1 = full matrix. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Row of sequence1. |
| `alignedSequence2` | string | Row of sequence2. |
| `score` | integer | Alignment score. |
| `alignmentType` | string | Global or Local. |
| `startPosition1` | integer | 0-based start. |
| `startPosition2` | integer | 0-based start. |
| `endPosition1` | integer | End (exclusive). |
| `endPosition2` | integer | End (exclusive). |
| `exactScore` | number | Score as a real number. |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequence is null |
| 1002 | mode is not global / local; a gap score is positive; band < -1; a global band narrower than |m - n| |

## Examples

### Example 1: Identical sequences, band 0

**Tool Call:**
```json
{
  "tool": "dual_affine_align",
  "arguments": {
    "sequence1": "ACGTACGT",
    "sequence2": "ACGTACGT",
    "band": 0
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "ACGTACGT",
  "alignedSequence2": "ACGTACGT",
  "score": 16,
  "alignmentType": "Global",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 8,
  "endPosition2": 8,
  "exactScore": 16.0
}
```

## See Also

- [global_align_affine](global_align_affine.md)
- [local_align_affine](local_align_affine.md)
