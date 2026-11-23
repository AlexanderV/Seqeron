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

Pairwise alignment with the TWO-PIECE (dual) affine gap model of ksw2 / minimap2 (ksw2_extd2): a gap of length k scores max(gapOpen1 + k*gapExtend1, gapOpen2 + k*gapExtend2) (a short-gap and a long-gap regime; all four gap numbers are scores <= 0, ksw2's q + k*e cost is gapOpen = -q, gapExtend = -e; defaults = minimap2 -A 2 -B 4 -O 4,24 -E 2,1), Gotoh's O(m n) program with two gap states per direction, mode 'global' (end to end) or 'local' (Smith-Waterman), optional diagonal band (ksw2 w: only the cells |i - j| <= band; a global alignment needs |m - n| <= band). VERIFICATION: 1 600 committed cases (400 global, 400 local, 400 banded global, 400 banded local; scores equal Biopython 1.88 PairwiseAligner with the callable gap score and an independent Waterman-Smith-Beyer program, every alignment re-scores to its score). Z-drop / end-bonus extension alignment is the separate tool ksw2_extend.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.DualAffine.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.DualAffine.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md) (section 4)
- Oracle: `scripts/alignment/generate_dual_affine_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/dual_affine_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `mode` | string | No | 'global' (default) or 'local'. |
| `match` | integer | No | Score of an identical pair (minimap2 -A). (default: 2) |
| `mismatch` | integer | No | Score of a different pair (minimap2 -B is the penalty: -4). (default: -4) |
| `gapOpen1` | integer | No | Short regime: gap open score (<= 0; minimap2 -O first value as -4). (default: -4) |
| `gapExtend1` | integer | No | Short regime: gap extension score per position (<= 0; -E first value as -2). (default: -2) |
| `gapOpen2` | integer | No | Long regime: gap open score (<= 0; -O second value as -24). (default: -24) |
| `gapExtend2` | integer | No | Long regime: gap extension score per position (<= 0; -E second value as -1). (default: -1) |
| `band` | integer | No | Band half-width (cells with \|i - j\| <= band); -1 = full matrix. (default: -1) |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. This pairwise tool (it knows its sequences) also writes the alignment in the COORDINATE formats of Biopython 1.88 Alignment.format, byte for byte: bed (bed3 ... bed12), psl, chain, maf, sam (sam:md adds the MD tag), exonerate (exonerate:cigar), and the strings cigar (cigar:eqx = extended = / X; soft clips for the unaligned query), md (SAM MD tag) and cs (cs:long; minimap2 difference string). In those formats sequence1 is the target (reference) and sequence2 the query, except for the semi-global tools where sequence2 is the reference (target); the plain name sam is the Biopython SAM line (the EMBOSS one is 'emboss:sam'); a 'bio:' prefix forces the coordinate writer. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |
| `trimEndGaps` | boolean | No | Coordinate formats (bed, psl, chain, maf, sam, exonerate, cigar, md, cs): drop the leading and trailing columns that have a gap in one row (the free end gaps of a semi-global result), so that the record starts and ends with an aligned pair (Biopython alignment[:, first:last]); false (default) keeps them as D / I blocks as Biopython's PairwiseAligner does (2D4M2D). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Row of sequence1. |
| `alignedSequence2` | string | Row of sequence2. |
| `score` | integer | Alignment score. |
| `alignmentType` | string | Global or Local. |
| `startPosition1` | integer | 0-based start. |
| `startPosition2` | integer | 0-based start. |
| `endPosition1` | integer | 0-based inclusive end (length − 1 for a global result). |
| `endPosition2` | integer | 0-based inclusive end (length − 1 for a global result). |
| `exactScore` | number | Score as a real number. |
| `formatted` | string | No | The alignment written in the requested format (null when format is omitted). |

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
  "endPosition1": 7,
  "endPosition2": 7,
  "exactScore": 16.0
}
```

## See Also

- [global_align_affine](global_align_affine.md)
- [local_align_affine](local_align_affine.md)
