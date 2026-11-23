# global_align_affine

Global pairwise alignment with affine gap costs (Gotoh 1982): a gap of length k scores gapOpen + k·gapExtend (the BLAST+ / Biopython open_gap_score = gapOpen + gapExtend convention), so BLAST scoring (match 2, mismatch -3, gapOpen -5, gapExtend -2) can be applied. Border gaps are charged an opening each (Flouri et al. 2015). Equals Biopython PairwiseAligner(mode='global') and parasail nw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `global_align_affine` |
| **Method ID** | `SequenceAligner.GlobalAlignAffine` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Global pairwise alignment with affine gap costs (Gotoh 1982): a gap of length k scores gapOpen + k·gapExtend (the BLAST+ / Biopython open_gap_score = gapOpen + gapExtend convention), so BLAST scoring (match 2, mismatch -3, gapOpen -5, gapExtend -2) can be applied. Border gaps are charged an opening each (Flouri et al. 2015). Equals Biopython PairwiseAligner(mode='global') and parasail nw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L387](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L387)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (DNA, uppercased internally). |
| `sequence2` | string | Yes | Second sequence. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E. (default: -2) |
| `gapExtend` | number | No | Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. This pairwise tool (it knows its sequences) also writes the alignment in the COORDINATE formats of Biopython 1.88 Alignment.format, byte for byte: bed (bed3 ... bed12), psl, chain, maf, sam (sam:md adds the MD tag), exonerate (exonerate:cigar), and the strings cigar (cigar:eqx = extended = / X; soft clips for the unaligned query), md (SAM MD tag) and cs (cs:long; minimap2 difference string). In those formats sequence1 is the target (reference) and sequence2 the query, except for the semi-global tools where sequence2 is the reference (target); the plain name sam is the Biopython SAM line (the EMBOSS one is 'emboss:sam'); a 'bio:' prefix forces the coordinate writer. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |
| `trimEndGaps` | boolean | No | Coordinate formats (bed, psl, chain, maf, sam, exonerate, cigar, md, cs): drop the leading and trailing columns that have a gap in one row (the free end gaps of a semi-global result), so that the record starts and ends with an aligned pair (Biopython alignment[:, first:last]); false (default) keeps them as D / I blocks as Biopython's PairwiseAligner does (2D4M2D). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Aligned sequence1 (with '-' gaps). |
| `alignedSequence2` | string | Aligned sequence2 (with '-' gaps). |
| `score` | integer | Alignment score. |
| `alignmentType` | string | Alignment type (Global, Local or SemiGlobal). |
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

### Example 1: Long gap, blastn scoring

**Tool Call:**
```json
{
  "tool": "global_align_affine",
  "arguments": {
    "sequence1": "TTTTTTTTTTACGT",
    "sequence2": "ACGT",
    "match": 2,
    "mismatch": -3,
    "gapOpen": -5,
    "gapExtend": -2
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "TTTTTTTTTTACGT",
  "alignedSequence2": "----------ACGT",
  "score": -17,
  "alignmentType": "Global",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 13,
  "endPosition2": 3
}
```

## Worked Example

TTTTTTTTTTACGT vs ACGT with blastn scoring (2/-3/-5/-2): 4 matches (+8) plus one 10-base gap (-5 + 10·-2 = -25) = -17 (Biopython = parasail). The linear `global_align` ignores gapOpen and scores -12.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L387](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L387)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
