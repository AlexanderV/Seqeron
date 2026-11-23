# semi_global_align

Performs semi-global (fitting / glocal) pairwise alignment with free end gaps in sequence2 by default (endGaps = 'Fitting'; useful for fitting a shorter query into a longer reference); the optional endGaps parameter selects any of the 16 end-gap policies (overlap / suffix-prefix in both orientations, sequence2 inside sequence1, all four ends free = EMBOSS needle's default -endweight N, none = global; Biopython PairwiseAligner left/right_insertion/deletion_score = 0, parasail sg_* variants). Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use semi_global_align_affine for gapOpen + k·gapExtend). With the default endGaps startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1) and startPosition1/endPosition1 span the whole query; for any endGaps they are the consumed intervals of each sequence. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Matrix mode accepts protein sequences (no DNA validation). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `semi_global_align` |
| **Method ID** | `SequenceAligner.SemiGlobalAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs semi-global (fitting / glocal) pairwise alignment with free end gaps in sequence2 by default (endGaps = 'Fitting'; useful for fitting a shorter query into a longer reference); the optional endGaps parameter selects any of the 16 end-gap policies (overlap / suffix-prefix in both orientations, sequence2 inside sequence1, all four ends free = EMBOSS needle's default -endweight N, none = global; Biopython PairwiseAligner left/right_insertion/deletion_score = 0, parasail sg_* variants). Linear gap model: every gap position costs gapExtend; gapOpen is IGNORED by this tool (use semi_global_align_affine for gapOpen + k·gapExtend). With the default endGaps startPosition2/endPosition2 report the fitted interval of sequence2 (0-based inclusive; the free flanks lie outside it; an empty interval has endPosition2 = startPosition2 - 1) and startPosition1/endPosition1 span the whole query; for any endGaps they are the consumed intervals of each sequence. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Matrix mode accepts protein sequences (no DNA validation). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L873](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L873)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | Query sequence (typically shorter; any letters - N, IUPAC, protein - uppercased internally). |
| `sequence2` | string | Yes | Reference sequence (typically longer; any letters, uppercased internally). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open score (ignored by this linear-gap tool; kept for schema compatibility). May be fractional when matrix is given. (default: -2) |
| `gapExtend` | number | No | Gap score charged per gap position (<= 0). May be fractional (e.g. -0.5) when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `endGaps` | string | No | Which end gaps are free (default 'Fitting' = the classic semi-global: free leading/trailing sequence2 residues). Comma-separated, case-insensitive names of the flags FreeStart1 (leading sequence1 residues free), FreeEnd1 (trailing sequence1), FreeStart2 (leading sequence2), FreeEnd2 (trailing sequence2) or of a preset: None (= global, all end gaps penalised; result tagged Global), Fitting (FreeStart2+FreeEnd2; sequence1 inside sequence2; parasail sg_dx, Biopython end_insertion_score=0), FittingSequence2 (FreeStart1+FreeEnd1; sequence2 inside sequence1; sg_qx), OverlapSuffix1Prefix2 (FreeStart1+FreeEnd2: a suffix of sequence1 overlaps a prefix of sequence2; sg_qb_de), OverlapSuffix2Prefix1 (FreeEnd1+FreeStart2; sg_qe_db), All (all four free; parasail sg, EMBOSS needle's default -endweight N). Example: 'FreeStart1,FreeEnd2'. startPosition1/endPosition1 and startPosition2/endPosition2 are the consumed intervals of each sequence (the free flanks lie outside them; an empty interval has end = start - 1); a free flank is shown against '-' in the aligned rows. |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Aligned query (leading/trailing '-' for free end gaps). |
| `alignedSequence2` | string | Aligned reference (reproduces the full reference). |
| `score` | integer | Optimal fitting score. |
| `alignmentType` | string | Alignment type ("SemiGlobal"; "Global" when endGaps is None). |
| `startPosition1` | integer | First sequence1 index consumed by the alignment path: 0 unless leading sequence1 residues are free (FreeStart1 / overlap policies). |
| `startPosition2` | integer | First sequence2 index consumed by the alignment (fitted-interval start = number of free leading reference residues). |
| `endPosition1` | integer | Last sequence1 index consumed by the alignment path: query length - 1 unless trailing sequence1 residues are free. |
| `endPosition2` | integer | Last sequence2 index consumed by the alignment (fitted-interval end); an empty interval has endPosition2 = startPosition2 - 1. |
| `formatted` | string | No | The alignment written in the requested format (null when format is omitted). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required (null). Any letters are accepted (N, IUPAC codes, protein); an empty string is accepted |

## Examples

### Example 1: Fit query in reference

**Tool Call:**
```json
{
  "tool": "semi_global_align",
  "arguments": {
    "sequence1": "ACGT",
    "sequence2": "TTACGTTT"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "--ACGT--",
  "alignedSequence2": "TTACGTTT",
  "score": 4,
  "alignmentType": "SemiGlobal",
  "startPosition1": 0,
  "startPosition2": 2,
  "endPosition1": 3,
  "endPosition2": 5
}
```

### Example 2: Overlap (suffix of sequence1 = prefix of sequence2)

**Tool Call:**
```json
{
  "tool": "semi_global_align",
  "arguments": {
    "sequence1": "GGGGACGTAC",
    "sequence2": "ACGTACTTTT",
    "match": 1,
    "mismatch": -1,
    "gapOpen": 0,
    "gapExtend": -1,
    "endGaps": "OverlapSuffix1Prefix2"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "GGGGACGTAC----",
  "alignedSequence2": "----ACGTACTTTT",
  "score": 6,
  "alignmentType": "SemiGlobal",
  "startPosition1": 4,
  "startPosition2": 0,
  "endPosition1": 9,
  "endPosition2": 5
}
```

## Worked Example

Free end gaps let the query `ACGT` match its exact occurrence in `TTACGTTT` (score 4); the reference flanks are carried as end gaps.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L873](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L873)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
