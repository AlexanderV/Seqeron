# local_align_affine

Local pairwise alignment with affine gap costs (Smith–Waterman 1981 with Gotoh 1982 gaps): a gap of length k scores gapOpen + k·gapExtend; zero floor. Equals Biopython PairwiseAligner(mode='local') and parasail sw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `local_align_affine` |
| **Method ID** | `SequenceAligner.LocalAlignAffine` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Local pairwise alignment with affine gap costs (Smith–Waterman 1981 with Gotoh 1982 gaps): a gap of length k scores gapOpen + k·gapExtend; zero floor. Equals Biopython PairwiseAligner(mode='local') and parasail sw_trace. With matrix = BLOSUM62 | PAM250 | EDNAFULL the alignment is scored by that published residue-pair table (protein or IUPAC DNA, any symbol of the matrix) instead of match/mismatch and gapOpen/gapExtend may be fractional (EMBOSS needle/water 6.6.0 = Biopython PairwiseAligner = parasail on >1000 oracle rows; EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E); exactScore is then the real-valued score (score = rounded). Memory: inputs above 25,000,000 dynamic-programming cells ((len1+1)*(len2+1), e.g. two 5 kb sequences) are aligned automatically in exact LINEAR SPACE (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991-style locate-then-solve) - the returned score, rows and coordinates are identical to the full-matrix result (ties included), only the running time grows (about 1.3-2x) and the memory stays O(len1 + len2) instead of O(len1 * len2) (20 kb x 20 kb affine: 111 MB peak instead of 4.8 GB).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L807](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L807)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open score added once per gap run (negative = penalty). May be fractional when matrix is given; EMBOSS gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E. (default: -2) |
| `gapExtend` | number | No | Gap score per gap position, in addition to gapOpen for the first one. May be fractional (e.g. -0.5) when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. This pairwise tool (it knows its sequences) also writes the alignment in the COORDINATE formats of Biopython 1.88 Alignment.format, byte for byte: bed (bed3 ... bed12), psl, chain, maf, sam (sam:md adds the MD tag), exonerate (exonerate:cigar), and the strings cigar (cigar:eqx = extended = / X; soft clips for the unaligned query), md (SAM MD tag) and cs (cs:long; minimap2 difference string). In those formats sequence1 is the target (reference) and sequence2 the query, except for the semi-global tools where sequence2 is the reference (target); the plain name sam is the Biopython SAM line (the EMBOSS one is 'emboss:sam'); a 'bio:' prefix forces the coordinate writer. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |
| `significance` | boolean | No | true: add the BLAST+ statistics of the score (bitScore, eValue, pValue and the Karlin-Altschul parameters) to the result. Protein matrices of the NCBI tables (BLOSUM45 / 50 / 62 / 90, PAM30 / 70 / 250 with the gap costs the tables list, e.g. BLOSUM62 existence 11 extension 1 = gapOpen -11, gapExtend -1) get Spouge's finite-size E-value exactly as blastp -comp_based_stats 0; match / mismatch with the blastn reward / penalty pairs and gap costs of the tables (e.g. 2 / -3 with gapOpen -5, gapExtend -2) the length-adjusted Karlin-Altschul E-value as blastn -task blastn. Costs outside the tables are refused (default false). |
| `dbLength` | integer | No | With significance: the length of the database the subject (sequence2) belongs to, for the E-value (default: the length of sequence2, a pairwise comparison, as blastp / blastn -subject). |
| `dbSequences` | integer | No | With significance on nucleotides: the number of sequences of the database (default 1). |

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
| `significance` | object | No | BLAST+ statistics of the score (only with significance=true): bitScore, eValue, pValue, method, lambda, k, h, lengthAdjustment, searchSpace |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required (null). An empty string is accepted: the library aligns it to all gaps (the local tools return no alignment: score 0 / empty list) |

## Examples

### Example 1: Local window, blastn scoring

**Tool Call:**
```json
{
  "tool": "local_align_affine",
  "arguments": {
    "sequence1": "GGGACGTTTACGTCCC",
    "sequence2": "TTACGTACGTAA",
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
  "alignedSequence1": "TTACGT",
  "alignedSequence2": "TTACGT",
  "score": 12,
  "alignmentType": "Local",
  "startPosition1": 7,
  "startPosition2": 0,
  "endPosition1": 12,
  "endPosition2": 5
}
```

## Worked Example

The best local alignment is TTACGT (6 matches x 2 = 12) at sequence1[7..12] / sequence2[0..5] (Biopython = parasail).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L807](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L807)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
