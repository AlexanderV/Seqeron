# water_align

EMBOSS 6.6 water: local (Smith-Waterman-Gotoh) alignment as the program reports it.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `water_align` |
| **Method ID** | `SequenceAligner.WaterAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6 water: local alignment (Smith-Waterman with Gotoh's affine gaps) as the program computes and reports it: a gap of length k costs gapOpen + (k-1)*gapExtend (positive penalties, EMBOSS -gapopen / -gapextend, rounded to 16 mantissa bits), single-precision scores, the cell of the best score is the LAST one scanned within the float epsilon of the maximum, and the walk-back follows the program's own rules; rows (ties included), start / end coordinates and score are those of the installed water 6.6.0. startPosition1 / startPosition2 are the 0-based first aligned residues. Needs sequences of at least two residues (shorter input returns the empty alignment, as the program's maximum is undefined for a one-row matrix); more than 2^28 matrix cells are refused. Default matrix EDNAFULL (both sequences nucleotide) or EBLOSUM62. VERIFICATION: rows, start coordinates and printed score identical on 860 / 860 committed random DNA / IUPAC / protein pairs (10 matrices). Use format / write_alignment (program 'water') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Water.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Water.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: `scripts/alignment/generate_water_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (water_oracle.tsv)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (-asequence). |
| `sequence2` | string | Yes | Second sequence (-bsequence). |
| `matrix` | string | No | Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62 (protein). |
| `gapOpen` | number | No | Gap opening penalty (-gapopen, 1 - 100): the cost of a gap of length 1. (default: 10) |
| `gapExtend` | number | No | Gap extension penalty (-gapextend, 0 - 10): the cost of each further gap position. (default: 0.5) |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers as on the command line, with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 -sprotein1 and the same with 2 for the second sequence (begin / end: 1-based, negative = from the end; -sreverse: reverse-complement a nucleotide sequence, the range then refers to the original; -snucleotide / -sprotein force the type). Measured on EMBOSS 6.6.0: the alignment covers the range, and the formatted alignment (default format 'pair') reports the positions of the ORIGINAL sequence, counting down for a reversed one. Omit for the whole sequences. The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. -sformat selects the file format of the text; -ufo -fformat -fopenfile name the feature file of a program that takes the features of its sequences (diffseq: the files are given by the tool); -sask -squick -sdbname -sid -iquery are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Row of sequence1. |
| `alignedSequence2` | string | Row of sequence2. |
| `score` | integer | Alignment score (rounded). |
| `alignmentType` | string | Global (needle) / Local (water). |
| `startPosition1` | integer | 0-based first aligned residue of sequence1. |
| `startPosition2` | integer | 0-based first aligned residue of sequence2. |
| `endPosition1` | integer | 0-based inclusive end of the aligned part in sequence1. |
| `endPosition2` | integer | 0-based inclusive end of the aligned part in sequence2. |
| `exactScore` | number | The single-precision score of the program. |
| `formatted` | string | The alignment written in the requested format (null when format is omitted). |

## Errors

| Code | Message |
|------|---------|
| 1001 | a sequence argument is null or empty / names and matrix have the wrong shape |
| 1002 | a gap penalty is outside the range of the program / unknown matrix |

## Examples

### Example 1: Shared 16-mer inside unrelated flanks (water a.fa b.fa)

**Tool Call:**
```json
{
  "tool": "water_align",
  "arguments": {
    "sequence1": "GGGGACGTACGTAGCTAGCTAAAA",
    "sequence2": "TTTTACGTACGTAGCTAGCTCCCC"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "ACGTACGTAGCTAGCT",
  "alignedSequence2": "ACGTACGTAGCTAGCT",
  "score": 80,
  "alignmentType": "Local",
  "startPosition1": 4,
  "startPosition2": 4,
  "endPosition1": 19,
  "endPosition2": 19,
  "exactScore": 80.0,
  "formatted": null
}
```

## See Also

- [local_align_affine](local_align_affine.md)
- [supermatcher_align](supermatcher_align.md)
- [write_alignment](write_alignment.md)
