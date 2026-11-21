# needle_align

EMBOSS 6.6 needle: global alignment with the program's single-precision arithmetic and end-gap options.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `needle_align` |
| **Method ID** | `SequenceAligner.NeedleAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6 needle: global alignment (Needleman-Wunsch with Gotoh's three matrices) in the arithmetic and with the conventions of the program: a gap of length k costs gapOpen + (k-1)*gapExtend (positive penalties, EMBOSS -gapopen / -gapextend, rounded to 16 mantissa bits as EMBOSS does), scores are single-precision floats, and end gaps are FREE unless endWeight (-endweight): then an end gap of length k costs endOpen + (k-1)*endExtend (-endopen / -endextend). Rows (ties included), score and the printed score text are those of the installed needle 6.6.0. Default matrix EDNAFULL (both sequences nucleotide) or EBLOSUM62. exactScore is the single-precision score of the program, score its rounding. VERIFICATION: rows and printed score identical on 887 / 887 committed random DNA / IUPAC / protein pairs (1-200 residues, integer / dyadic / non-dyadic gap penalties, 445 without and 442 with end-gap weights, 10 matrices). Sequences whose matrices would exceed 2^27 = 134 217 728 cells are refused (the program keeps three float matrices and a compass matrix, 16 bytes per cell): use stretcher_align (linear space). Use format / write_alignment (program 'needle') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Needle.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Needle.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: `scripts/alignment/generate_needle_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (needle_oracle.tsv)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (-asequence), at least one residue. |
| `sequence2` | string | Yes | Second sequence (-bsequence), at least one residue. |
| `matrix` | string | No | Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62 (protein). |
| `gapOpen` | number | No | Gap opening penalty (-gapopen, 1 - 100): the cost of a gap of length 1. (default: 10) |
| `gapExtend` | number | No | Gap extension penalty (-gapextend, 0 - 10): the cost of each further gap position. (default: 0.5) |
| `endWeight` | boolean | No | -endweight: apply the end-gap penalties (default false: end gaps are free). |
| `endOpen` | number | No | -endopen (1 - 100): cost of an end gap of length 1 (used only with endWeight). (default: 10) |
| `endExtend` | number | No | -endextend (0 - 10): cost of each further end-gap position (used only with endWeight). (default: 0.5) |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers as on the command line, with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 -sprotein1 and the same with 2 for the second sequence (begin / end: 1-based, negative = from the end; -sreverse: reverse-complement a nucleotide sequence, the range then refers to the original; -snucleotide / -sprotein force the type). Measured on EMBOSS 6.6.0: the alignment covers the range, and the formatted alignment (default format 'pair') reports the positions of the ORIGINAL sequence, counting down for a reversed one. Omit for the whole sequences. The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. -sformat -sask -squick -sdbname -sid -iquery -ufo -fformat -fopenfile are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Row of sequence1. |
| `alignedSequence2` | string | Row of sequence2. |
| `score` | integer | Alignment score (rounded). |
| `alignmentType` | string | Global (needle) / Local (water). |
| `startPosition1` | integer | 0-based first aligned residue of sequence1. |
| `startPosition2` | integer | 0-based first aligned residue of sequence2. |
| `endPosition1` | integer | End of the aligned part in sequence1. |
| `endPosition2` | integer | End of the aligned part in sequence2. |
| `exactScore` | number | The single-precision score of the program. |
| `formatted` | string | The alignment written in the requested format (null when format is omitted). |

## Errors

| Code | Message |
|------|---------|
| 1001 | a sequence argument is null or empty / names and matrix have the wrong shape |
| 1002 | a gap penalty is outside the range of the program / unknown matrix |

## Examples

### Example 1: DNA pair (needle a.fa b.fa)

**Tool Call:**
```json
{
  "tool": "needle_align",
  "arguments": {
    "sequence1": "ACGTACGTAGCTAGCTAGCATCGATCGATCGAT",
    "sequence2": "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "ACGTACGTAGCTAGCTAGCATCG-ATCGATCGAT",
  "alignedSequence2": "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT",
  "score": 146,
  "alignmentType": "Global",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 33,
  "endPosition2": 34,
  "exactScore": 146.0,
  "formatted": null
}
```

## See Also

- [stretcher_align](stretcher_align.md)
- [needleall_align](needleall_align.md)
- [write_alignment](write_alignment.md)
- [global_align_affine](global_align_affine.md)
