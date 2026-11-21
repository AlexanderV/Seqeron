# stretcher_align

EMBOSS 6.6 stretcher: Myers-Miller linear-space global alignment with integer scores.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `stretcher_align` |
| **Method ID** | `SequenceAligner.StretcherAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6 stretcher: global alignment by the Myers–Miller linear-space algorithm (Myers & Miller 1988) with affine gaps and INTEGER scores; a gap of length k costs gapOpen + (k-1)·gapExtend (the needle/EMBOSS convention; EMBOSS -gapopen / -gapextend, positive penalties), end gaps are charged like any gap, residues are upper-cased and a residue that is not in the matrix scores 0 (EMBOSS 'bad character'). The rows are the rows of the installed stretcher 6.6.0 (its midpoint tie rules differ from global_align_affine; the score is the same): VERIFICATION: rows and score identical on 689 / 689 committed random DNA / IUPAC / protein pairs (1–260 residues, gap penalties 0–20 / 0–6, 10 matrices). Default matrix EDNAFULL (nucleotide) or EBLOSUM62 (protein) with gapOpen 16 / gapExtend 4 (nucleotide) or 12 / 2 (protein). Use write_alignment (program 'stretcher') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Stretcher.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Stretcher.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: `scripts/alignment/generate_stretcher_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (stretcher_oracle.tsv)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (-asequence). |
| `sequence2` | string | Yes | Second sequence (-bsequence). |
| `matrix` | string | No | Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62 (protein). |
| `gapOpen` | integer | No | Gap opening penalty (-gapopen, integer >= 0): the cost of a gap of length 1. Omit for 16 (nucleotide) / 12 (protein). |
| `gapExtend` | integer | No | Gap extension penalty (-gapextend, integer >= 0): the cost of each further gap position. Omit for 4 (nucleotide) / 2 (protein). |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers as on the command line, with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 -sprotein1 and the same with 2 for the second sequence (begin / end: 1-based, negative = from the end; -sreverse: reverse-complement a nucleotide sequence, the range then refers to the original; -snucleotide / -sprotein force the type). Measured on EMBOSS 6.6.0: the alignment covers the range, and the formatted alignment (default format 'pair') reports the positions of the ORIGINAL sequence, counting down for a reversed one. Omit for the whole sequences. The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. -sformat -sask -squick -sdbname -sid -iquery -ufo -fformat -fopenfile are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequence1` | string | Row of sequence1. |
| `alignedSequence2` | string | Row of sequence2. |
| `score` | integer | Alignment score. |
| `alignmentType` | string | Global. |
| `startPosition1` | integer | 0. |
| `startPosition2` | integer | 0. |
| `endPosition1` | integer | Length of sequence1. |
| `endPosition2` | integer | Length of sequence2. |
| `exactScore` | number | The score as a real number. |
| `formatted` | string | No | The alignment written in the requested format (null when format is omitted). |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequence1 / sequence2 is null |
| 1002 | unknown matrix; gapOpen or gapExtend negative |

## Examples

### Example 1: DNA pair (stretcher a.fa b.fa)

**Tool Call:**
```json
{
  "tool": "stretcher_align",
  "arguments": {
    "sequence1": "ACGTACGTAGCTAGCTAGCATCGATCGATCGAT",
    "sequence2": "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT"
  }
}
```

**Response:**
```json
{
  "alignedSequence1": "ACGTACGTAGCTAGCTAGCATCGA-TCGATCGAT",
  "alignedSequence2": "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT",
  "score": 140,
  "alignmentType": "Global",
  "startPosition1": 0,
  "startPosition2": 0,
  "endPosition1": 33,
  "endPosition2": 34,
  "exactScore": 140.0
}
```

## See Also

- [global_align_affine](global_align_affine.md)
- [write_alignment](write_alignment.md)
