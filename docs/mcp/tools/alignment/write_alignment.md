# write_alignment

Write an alignment in an EMBOSS -aformat output format (pair, markx0-10, srspair, msf, clustal, mega, nexus, phylip, selex, sam, ...).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `write_alignment` |
| **Method ID** | `EmbossAlignmentWriter.Format` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Writes an alignment (gapped rows) in an EMBOSS 6.6.0 -aformat output format, byte for byte as the installed programs print it: PAIRWISE formats (exactly two rows; header with program, run date, command line, matrix, gap penalties and Length / Identity / Similarity / Gaps / Score): pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee, unknown, multiple; SEQUENCE-SET formats (any number of rows): fasta, a2m, msf, clustal, mega, meganon, nexus, paup, nexusnon, paupnon, phylip, phylipnon, selex, treecon. The pair formats take the score, the matrix (decides the similarity marks and count), the gap penalties and the program name (needle, water, stretcher, matcher, supermatcher: stretcher prints integer scores). VERIFICATION: 3 303 needle / water / stretcher files (29 formats x 114 alignment sets, faked run dates) and 1 309 seqret files (11 sequence-set formats, 1–9 rows) are identical to the installed EMBOSS 6.6.0. Offsets (residues of the original sequences outside the aligned part) are given per row for the local alignments.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/EmbossAlignmentWriter.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/EmbossAlignmentWriter.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: `scripts/alignment/generate_emboss_aformat_oracle.py + generate_emboss_seqout_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (emboss_aformat_oracle.tsv, emboss_seqout_oracle.tsv)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedSequences` | array<string> | Yes | The aligned rows ('-' = gap); two for the pairwise formats, one or more for the sequence-set formats. |
| `format` | string | No | The EMBOSS -aformat name: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee, unknown, multiple, fasta, a2m, msf, clustal, mega, meganon, nexus, paup, nexusnon, paupnon, phylip, phylipnon, selex, treecon. (default: srspair) |
| `names` | array<string> | No | Sequence names (the first word of the FASTA header); default seq1, seq2, ... |
| `descriptions` | array<string> | No | Optional descriptions (one per row; shown by the FASTA writer and with showDescription). |
| `offsets` | array<integer> | No | Residues of each original sequence before the aligned part (water: the local start - 1); default 0. |
| `offends` | array<integer> | No | Residues of each original sequence after the aligned part; default 0. |
| `score` | number | No | The alignment score (pairwise formats). (default: 0) |
| `matrix` | string | No | Substitution matrix deciding the Similarity count and the marks (EMBOSS name or NCBI name). Omit: EDNAFULL for nucleotide rows, EBLOSUM62 otherwise. |
| `gapOpen` | number | No | Gap open penalty printed in the header (Gap_penalty); omit for no line. |
| `gapExtend` | number | No | Gap extension penalty printed in the header (Extend_penalty); omit for no line. |
| `program` | string | No | The program that produced the alignment: needle, water, stretcher, matcher, supermatcher (header, SAM @PG; stretcher = integer scores, -aglobal Y; water = -aglobal N). (default: needle) |
| `width` | integer | No | -awidth: residues per line (0 = 50). (default: 0) |
| `runDate` | string | No | Run date (ISO 8601, e.g. 2026-11-19T21:43:53); default: now (UTC). |
| `commandLine` | array<string> | No | The command line words printed in the header, one per line (e.g. '-asequence a.fa'). |
| `showAccession` | boolean | No | -aaccshow: print (accession) in the header. (default: false) |
| `showDescription` | boolean | No | -adesshow: print the description in the header. (default: false) |
| `dialect` | string | No | 'emboss' (default): the EMBOSS 6.6 writers named above. 'biopython': the Biopython 1.88 Bio.Align writers fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; the descriptions are written, the Clustal header says Biopython 1.88). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The alignment file. |
| `format` | string | The EMBOSS format name written. |
| `dialect` | string | emboss. |

## Errors

| Code | Message |
|------|---------|
| 1001 | alignedSequences is null or empty, a pairwise format without exactly two rows |
| 1002 | unknown format, matrix or malformed runDate |

## Examples

### Example 1: Three rows in FASTA

**Tool Call:**
```json
{
  "tool": "write_alignment",
  "arguments": {
    "alignedSequences": [
      "ACGT-ACGTA",
      "AC-TTACG-A",
      "ACGTTAC--A"
    ],
    "format": "fasta",
    "names": [
      "s1",
      "s2",
      "s3"
    ]
  }
}
```

**Response:**
```json
{
  "text": ">s1\nACGT-ACGTA\n>s2\nAC-TTACG-A\n>s3\nACGTTAC--A\n",
  "format": "fasta",
  "dialect": "emboss"
}
```

## See Also

- [stretcher_align](stretcher_align.md)
- [format_alignment](format_alignment.md)
