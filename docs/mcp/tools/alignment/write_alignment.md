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

Writes an alignment (gapped rows) in an EMBOSS 6.6.0 -aformat output format, byte for byte as the installed programs print it: PAIRWISE formats (two rows; header with program, run date, command line, matrix, gap penalties and Length / Identity / Similarity / Gaps / Score): pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee, unknown, multiple; SEQUENCE-SET formats (any number of rows): fasta, a2m, msf, clustal, mega, meganon, nexus, paup, nexusnon, paupnon, phylip, phylipnon, selex, treecon. VERIFICATION: 3 303 needle / water / stretcher files (29 formats x 114 alignment sets, faked run dates) and 1 309 seqret files (11 sequence-set formats, 1-9 rows) identical to the installed EMBOSS 6.6.0.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/EmbossAlignmentWriter.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/EmbossAlignmentWriter.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: `scripts/alignment/generate_emboss_aformat_oracle.py + generate_emboss_seqout_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (emboss_aformat_oracle.tsv, emboss_seqout_oracle.tsv)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedSequences` | array<string> | Yes | Aligned rows ('-' = gap): two for the pairwise formats, one or more for the sequence-set formats. |
| `format` | string | No | EMBOSS -aformat name. |
| `names` | array<string> | No | Sequence names; default seq1, seq2, ... |
| `descriptions` | array<string> | No | Descriptions (FASTA writer, showDescription). |
| `offsets` | array<integer> | No | Residues of each original sequence before the aligned part. |
| `offends` | array<integer> | No | Residues of each original sequence after the aligned part. |
| `score` | number | No | Alignment score (pairwise formats). |
| `matrix` | string | No | Matrix deciding the similarity marks (EMBOSS / NCBI name); default by sequence type. |
| `gapOpen` | number | No | Gap_penalty printed in the header. |
| `gapExtend` | number | No | Extend_penalty printed in the header. |
| `program` | string | No | needle, water, stretcher, matcher, supermatcher. |
| `width` | integer | No | -awidth (0 = 50). |
| `runDate` | string | No | ISO 8601 run date; default now (UTC). |
| `commandLine` | array<string> | No | Command-line arguments printed in the header. |
| `showAccession` | boolean | No | -aaccshow. |
| `showDescription` | boolean | No | -adesshow. |

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
