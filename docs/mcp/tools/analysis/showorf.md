# showorf

EMBOSS 6.6.0 showorf: Displays a nucleotide sequence and its translations in pretty format, as showorf of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `showorf` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 showorf: the sequence in blocks of width bases with a ruler and the numbers of the bases, and under it the translations of the frames (F1 F2 F3 forward, R1 R2 R3 reverse: the amino acid under the first base of every codon, the numbers of the first and last residue of the block); the genetic code is the EMBOSS table of the code number. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Report.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Report.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_report_oracle.py](../../../../scripts/alignment/generate_emboss_prog_report_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (the first sequence is used), or a sequence as letters. |
| `frames` | string | No | The frames to translate (-frames): a comma list of 0 none, 1 F1, 2 F2, 3 F3, 4 R1, 5 R2, 6 R3 (default 1,2,3,4,5,6). |
| `table` | string | No | The genetic code (-table): 0 .. 6, 9 .. 16, 21 .. 23. |
| `width` | integer | No | The width of the screen (-width, default 50; at least 10). |
| `ruler` | boolean | No | Add a ruler (-ruler, default true). (default: true) |
| `plabel` | boolean | No | Number the translations (-plabel, default true). (default: true) |
| `nlabel` | boolean | No | Number the sequence (-nlabel, default true). (default: true) |
| `sequenceFormat` | string | No | The input format (-sformat1); default: the program finds the format. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that the program may read (a file named by a USA, @file). |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The standard output of the program (the sequences) or the text of its output file. |
| `files` | object | The files the program wrote, by name. |
| `messages` | string[] | What the library printed on the error stream (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequence cannot be empty |
| 1001 | The program died: the messages of the library (a value it refuses, an input no reader reads, a sequence of the wrong kind) |
| 1001 | An input on which the reference program never stops or dies by a fault is refused |

## Examples

### Example 1: A run on a small input

**User Prompt:**
> Run showorf on this input.

**Expected Tool Call:**
```json
{"tool": "showorf", "arguments": {"sequence": ">s1 d\nATGAAATTTGGGCCCTAGATGCCCAAATTT\n", "frames": "1,4", "width": 30}}
```

**Response:**
```json
{"text": "SHOWORF of s1 from 1 to 30\n\n           ---------|---------|---------|\n         1 ATGAAATTTGGGCCCTAGATGCCCAAATTT 30\nF1       1 M  K  F  G  P  *  M  P  K  F   4\nR1      10   H  F  K  P  G  L  H  G  F  K 1\n\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
