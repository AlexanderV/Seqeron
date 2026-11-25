# transeq

EMBOSS 6.6.0 transeq: Translates nucleotide sequences to protein in the chosen frames, as transeq of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `transeq` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 transeq: the frames 1 2 3 -1 -2 -3 (F = forward three, R = reverse three, 6 = all; a comma list), the genetic code, regions (translated joined together), trim of trailing X and *, clean (stops to X), the alternative definition of the reverse frames, and the first START codon as M for regions; the translation is the ajtrn.c model (EmbossTrn: ambiguity codes resolved as the EMBOSS genetic code tables do). Equals the installed program on 320 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP61, F280.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.SequenceFormats/EmbossSequencePrograms.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_seq_oracle.py](../../../../scripts/alignment/generate_emboss_prog_seq_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `frame` | string | No | The frames (-frame): 1 2 3 F -1 -2 -3 R 6, separated by commas; default 1. |
| `table` | string | No | The genetic code (-table): 0 .. 6, 9 .. 16, 21 .. 23. |
| `regions` | string | No | The regions to translate (-regions): '4-57,78-94'; default the whole sequence. |
| `trim` | boolean | No | Remove the 'X' and '*' at the end of the translation (-trim). (default: false) |
| `clean` | boolean | No | Change the '*' of the STOP codons to 'X' (-clean). (default: false) |
| `alternative` | boolean | No | Frame -1 starts in the last codon (-alternative). (default: false) |
| `methionine` | boolean | No | A START codon at the beginning of a region translation is an M (-methionine, default true). (default: true) |
| `sequenceFormat` | string | No | The input format (-sformat1); default: the program finds the format. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `outputFormat` | string | No | The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret. |
| `files` | object | No | More files by name (file name -> text) that the program may read (a file named by a USA, @file). |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The standard output of the program: the sequences written to the standard output. |
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
> Run transeq on this sequence.

**Expected Tool Call:**
```json
{"tool": "transeq", "arguments": {"sequence": ">s1 desc here\nATGAAATTTGGGCCCTAGATGCCCTAAGGGATGAAATAA\n", "frame": "1"}}
```

**Response:**
```json
{"text": ">s1_1 desc here\nMKFGP*MP*GMK*\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
