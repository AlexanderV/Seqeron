# seqxref

EMBOSS 6.6.0 seqxref: Reports the database cross-references of sequence entries, as seqxref of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `seqxref` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 seqxref: for every entry '#USA: count' and the cross-references db:id (the DR lines of EMBL and Swiss-Prot entries ...) sorted by database name and identifier, equal ones once. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Report.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.SequenceFormats/EmbossSequencePrograms.Report.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_report_oracle.py](../../../../scripts/alignment/generate_emboss_prog_report_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
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
> Run seqxref on this input.

**Expected Tool Call:**
```json
{"tool": "seqxref", "arguments": {"sequence": "ID   TEST_HUMAN    Reviewed;    10 AA.\nAC   P12345;\nDR   PDB; 1ABC; X-ray; 1.9 A; A=1-10.\nDR   EMBL; X1234; -.\nSQ   SEQUENCE   10 AA;  1000 MW;  ABCD CRC64;\n     MKVLAAGIVG\n//\n"}}
```

**Response:**
```json
{"text": "#swiss::sequence.dat:TEST_HUMAN: 2\nEMBL:X1234\nPDB:1ABC\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
