# infoseq

EMBOSS 6.6.0 infoseq: Displays basic information about sequences, as infoseq of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `infoseq` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 infoseq: for every sequence the USA, database, name, accession, GI, version, type, length, percent GC (nucleic), organism and description, in columns (or delimited, or as an HTML table), with a heading line. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Report.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.SequenceFormats/EmbossSequencePrograms.Report.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_report_oracle.py](../../../../scripts/alignment/generate_emboss_prog_report_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `columns` | boolean | No | Print the information in aligned columns (-columns, default true). (default: true) |
| `delimiter` | string | No | The delimiter of the records without columns (-delimiter, default \|). |
| `html` | boolean | No | Format the output as an HTML table (-html). (default: false) |
| `only` | boolean | No | Display only the columns that are asked for (-only). (default: false) |
| `heading` | boolean | No | Display the column headings (-heading, default true). (default: true) |
| `usa` | boolean | No | The USA column (-usa). (default: true) |
| `database` | boolean | No | The database column (-database). (default: true) |
| `name` | boolean | No | The name column (-name). (default: true) |
| `accession` | boolean | No | The accession column (-accession). (default: true) |
| `gi` | boolean | No | The GI column (-gi). (default: false) |
| `seqversion` | boolean | No | The version column (-seqversion). (default: false) |
| `type` | boolean | No | The type column (-type). (default: true) |
| `length` | boolean | No | The length column (-length). (default: true) |
| `pgc` | boolean | No | The percent GC column (-pgc). (default: true) |
| `organism` | boolean | No | The organism column (-organism). (default: true) |
| `description` | boolean | No | The description column (-description). (default: true) |
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
> Run infoseq on this input.

**Expected Tool Call:**
```json
{"tool": "infoseq", "arguments": {"sequence": ">a desc here\nACGTACGTAC\n>b\nMKVLAAGIV\n"}}
```

**Response:**
```json
{"text": "USA                      Database  Name           Accession      Type Length %GC    Organism            Description \nfasta::sequence.dat:a    -              a              -              N    10     50.00                      desc here\nfasta::sequence.dat:b    -              b              -              P    9                          \n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
