# checktrans

EMBOSS 6.6.0 checktrans: Reports the STOP codons and the open reading frames of a protein translation, as checktrans of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `checktrans` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 checktrans: for every protein sequence the stretches between its stops (*) that are at least orfml residues long: the report out.txt lists their positions and the number of stops, the sequences name_N are returned in the file seqout.txt (sequence format outputFormat) and the protein feature table in the file features.txt (feature format featureFormat, default gff3); with addlast (default) a stop is forced at the end (it is not counted). Equals the installed program on 300 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Checktrans.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Checktrans.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_wp63_oracle.py](../../../../scripts/alignment/generate_emboss_prog_wp63_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `orfml` | integer | No | The minimum length of an ORF to report (-orfml, default 100). |
| `addlast` | boolean | No | Force the sequence to end with a stop (-addlast, default true). (default: true) |
| `outputFormat` | string | No | The sequence output format (-osformat, default fasta). |
| `featureFormat` | string | No | The feature output format (-offormat, default gff3): gff, gff2, gff3, bed, pir, swiss, dasgff, draw, debug, refseqp. |
| `sequenceFormat` | string | No | The input format (-sformat1); default: the program finds the format. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that the program may read (a file named by a USA, @file). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text of the main output of the program (the report, or its standard output). |
| `files` | object | The other files the program wrote (sequences, features, graph data), by name. |
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
> Run checktrans on this input.

**Expected Tool Call:**
```json
{"tool": "checktrans", "arguments": {"sequence": ">p1 test\nMKVLAAGIV*ACDEFGHIK*PQRS\n", "orfml": 3}}
```

**Response:**
```json
{"text": "\n\nCHECKTRANS of p1 from 1 to 25\n\n\tORF#\tPos\tLen\tORF Range\tSequence name\n\n\t1\t10\t9\t1-9\tp1_1\n\t2\t20\t9\t11-19\tp1_2\n\t3\t25\t4\t21-24\tp1_3\n\n\tTotal STOPS:     2\n\n ", "files": {"features.txt": "##gff-version 3\n##sequence-region p1 1 24\n#!Date 2023-11-14\n#!Type Protein\n#!Source-version EMBOSS 6.6.0.0\np1\tchecktrans\tpolypeptide_region\t1\t9\t.\t.\t.\tID=p1.1\np1\tchecktrans\tpolypeptide_region\t11\t19\t.\t.\t.\tID=p1.2\np1\tchecktrans\tpolypeptide_region\t21\t24\t.\t.\t.\tID=p1.3\n", "seqout.txt": ">p1_1\nMKVLAAGIV\n>p1_2\nACDEFGHIK\n>p1_3\nPQRS\n"}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
