# patmatdb

EMBOSS 6.6.0 patmatdb: Searches protein sequences with a PROSITE sequence motif, as patmatdb of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `patmatdb` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 patmatdb: every protein sequence is searched with the PROSITE pattern motif (case-independent; < > anchor, x any residue, [ ] a class, { } excluded residues, (n) (n,m) counts) turned into a regular expression of the library (PCRE 7); the hits are the report out.txt (default format dbmotif: every report format of EMBOSS, the report qualifiers of reportQualifiers); a text that the PROSITE syntax refuses (an unmatched bracket, a character that is no residue) stops the program with the message of the library, a pattern that the regular expression library refuses (a range out of order, a count above 65535, nothing to repeat) makes the reference program die by a segmentation fault and is refused with the message of the library. Equals the installed program on 300 runs (every input format of seqret, ranges, the report options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP67, F306.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Prosite.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Prosite.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_prosite_oracle.py](../../../../scripts/alignment/generate_emboss_prog_prosite_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the protein sequence file (any input format of seqret), or a sequence as letters. |
| `motif` | string | Yes | The PROSITE pattern to search for (-motif, at least 2 characters), e.g. [DE](2)HS{P}X(2)PX(2,4)C. |
| `reportFormat` | string | No | The report format (-rformat, default dbmotif): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq. |
| `reportQualifiers` | string | No | Report qualifiers as typed on the command line: -raccshow -rdesshow -rusashow -norscoreshow -norstrandshow -rmaxall N -rmaxseq N (-rformat is the parameter reportFormat). |
| `sequenceFormat` | string | No | The input format (-sformat1); default: the program finds the format. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that the program may read (a file named by a USA, @file, the PROSITE data files of the data directory). |

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
> Run patmatdb on this input.

**Expected Tool Call:**
```json
{"tool": "patmatdb", "arguments": {"sequence": ">p1 test\nMAAAAKKCAGGCLLAKKC\n", "motif": "A-x(2)-C"}}
```

**Response:**
```json
{"text": "########################################\n# Program: patmatdb\n# Rundate: Tue 14 Nov 2023 22:13:20\n# Commandline: patmatdb\n#    -auto\n#    -sequence sequence.dat\n#    -motif A-x(2)-C\n#    -outfile out.txt\n# Report_format: dbmotif\n# Report_file: out.txt\n########################################\n\n#=======================================\n#\n# Sequence: p1     from: 1   to: 18\n# HitCount: 3\n#\n# Motif: A-X(2)-C\n#\n#=======================================\n\nLength = 4\nStart = position 15 of sequence\nEnd = position 18 of sequence\n\n\nGGCLLAKKC\n     |  |\n    15  18\n\nLength = 4\nStart = position 9 of sequence\nEnd = position 12 of sequence\n\n\nAAKKCAGGCLLAKK\n     |  |\n     9  12\n\nLength = 4\nStart = position 5 of sequence\nEnd = position 8 of sequence\n\n\nMAAAAKKCAGGCL\n    |  |\n    5  8\n\n\n#---------------------------------------\n#---------------------------------------\n\n#---------------------------------------\n# Total_sequences: 1\n# Total_length: 18\n# Reported_sequences: 1\n# Reported_hitcount: 3\n#---------------------------------------\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
