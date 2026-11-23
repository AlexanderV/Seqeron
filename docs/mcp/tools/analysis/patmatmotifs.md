# patmatmotifs

EMBOSS 6.6.0 patmatmotifs: Scans a protein sequence with the patterns of the PROSITE database, as patmatmotifs of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `patmatmotifs` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 patmatmotifs: the protein sequence is searched with every pattern of the PROSITE database in the form prosextract makes (prositeLines: the text of PROSITE/prosite.lines: identifier, accession, description word, pattern and regular expression of every PATTERN entry; the tool prosextract makes it from prosite.dat), the hits are the report out.txt (default format dbmotif: every report format of EMBOSS, the report qualifiers of reportQualifiers); the simple patterns of post-translational modification sites (myristyl, asn_glycosylation, camp_phospho_site, pkc_phospho_site, ck2_phospho_site, tyr_phospho_site) are skipped unless noprune; with full the documentation of every pattern that has a hit is appended (documentation: accession -> the text of the file PROSITE/PSnnnnn that prosextract makes). Equals the installed program on 300 runs (every input format of seqret, ranges, the report options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP67, F306.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Prosite.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Prosite.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_prosite_oracle.py](../../../../scripts/alignment/generate_emboss_prog_prosite_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the protein sequence file (any input format of seqret), or a sequence as letters. |
| `prositeLines` | string | No | The text of the PROSITE data file prosite.lines made by prosextract (PROSITE/prosite.lines): required, the sandbox of the tool has no PROSITE database. |
| `documentation` | object | No | The documentation files of the PROSITE entries made by prosextract: accession (PS00001) -> text; read with full. |
| `full` | boolean | No | Provide the full documentation of the patterns that match (-full). (default: false) |
| `prune` | boolean | No | Ignore the simple patterns (-prune, default true; false = -noprune). (default: true) |
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
> Run patmatmotifs on this input.

**Expected Tool Call:**
```json
{"tool": "patmatmotifs", "arguments": {"sequence": ">s1 test\nMAKKCAGGCLLAKKCDEFGHIKLMNPQRSTVWY\n", "prositeLines": "MYPAT PS99999\n Test\n A-x(2)-C\n^A[^BJOUXZ]{2}C\n\n", "full": true, "documentation": {"PS99999": "MYPAT documentation line 1\nline 2\n"}}}
```

**Response:**
```json
{"text": "########################################\n# Program: patmatmotifs\n# Rundate: Tue 14 Nov 2023 22:13:20\n# Commandline: patmatmotifs\n#    -auto\n#    -sequence sequence.dat\n#    -full\n#    -outfile out.txt\n# Report_format: dbmotif\n# Report_file: out.txt\n########################################\n\n#=======================================\n#\n# Sequence: s1     from: 1   to: 33\n# HitCount: 3\n#\n# Full: Yes\n# Prune: Yes\n# Data_file: PROSITE/prosite.lines\n#\n#=======================================\n\nLength = 4\nStart = position 12 of sequence\nEnd = position 15 of sequence\n\nMotif = MYPAT\n\nGGCLLAKKCDEFGH\n     |  |\n    12  15\n\nLength = 4\nStart = position 6 of sequence\nEnd = position 9 of sequence\n\nMotif = MYPAT\n\nMAKKCAGGCLLAKK\n     |  |\n     6  9\n\nLength = 4\nStart = position 2 of sequence\nEnd = position 5 of sequence\n\nMotif = MYPAT\n\nMAKKCAGGCL\n |  |\n 2  5\n\n\n#---------------------------------------\n#\n# Motif: MYPAT\n# Count: 3\n# \n# MYPAT documentation line 1\n# line 2\n# \n# ***************\n# \n#\n#---------------------------------------\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
