# prosextract

EMBOSS 6.6.0 prosextract: Builds the PROSITE motif database for patmatmotifs from prosite.dat and prosite.doc, as prosextract of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `prosextract` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 prosextract: for every PATTERN entry of prosite.dat (the text of the SIB file) the file PROSITE/prosite.lines (returned as the text) holds the identifier, the accession, the first word of the description, the pattern and its regular expression; for every {PSnnnnn; NAME} entry of prosite.doc one documentation file PROSITE/PSnnnnn (returned in files, by the path of the data directory /usr/share/EMBOSS/data/PROSITE/PSnnnnn) holds the text between {BEGIN} and {END}. Equals the installed program on 300 runs (every input format of seqret, ranges, the report options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP67, F306.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Prosite.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.Motifs/EmbossMotifPrograms.Prosite.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_prosite_oracle.py](../../../../scripts/alignment/generate_emboss_prog_prosite_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `prositeDat` | string | Yes | The text of the PROSITE data file prosite.dat (ftp.expasy.org/databases/prosite). |
| `prositeDoc` | string | Yes | The text of the PROSITE documentation file prosite.doc. |
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
> Run prosextract on this input.

**Expected Tool Call:**
```json
{"tool": "prosextract", "arguments": {"prositeDat": "ID   A_PAT; PATTERN.\nAC   PS99999;\nDT   01-APR-1990 CREATED.\nDE   Test pattern site.\nPA   A-x(2)-C.\n//\n", "prositeDoc": "{PDOC99999}\n{PS99999; A_PAT}\n{BEGIN}\nDocumentation of A_PAT.\n{END}\n"}}
```

**Response:**
```json
{"text": "A_PAT PS99999\n Test\n A-x(2)-C\n^A[^BJOUXZ]{2}C\n\n", "files": {"/usr/share/EMBOSS/data/PROSITE/PS99999": "Documentation of A_PAT.\n"}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
