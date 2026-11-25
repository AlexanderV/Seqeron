# listor

EMBOSS 6.6.0 listor: Writes a list file of the logical OR / AND / XOR / NOT of two sets of sequences, as listor of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `listor` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 listor: the USAs of the sequences of the first set and the second set that are (OR: all, AND: in both, XOR: in one only, NOT: in the first and not in the second); sequences are equal when their residues are equal ignoring the case. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Sets.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.SequenceFormats/EmbossSequencePrograms.Sets.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_report_oracle.py](../../../../scripts/alignment/generate_emboss_prog_report_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `firstSequences` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `secondSequences` | string | No | The text of the second sequence file (any input format of seqret), or a sequence as letters. (default: ) |
| `operator` | string | No | OR, AND, XOR or NOT (-operator, default OR; a unique abbreviation is accepted). |
| `sequenceFormat` | string | No | The input format (-sformat1); default: the program finds the format. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that the program may read (a file named by a USA, @file). |

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
> Run listor on this input.

**Expected Tool Call:**
```json
{"tool": "listor", "arguments": {"firstSequences": ">a\nAC\n>b\nGT\n", "secondSequences": ">c\nac\n>d\nTT\n", "operator": "AND"}}
```

**Response:**
```json
{"text": "fasta::firstsequences.dat:a\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
