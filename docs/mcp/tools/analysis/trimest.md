# trimest

EMBOSS 6.6.0 trimest: Removes the poly-A tails of EST sequences, as trimest of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `trimest` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 trimest: the longer of the 3' poly-A tail and the 5' poly-T head (with mismatches allowed; N is ignored) is removed from every sequence, or written in lower case with tolower; a removed poly-T head makes the sequence reverse-complemented unless reverse is false; the description tells what was done. Equals the installed program on 320 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP61, F280.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_seq_oracle.py](../../../../scripts/alignment/generate_emboss_prog_seq_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `minLength` | integer | No | The minimum length of a tail (-minlength, default 4). |
| `mismatches` | integer | No | The number of contiguous mismatches allowed in a tail (-mismatches, default 1). |
| `reverse` | boolean | No | Reverse-complement a sequence whose poly-T head was removed (-reverse, default true). (default: true) |
| `toLower` | boolean | No | Change the tail to lower case instead of removing it (-tolower). (default: false) |
| `fivePrime` | boolean | No | Also look for a 5' poly-T head (-fiveprime, default true). (default: true) |
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
> Run trimest on this sequence.

**Expected Tool Call:**
```json
{"tool": "trimest", "arguments": {"sequence": ">e1 est\nTTTTTTTTACGTACGTNNACGTAAAAAAAAAAAA\n"}}
```

**Response:**
```json
{"text": ">e1 est [poly-A tail removed]\nTTTTTTTTACGTACGTNNACGT\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
