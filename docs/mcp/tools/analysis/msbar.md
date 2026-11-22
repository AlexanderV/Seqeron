# msbar

EMBOSS 6.6.0 msbar: Mutates sequences with point, codon and block mutations, as msbar of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `msbar` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 msbar: count times the block, codon (nucleic sequences) and point mutations of the lists are done at random places: insertions, deletions, changes, duplications and moves (1 = a random one of them, 0 = none; a comma list of at most 4); the result is tried again (10 times) when it equals a sequence of otherSequence. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Random.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Random.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_random_oracle.py](../../../../scripts/alignment/generate_emboss_prog_random_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `count` | integer | No | The number of times the mutation operations are done (-count, default 1). |
| `point` | string | No | The point mutations (-point): 0 none, 1 any, 2 insertions, 3 deletions, 4 changes, 5 duplications, 6 moves; a comma list of at most 4 (default 0). |
| `block` | string | No | The block mutations (-block), as point (default 0). |
| `codon` | string | No | The codon mutations (-codon), as point (default 0; nucleic sequences only). |
| `inframe` | boolean | No | Do the codon and block operations in frame (-inframe). (default: false) |
| `minimum` | integer | No | The minimum size of a block mutation (-minimum, default 1). |
| `maximum` | integer | No | The maximum size of a block mutation (-maximum, default 10). |
| `otherSequence` | string | No | Sequences the result must not equal (-othersequence: the text of a sequence file, or a USA such as asis:ACGT). |
| `microseconds` | integer | No | The microseconds of the clock when the random number generator is seeded (the seed is microseconds % 9999 + 1; the program takes them from the clock of the machine): the same value gives the same sequences. Default: the clock. |
| `sequenceFormat` | string | No | The input format (-sformat1); default: the program finds the format. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `outputFormat` | string | No | The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret. |
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
> Run msbar on this input.

**Expected Tool Call:**
```json
{"tool": "msbar", "arguments": {"sequence": ">s1 d\nACGTACGTACGTAAAACCCGGGTTT\n", "point": "2,3", "count": 3, "microseconds": 7}}
```

**Response:**
```json
{"text": ">s1 d\nACGTAGTACGCATAAACCCGGGTTT\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
