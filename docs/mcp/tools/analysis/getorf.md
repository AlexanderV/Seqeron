# getorf

EMBOSS 6.6.0 getorf: Finds and extracts the open reading frames (ORFs) of nucleotide sequences, as getorf of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `getorf` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 getorf: the regions between STOP codons or between START and STOP codons of the six frames (the circular sequence and the nucleotides around the codons too) are returned as protein or nucleotide sequences named name_1, name_2, ... with the header [start - end] (REVERSE SENSE) and the description of the sequence; the scan is Translator.FindOrfsGetorfProgram (the model of getorf_FindORFs over the EMBOSS genetic code files). Equals the installed program on 800 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP61, F280.

## Core Documentation Reference

- Source: [EmbossGetorfProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossGetorfProgram.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_getorf_oracle.py](../../../../scripts/alignment/generate_emboss_prog_getorf_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `table` | string | No | The genetic code (-table): 0 standard with AUG start only, 1 .. 6, 9 .. 16, 21 .. 23 (the EMBOSS EGC files). |
| `minSize` | integer | No | The minimum size of an ORF in nucleotides (-minsize, default 30). |
| `maxSize` | integer | No | The maximum size of an ORF in nucleotides (-maxsize, default 1000000). |
| `find` | string | No | The type of output (-find): 0 protein between STOP codons, 1 protein between START and STOP, 2 nucleotides between STOP codons, 3 nucleotides between START and STOP, 4 nucleotides around START codons, 5 around initial STOP codons, 6 around ending STOP codons. |
| `methionine` | boolean | No | A START codon at the beginning of a protein is an M (-methionine, default true). (default: true) |
| `circular` | boolean | No | The sequence is circular (-circular). (default: false) |
| `reverse` | boolean | No | ORFs of the reverse complement are found as well (-reverse, default true). (default: true) |
| `flanking` | integer | No | The number of nucleotides on each side of the codon for find 4 .. 6 (-flanking, default 100). |
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
> Run getorf on this sequence.

**Expected Tool Call:**
```json
{"tool": "getorf", "arguments": {"sequence": ">s1 desc here\nATGAAATTTGGGCCCTAGATGCCCTAAGGGATGAAATAA\n", "minSize": 6}}
```

**Response:**
```json
{"text": ">s1_1 [1 - 15] desc here\nMKFGP\n>s1_2 [19 - 24] desc here\nMP\n>s1_3 [5 - 31] desc here\nNLGPRCPKG\n>s1_4 [3 - 38] desc here\nEIWALDALRDEI\n>s1_5 [28 - 36] desc here\nGMK\n>s1_6 [38 - 27] (REVERSE SENSE) desc here\nYFIP\n>s1_7 [23 - 18] (REVERSE SENSE) desc here\nGI\n>s1_8 [14 - 3] (REVERSE SENSE) desc here\nGPNF\n>s1_9 [37 - 2] (REVERSE SENSE) desc here\nISSLRASRAQIS\n>s1_10 [39 - 1] (REVERSE SENSE) desc here\nLFHPLGHLGPKFH\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
