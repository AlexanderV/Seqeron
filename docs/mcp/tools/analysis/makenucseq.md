# makenucseq

EMBOSS 6.6.0 makenucseq: Creates random nucleotide sequences, as makenucseq of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `makenucseq` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 makenucseq: amount sequences of length lower case nucleotides, each residue picked at random from ACGT (or from the codons of the codon usage file with the frequencies of the file, the end cut to the length), with the insert put at start; the names are EMBOSS_001, EMBOSS_002 .... Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Random.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Random.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_random_oracle.py](../../../../scripts/alignment/generate_emboss_prog_random_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `amount` | integer | No | The number of sequences (-amount, default 100). |
| `length` | integer | No | The length of every sequence (-length, default 100). |
| `insert` | string | No | A string put into every sequence (-insert). |
| `start` | integer | No | The position of the insert (-start, default 1). |
| `codonFile` | string | No | The text of a codon usage file in the EMBOSS (cusp) format (-codonfile). |
| `microseconds` | integer | No | The microseconds of the clock when the random number generator is seeded (the seed is microseconds % 9999 + 1; the program takes them from the clock of the machine): the same value gives the same sequences. Default: the clock. |
| `outputFormat` | string | No | The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret. |

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
> Run makenucseq on this input.

**Expected Tool Call:**
```json
{"tool": "makenucseq", "arguments": {"amount": 2, "length": 20, "microseconds": 7}}
```

**Response:**
```json
{"text": ">EMBOSS_001\ntcactctggtccagtttgcc\n>EMBOSS_002\ncgcaagcttagactaaaatt\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
