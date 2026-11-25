# newseq

EMBOSS 6.6.0 newseq: Makes a sequence from a name, a description, a type and the residues, as newseq of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `newseq` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 newseq: the sequence is written in the chosen format with the gaps and the ambiguity codes of its type (nucleic: gapnucleotide; protein: gapprotein) as the output sequence check leaves them. Equals the installed program on 320 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP61, F280.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.SequenceFormats/EmbossSequencePrograms.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_seq_oracle.py](../../../../scripts/alignment/generate_emboss_prog_seq_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `name` | string | Yes | The name of the sequence (-name). |
| `description` | string | No | The description of the sequence (-description). (default: ) |
| `type` | string | No | N nucleic or P protein (-type). (default: N) |
| `residues` | string | No | The residues (-sequence). (default: ) |
| `outputFormat` | string | No | The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret. |

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
> Run newseq on this sequence.

**Expected Tool Call:**
```json
{"tool": "newseq", "arguments": {"name": "nn", "description": "d d", "type": "N", "residues": "ACGTRXYZ"}}
```

**Response:**
```json
{"text": ">nn d d\nACGTRXYZ\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
