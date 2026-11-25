# noreturn

EMBOSS 6.6.0 noreturn: Changes the line ends of a text file, as noreturn of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `noreturn` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 noreturn: the line ends become those of the system (unix LF, pc CRLF, mac CR); a carriage return inside a line is a line end unless the output is mac. Equals the installed program on 320 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP61, F280.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.SequenceFormats/EmbossSequencePrograms.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_seq_oracle.py](../../../../scripts/alignment/generate_emboss_prog_seq_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text` | string | Yes | The text of the file to process. |
| `system` | string | No | unix, pc or mac (-system, default unix). |

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
> Run noreturn on this sequence.

**Expected Tool Call:**
```json
{"tool": "noreturn", "arguments": {"text": "a\rb\r\nc\n"}}
```

**Response:**
```json
{"text": "a\nb\nc\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
