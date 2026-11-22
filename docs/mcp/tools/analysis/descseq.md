# descseq

EMBOSS 6.6.0 descseq: the name and / or the description of the FIRST sequence of a file (any input format of seqret) is replaced by name / description, or appended to the existing ones (append), and the sequence is returned in FASTA (the header is the name, the sequence version or accession, and the description). Equals the installed program on 1 000 runs (19 input formats, names and descriptions with blanks, tabs and punctuation, damaged files). A5-8c3g, B09 F250.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `descseq` |
| **Method ID** | `EmbossDescseqPrograms.Descseq` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 descseq: the name and / or the description of the FIRST sequence of a file (any input format of seqret) is replaced by name / description, or appended to the existing ones (append), and the sequence is returned in FASTA (the header is the name, the sequence version or accession, and the description). Equals the installed program on 1 000 runs (19 input formats, names and descriptions with blanks, tabs and punctuation, damaged files). A5-8c3g, B09 F250.

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `name` | string | No | The new name (-name); empty = the name is kept. (default: ) |
| `description` | string | No | The new description (-description); empty = the description is kept. (default: ) |
| `append` | boolean | No | Append the name and the description to the existing ones (-append). (default: false) |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -slower1 -supper1 -snucleotide1 -sprotein1. |
| `files` | object | No | More files by name (file name -> text) that sequence may name with a 'format::' prefix. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceText` | string | The text of the sequence output file (FASTA). |
| `sequences` | integer | The number of sequences written. |
| `messages` | array | What the library printed on its error stream (Warning: ..., Error: ...), in order. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty sequence, a text no sequence reader reads, an input of an alignment or binary format this port does not read). |

## Examples

### Example 1: Replace the name and the description

**User Prompt:**
> Rename the first sequence

**Expected Tool Call:**
```json
{"tool": "descseq", "arguments": {"sequence": ">seq1 old words\nACGTACGT\n>seq2\nTTTT\n", "name": "New", "description": "new text"}}
```

**Response:**
```json
{"sequenceText": ">New new text\nACGTACGT\n", "sequences": 1, "messages": []}
```

## See Also

- [seqret](seqret.md) - read and write sequences
- [maskfeat](maskfeat.md) - mask the regions of features
