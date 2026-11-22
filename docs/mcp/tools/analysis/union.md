# union

EMBOSS 6.6.0 union: the sequences of a file (any input format of seqret) are joined into ONE sequence that carries the name, description, accession and type of the first; with findOverlap the bases by which a sequence repeats the end of the one before (the longest end of the previous sequence that is the start of this one) are joined once, and with overlapFile the overlap of every sequence is returned (one number per sequence, the first is 0). The result is written in FASTA or in the output format of outputFormat (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig, gde, fitch, acedb, mase). The features (-feature, -source) are not joined by this tool yet. Equals the installed program on 200 runs (fragments that overlap by random lengths, the reader inputs, damaged files, every output format, the format-for-sequence-type errors). A5-8c3h, B09 F253.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `union` |
| **Method ID** | `EmbossUnionProgram.Union` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 union: the sequences of a file (any input format of seqret) are joined into ONE sequence that carries the name, description, accession and type of the first; with findOverlap the bases by which a sequence repeats the end of the one before (the longest end of the previous sequence that is the start of this one) are joined once, and with overlapFile the overlap of every sequence is returned (one number per sequence, the first is 0). The result is written in FASTA or in the output format of outputFormat (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig, gde, fitch, acedb, mase). The features (-feature, -source) are not joined by this tool yet. Equals the installed program on 200 runs (fragments that overlap by random lengths, the reader inputs, damaged files, every output format, the format-for-sequence-type errors). A5-8c3h, B09 F253.

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `findOverlap` | boolean | No | Look for the overlaps of consecutive sequences when joining (-findoverlap). (default: false) |
| `overlapFile` | boolean | No | Return the overlap of every sequence (-overlapfile). (default: false) |
| `outputFormat` | string | No | The sequence output format (-osformat); default fasta. (default: fasta) |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -slower1 -supper1 -snucleotide1 -sprotein1. |
| `files` | object | No | More files by name (file name -> text) that sequence may name with a 'format::' prefix. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceText` | string | The text of the sequence output file. |
| `overlapText` | string | The overlap of every sequence, one number per line (empty without overlapFile). |
| `messages` | array | What the library printed on its error stream (Warning: ..., Error: ...), in order. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty sequence, a text no sequence reader reads, an input of an alignment or binary format this port does not read). |

## Examples

### Example 1: Join overlapping fragments

**User Prompt:**
> Join these fragments and count the overlaps

**Expected Tool Call:**
```json
{"tool": "union", "arguments": {"sequence": ">f1\nACGTACGTAC\n>f2\nCGTACGTACGGG\n>f3\nGGGTTTAAA\n", "findOverlap": true, "overlapFile": true}}
```

**Response:**
```json
{"sequenceText": ">f1\nACGTACGTACGGGTTTAAA\n", "overlapText": "0\n9\n3\n", "messages": []}
```

## See Also

- [seqret](seqret.md) - read and write sequences
- [maskfeat](maskfeat.md) - mask the regions of features
- [descseq](descseq.md) - alter the name or description of a sequence
