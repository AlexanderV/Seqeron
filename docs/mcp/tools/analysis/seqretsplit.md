# seqretsplit

EMBOSS 6.6.0 seqretsplit: the sequences of a file (any input format of seqret) are written each to a file of its own, name.extension in lower case (the extension is the output format as typed, or outputExtension; outputDirectory is the directory), in FASTA or in the output format of outputFormat (the formats of the union tool). A sequence whose name an earlier sequence used (in any case) gets the name name.002, name.003 ... and the library warns 'Duplicate name ... changed to ...'. With outputSingle false the sequences go to one file (sequenceText), as seqret writes them. Equals the installed program on 200 runs (duplicate names in any case, the reader inputs, damaged files, every output format, directories and extensions). A5-8c3h, B09 F254.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `seqretsplit` |
| **Method ID** | `EmbossSeqretPrograms.Seqretsplit` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 seqretsplit: the sequences of a file (any input format of seqret) are written each to a file of its own, name.extension in lower case (the extension is the output format as typed, or outputExtension; outputDirectory is the directory), in FASTA or in the output format of outputFormat (the formats of the union tool). A sequence whose name an earlier sequence used (in any case) gets the name name.002, name.003 ... and the library warns 'Duplicate name ... changed to ...'. With outputSingle false the sequences go to one file (sequenceText), as seqret writes them. Equals the installed program on 200 runs (duplicate names in any case, the reader inputs, damaged files, every output format, directories and extensions). A5-8c3h, B09 F254.

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `firstOnly` | boolean | No | Read one sequence and stop (-firstonly). (default: false) |
| `outputFormat` | string | No | The sequence output format (-osformat); default fasta. (default: fasta) |
| `outputSingle` | boolean | No | Write every sequence to a file of its own (-ossingle2; the default of the program). (default: true) |
| `outputDirectory` | string | No | The directory of the files (-osdirectory2). |
| `outputExtension` | string | No | The extension of the files (-osextension2); default: the output format as typed. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that sequence may name with a 'format::' prefix. |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence) as the program reads it, instead of the file the text is put in: [format::]file[begin:end:r] (file = the name of an entry of files; the text of sequence is the file 'sequence.dat'), file:id / file:field:query / file{a\|b} (id, acc, sv, gi, des, key, org; wildcards * and ?), @list / list::list (a list file of files; # comments; every entry has its own range and format), asis::SEQUENCE; the range qualifiers -sbegin1 -send1 -sreverse1 of qualifiers are applied to the sequences and to their features. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceText` | string | The text of the sequence output file (empty when every sequence has a file of its own). |
| `singleFiles` | object | The file of every sequence by path (directory/name.extension). |
| `sequences` | integer | The number of sequences written. |
| `messages` | array | What the library printed on its error stream (Warning: ..., Error: ...), in order. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty sequence, a text no sequence reader reads, an input of an alignment or binary format this port does not read). |

## Examples

### Example 1: Split a file with a repeated name

**User Prompt:**
> Write every sequence to its own file

**Expected Tool Call:**
```json
{"tool": "seqretsplit", "arguments": {"sequence": ">a\nACGT\n>A\nGG\n>b\nTT\n"}}
```

**Response:**
```json
{"sequenceText": "", "singleFiles": {"a.fasta": ">a\nACGT\n", "a.002.fasta": ">A.002\nGG\n", "b.fasta": ">b\nTT\n"}, "sequences": 3, "messages": ["Warning: Duplicate name 'A' changed to 'A.002'"]}
```

## See Also

- [seqret](seqret.md) - read and write sequences
- [maskfeat](maskfeat.md) - mask the regions of features
- [descseq](descseq.md) - alter the name or description of a sequence
