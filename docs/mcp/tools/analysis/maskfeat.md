# maskfeat

EMBOSS 6.6.0 maskfeat: the regions of the features of the wanted types are masked in the sequences of a file (any input format of seqret that carries features: embl, genbank, refseq, refseqp, swiss, nbrf, gff2, gff3; the feature types are those of the vocabulary, wildcards * and ? allowed, several types separated by blanks, commas, semicolons or bars; the default repeat* masks repeat_region and repeat_unit): the residues become the mask character (default N for nucleotide and X for protein sequences) or lower case (toLower, or a blank or no mask character), and the sequences are returned in FASTA (the sequence type of the program is any: every character that is no letter and no * is dropped). Equals the installed program on 1 200 runs (every input format with features, damaged files, the runs on which the program refuses an input or dies by a fault). A5-8c3g, B09 F249.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `maskfeat` |
| **Method ID** | `EmbossMaskfeatPrograms.Maskfeat` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 maskfeat: the regions of the features of the wanted types are masked in the sequences of a file (any input format of seqret that carries features: embl, genbank, refseq, refseqp, swiss, nbrf, gff2, gff3; the feature types are those of the vocabulary, wildcards * and ? allowed, several types separated by blanks, commas, semicolons or bars; the default repeat* masks repeat_region and repeat_unit): the residues become the mask character (default N for nucleotide and X for protein sequences) or lower case (toLower, or a blank or no mask character), and the sequences are returned in FASTA (the sequence type of the program is any: every character that is no letter and no * is dropped). Equals the installed program on 1 200 runs (every input format with features, damaged files, the runs on which the program refuses an input or dies by a fault). A5-8c3g, B09 F249.

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file with features (see the formats above), or a sequence as letters (nothing to mask). |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `type` | string | No | The feature types to mask (-type): names separated by blanks, commas, semicolons or bars; * and ? wildcards; default repeat*. (default: repeat*) |
| `toLower` | boolean | No | Make the masked regions lower case (-tolower). (default: false) |
| `maskChar` | string | No | The mask character (-maskchar, one character; none or a blank = lower case); default N (nucleotide) / X (protein). |
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
| InvalidParams | Invalid arguments (an empty sequence, a mask character of more than one character, a text no sequence reader reads, an input of an alignment or binary format this port does not read, an input on which the reference program dies by a fault). |

## Examples

### Example 1: The default type masks the repeats

**User Prompt:**
> Mask the repeat regions of this EMBL entry

**Expected Tool Call:**
```json
{"tool": "maskfeat", "arguments": {"sequence": "ID   SEQA; SV 1; linear; genomic DNA; STD; UNC; 20 BP.\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   repeat_region   3..9\nFT   gene            12..18\nXX\nSQ   Sequence 20 BP; 0 A; 0 C; 0 G; 0 T; 0 other;\n     aacgtacgta cgtacgtacg                                                20\n//\n", "sequenceFormat": "embl"}}
```

**Response:**
```json
{"sequenceText": ">SEQA SEQA.1\naaNNNNNNNacgtacgtacg\n", "sequences": 1, "messages": []}
```

## See Also

- [seqret](seqret.md) - read and write sequences
- [featcopy](featcopy.md) - read and write feature tables
