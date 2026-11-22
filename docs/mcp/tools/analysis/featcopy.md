# featcopy

placeholder

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `featcopy` |
| **Method ID** | `EmbossFeaturePrograms.Featcopy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

placeholder

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `features` | string | Yes | The text of the feature file (a feature table of any input format, or an EMBL / GenBank entry; several tables may follow one another). |
| `featuresValue` | string | No | The value of -features: 'format::' prefix and [begin:end] suffix are read as the program reads them; default the file itself. |
| `outputFormat` | string | No | The output feature format (-offormat): gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug. |
| `featureBegin` | integer | No | -fbegin: the first table is cut to this range start (default 1; negative counts from the end). |
| `featureEnd` | integer | No | -fend: the first table is cut to this range end (default: its length). |
| `files` | object | No | More files by name (file name -> text) that featuresValue may name: a list file ('@list' or 'list::list': one query per word of a line, '#' starts a comment, the queries are files with an optional 'format::' prefix and [begin:end] suffix, or lists again), the files the list names; a file that is not here does not exist. 'asis::' is refused for features as the program refuses it. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fileName` | string | The name of the output file (the value of -outfeat without a format:: prefix). |
| `text` | string | The feature tables in the output format. |
| `format` | string | The output feature format as the program resolved it. |
| `tables` | integer | The number of feature tables written. |
| `messages` | array | What the library printed on its error stream (Warning: ..., Error: ...), in order. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty feature file, a file that no feature reader reads - 'Bad value for -features and no more retries' -, an output format this port does not write for a feature table, an input on which the reference program never stops). |

## Examples

### Example 1: Write a GFF3 feature table as EMBL

**User Prompt:**
> Write a GFF3 feature table as EMBL

**Expected Tool Call:**
```json
{"tool": "featcopy", "arguments": {"features": "##gff-version 3\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1;Name=abc\n", "outputFormat": "embl"}}
```

**Response:**
```json
{"fileName": "unknown.gff", "text": "FH   Key             Location/Qualifiers\nFH\nFT   gene            10..50\nFT                   /note=\"*Name abc\"\n", "format": "embl", "tables": 1, "messages": []}
```

## See Also

- [diffseq](diffseq.md) - differences between two sequences (reads features from separate files)
- [wordmatch](wordmatch.md) - exact word matches written as feature files
