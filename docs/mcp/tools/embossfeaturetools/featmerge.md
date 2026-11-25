# featmerge

placeholder

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `featmerge` |
| **Method ID** | `EmbossFeaturePrograms.Featmerge` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

placeholder

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `featuresA` | string | Yes | The text of the first feature file (-afeatures). |
| `featuresB` | string | Yes | The text of the second feature file (-bfeatures). |
| `featuresAValue` | string | No | The value of -afeatures: 'format::' prefix and [begin:end] suffix; default the file itself. |
| `featuresBValue` | string | No | The value of -bfeatures. |
| `outputFormat` | string | No | The output feature format (-offormat): gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug. |
| `featureBegin` | integer | No | -fbegin1: the range start of the first table (default 1). |
| `featureEnd` | integer | No | -fend1: the range end of the first table (default: its length). |
| `files` | object | No | More files by name (file name -> text) that featuresAValue / featuresBValue may name, list files included ('@list'). |

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

### Example 1: Merge two GFF3 tables and write BED

**User Prompt:**
> Merge two GFF3 tables and write BED

**Expected Tool Call:**
```json
{"tool": "featmerge", "arguments": {"featuresA": "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1;Name=abc\n", "featuresB": "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\texon\t20\t30\t.\t+\t.\tID=e1\n", "outputFormat": "bed"}}
```

**Response:**
```json
{"fileName": "unknown.gff", "text": "track name=seq1 description=\"EMBOSS:6.6.0.0\" useScore=0\nseq1\t9\t49\tgene\t0\t+\t9\t49\t0\t1\t0\t41\nseq1\t19\t29\texon\t0\t+\t19\t29\t0\t1\t0\t11\n", "format": "bed", "tables": 1, "messages": []}
```

## See Also

- [featcopy](featcopy.md) - read and write feature tables
- [wordmatch](wordmatch.md) - exact word matches written as feature files
