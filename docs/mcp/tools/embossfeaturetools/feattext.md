# feattext

placeholder

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `feattext` |
| **Method ID** | `EmbossFeaturePrograms.Feattext` |
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
| `files` | object | No | More files by name (file name -> text) that featuresValue may name, list files included ('@list'). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text of the entries the readers stored. |
| `tables` | integer | The number of feature tables read. |
| `messages` | array | What the library printed on its error stream (Warning: ..., Error: ...), in order. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty feature file, a file that no feature reader reads - 'Bad value for -features and no more retries' -, an output format this port does not write for a feature table, an input on which the reference program never stops). |

## Examples

### Example 1: The original text of an EMBL entry

**User Prompt:**
> The original text of an EMBL entry

**Expected Tool Call:**
```json
{"tool": "feattext", "arguments": {"features": "ID   SEQA; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   gene            1..30\nXX\n//\n"}}
```

**Response:**
```json
{"text": "ID   SEQA; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   gene            1..30\nXX\n//\n", "tables": 1, "messages": []}
```

## See Also

- [featcopy](featcopy.md) - read and write feature tables
- [featmerge](featmerge.md) - combine two feature tables
- [wordmatch](wordmatch.md) - exact word matches written as feature files
