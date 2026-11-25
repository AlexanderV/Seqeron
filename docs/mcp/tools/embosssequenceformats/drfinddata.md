# drfinddata

EMBOSS 6.6.0 drfinddata: finds the data resources of the DRCAT catalogue by the EDAM data term of their output.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `drfinddata` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 drfinddata: the data resources of the catalogue DRCAT whose EDAM data terms match the keywords of query (comma separated; the identifier, the accession or a word of the name or a synonym of a term of the ontology edam; * is a wildcard): every term found and, with subclasses (default), all its descendants (is_a); every catalogue entry that lists such a term as EDAMdat is written once, in the order the indexes of the installation return them, in the outputFormat drcat (the entries as in DRCAT.dat), basic (identifier, name and URL) or list (drcat:IDENTIFIER). The EDAM ontology and the DRCAT catalogue are those of the EMBOSS 6.6.0 installation, embedded in the library. Equals the installed program on 233 runs (standard output, error stream, exit status and written files byte for byte); an output directory, or an output file in a directory that does not exist (the reference program dies by a segmentation fault), is refused. B09 WP67, F309.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Drfinddata.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.SequenceFormats/EmbossSequencePrograms.Drfinddata.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_db_oracle.py](../../../../scripts/alignment/generate_emboss_prog_db_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | Yes | The keywords, comma separated: an EDAM identifier (1300 or EDAM_data:1300), a word of a term name (transcript) or a pattern with * (Al*); upper or lower case. |
| `subclasses` | boolean | No | Also find the descendants of the terms (-subclasses, default true; false = -nosubclasses). (default: true) |
| `outputFormat` | string | No | The output format (-oformat): drcat (default), basic or list. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text the program writes to standard output. |
| `files` | object | The other files the program wrote, by name (none). |
| `messages` | string[] | What the library printed on the error stream (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The program died: the messages of the library (a value it refuses) |
| 1001 | An input on which the reference program dies by a fault is refused |

## Examples

### Example 1: A run on the standard configuration

**User Prompt:**
> Find the data resources whose output is a microarray data type.

**Expected Tool Call:**
```json
{"tool": "drfinddata", "arguments": {"query": "microarray", "subclasses": false, "outputFormat": "list"}}
```

**Response:**
```json
{"text": "drcat:CleanEx\ndrcat:ArrayExpress\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
