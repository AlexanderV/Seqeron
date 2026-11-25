# showdb

EMBOSS 6.6.0 showdb: the table of the databases of the standard EMBOSS configuration.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `showdb` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 showdb: the table of the databases of the standard configuration of the installation (taxon, drcat and the ontologies), grouped by type and sorted by definition file and name, with the columns selected by qualifiers (-full all columns; -only with -type -id -query -all -comment -methods -taxscope -fields -numfields -aliases -numaliases -examples -numexamples -defined -release -heading -access to choose the columns; -noprotein -nonucleic -nosequence -nofeature -notext -notaxonomy -noresource -noobo -noxml and -nostandard -nouser -noinclude to leave databases out; -html for an HTML table); with database only that database (a name the configuration does not define stops the program). The databases of a user configuration file (.embossrc) are not read. Equals the installed program on 301 runs (standard output, error stream, exit status and written files byte for byte). B09 WP67, F309.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Dbtell.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.SequenceFormats/EmbossSequencePrograms.Dbtell.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_db_oracle.py](../../../../scripts/alignment/generate_emboss_prog_db_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `database` | string | No | The name of one database to give information on (-database); default: all. |
| `qualifiers` | string | No | Qualifiers as typed on the command line, e.g. "-full -nodefined" or "-only -type -comment -noheading". |

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
> Show the database taxon.

**Expected Tool Call:**
```json
{"tool": "showdb", "arguments": {"database": "taxon"}}
```

**Response:**
```json
{"text": "# Name         Type     Comment\n# ============ ======== =======\ntaxon          Taxonomy NCBI taxonomy\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
