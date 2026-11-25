# dbtell

EMBOSS 6.6.0 dbtell: the definition of a database of the standard EMBOSS configuration.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `dbtell` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 dbtell: the definition of the database as the standard configuration of the installation holds it (the ontologies chebi, eco, edam, edam_data, edam_format, edam_identifier, edam_operation, edam_topic, go, go_component, go_function, go_process, pw, ro, so, swo, the data resource catalogue drcat and the taxonomy taxon): the file that defines it, the access levels (id, query, all), the attributes; an alias (any case other than the lower case name) is reported as such; a name that only the DRCAT catalogue defines (embl, uniprot, pdb ...) gives the line 'defined in DRCAT'; any other name gives no output. The databases of a user configuration file (.embossrc) and of a server (-server) are not read. Equals the installed program on 302 runs (standard output, error stream, exit status and written files byte for byte). B09 WP67, F309.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Dbtell.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.SequenceFormats/EmbossSequencePrograms.Dbtell.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_db_oracle.py](../../../../scripts/alignment/generate_emboss_prog_db_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `database` | string | Yes | The name of the database (any case), e.g. go_process, drcat, EMBL. |
| `full` | boolean | No | Verbose output (-full): the aliases of the database follow its definition. (default: false) |

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

## Examples

### Example 1: A run on the standard configuration

**User Prompt:**
> Tell me the definition of the database go_process.

**Expected Tool Call:**
```json
{"tool": "dbtell", "arguments": {"database": "go_process"}}
```

**Response:**
```json
{"text": "# go_process is defined in standard\n# access levels id: Yes query: Yes all: Yes\n\nDBNAME go_process [\n   format:         \"obo\"\n   method:         \"emboss\"\n   type:           \"obo\"\n   comment:        \"Gene Ontology (biological processes)\"\n   directory:      \"/usr/share/EMBOSS/data/OBO\"\n   field:          \"id ! identifier without the prefix\"\n   field:          \"acc ! full name and any alternate identifier(s)\"\n   field:          \"nam ! words in the name\"\n   field:          \"isa ! parent identifier from is_a relation(s)\"\n   field:          \"des ! words in the description\"\n   field:          \"ns ! namespace\"\n   indexdirectory: \"/usr/share/EMBOSS/index\"\n   namespace:      \"biological_process\"\n   release:        \"1.1.3300\"\n   edamdat:        \"0966 ! Ontology term\"\n   edamfmt:        \"2196 ! OBO format\"\n   edamid:         \"1176 ! GO term ID\"\n   edamtpc:        \"0089 ! Ontologies\"\n   taxon:          \"1 all\"\n]\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
