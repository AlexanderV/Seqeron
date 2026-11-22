# marscan

EMBOSS 6.6.0 marscan: Finds matrix/scaffold recognition (MRS) signatures in DNA sequences, as marscan of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `marscan` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 marscan: the MAR/SAR recognition signatures: the 16 bp site awwrtaannwwgnnnc (one mismatch allowed) or its reverse gnnncwwnnttaywwt and the 8 bp site aataayaa or ttrttatt within 200 bp of each other; of a cluster of signatures the pair with the least distance is reported as the report out.txt (default format gff: the feature misc_signal with the positions of the two sites; every report format of EMBOSS, the report qualifiers of reportQualifiers). Equals the installed program on 300 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Marscan.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Marscan.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_wp63_oracle.py](../../../../scripts/alignment/generate_emboss_prog_wp63_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text of the main output of the program (the report, or its standard output). |
| `files` | object | The other files the program wrote (sequences, features, graph data), by name. |
| `messages` | string[] | What the library printed on the error stream (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequence cannot be empty |
| 1001 | The program died: the messages of the library (a value it refuses, an input no reader reads, a sequence of the wrong kind) |
| 1001 | An input on which the reference program never stops or dies by a fault is refused |

## Examples

### Example 1: A run on a small input

**User Prompt:**
> Run marscan on this input.

**Expected Tool Call:**
```json
{"tool": "marscan", "arguments": {"sequence": ">s1 test\nCCCCCCCCCCAATAACAACCCCCCCCCCCCCCAATGTAAGGTTGCCCCTTTTTTTTTT\n"}}
```

**Response:**
```json
{"text": "##gff-version 3\n##sequence-region s1 1 48\n#!Date 2023-11-14\n#!Type DNA\n#!Source-version EMBOSS 6.6.0.0\ns1\tmarscan\tregulatory_region\t11\t48\t.\t+\t.\tID=s1.1;note=*type MAR/SAR recognition site (MRS);note=*start8bp 11;note=*end8bp 18;note=*start16bp 33;note=*end16bp 48\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
