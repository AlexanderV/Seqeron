# delta_filter

MUMmer 3.23 delta-filter: keeps the alignments of a delta file that survive -i -l -u -q -r -g -m -1.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `delta_filter` |
| **Method ID** | `MummerDeltaFilter.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 delta-filter: reads a delta file (nucmer / promer) and keeps the alignments that survive the switches, in the order -i -l -u -q -r -g -m -1; the output is a delta file again (grouped by query record). minIdentity (-i, 0 - 100) and minLength (-l) drop weak alignments, minUnique (-u) those whose reference AND query part is not unique; queryLis (-q) maps each query position to its best hit in the reference (allowing reference overlaps), referenceLis (-r) each reference position to its best hit, globalLis (-g) the longest set of mutually consistent alignments (no rearrangements), manyToMany (-m) the union and oneToOne (-1) the intersection of the -r and -q chains; maxOverlap (-o) is the largest overlap of two chained alignments as a percent of their length; epsilon (-e, default -1) is the negligible score loss: a negative value breaks ties among equally good chains with glibc rand() seeded with 1, exactly as the program (reproducible). Byte for byte equal to the installed delta-filter on 300 + runs of random option sets over delta files of rearranged, duplicated and inverted multi-contig inputs. A percentage outside [0, 100] or a negative length is an error (the program's range errors).

## Core Documentation Reference

- Source: [MummerDeltaFilter.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerDeltaFilter.cs)
- Algorithm doc: [Nucmer_And_Delta_Programs.md](../../../algorithms/Alignment/Nucmer_And_Delta_Programs.md)
- Oracle: [generate_mummer_show_oracle.py](../../../../scripts/alignment/generate_mummer_show_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The filtered delta file text, grouped by query record. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Delta file text cannot be null or empty |
| 1001 | ERROR: Minimum identity must be within the range [0, 100] (and the other range errors) |

## Examples

### Example 1: 1-to-1 filter of one alignment

**User Prompt:**
> Keep only the 1-to-1 alignments of this delta file.

**Expected Tool Call:**
```json
{"tool": "delta_filter", "arguments": {"delta": "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n", "oneToOne": true}}
```

**Response:**
```json
{"text": "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n"}
```

delta-filter -1 d1.delta: the single alignment is its own best hit on both sides and is kept.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
