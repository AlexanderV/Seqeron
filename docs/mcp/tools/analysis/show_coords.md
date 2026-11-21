# show_coords

MUMmer 3.23 show-coords: the alignment coordinates of a delta file as a table, tab-delimited or btab.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `show_coords` |
| **Method ID** | `MummerShowCoords.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 show-coords: the coordinates of the alignments of a delta file as the human-readable table (header, [S1] [E1] | [S2] [E2] | [LEN 1] [LEN 2] | [% IDY] | tags), the tab-delimited table (tabular, -T) or btab (btab, -B). header = false is -H; brief (-b) merges overlapping alignments regardless of direction and drops the identity columns; coverage (-c) adds the percent coverage columns, lengths (-l) the sequence lengths, direction (-d) the FRM columns; minIdentity (-I) and minLength (-L) filter; sortByReference (-r) / sortByQuery (-q) sort; annotateOverlaps (-o) tags alignments that cover a whole sequence end ([IDENTITY] [CONTAINED] [CONTAINS] [END] [BEGIN]); the hidden warnings (-w) tag [DUPLICATE] [CONTAINS] [SHADOWED] [OVERLAPS]; las (-g) / weightedLas (-G) keep the longest ascending subset (ties broken with rand() seeded by randomSeed, the program seeds it with the current time); knockout (-k) is for promer files. Percentages are printed from the float values with %.2f. Byte for byte equal to the installed show-coords on 450 + runs; PROmer delta files are accepted.

## Core Documentation Reference

- Source: [MummerShowCoords.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerShowCoords.cs)
- Algorithm doc: [Nucmer_And_Delta_Programs.md](../../../algorithms/Alignment/Nucmer_And_Delta_Programs.md)
- Oracle: [generate_mummer_show_oracle.py](../../../../scripts/alignment/generate_mummer_show_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `delta` | string | Yes | Delta file text. |
| `header` | boolean | No | Print the header (default true; false = -H). |
| `brief` | boolean | No | -b: merge overlapping alignments, no identity columns (default false). |
| `btab` | boolean | No | -B: btab output (default false). |
| `coverage` | boolean | No | -c: percent coverage columns (default false). |
| `direction` | boolean | No | -d: alignment direction in the FRM columns (default false; always for promer). |
| `minIdentity` | number | No | -I: minimum alignment identity, 0 - 100 (default 0). |
| `knockout` | boolean | No | -k: knockout (promer data only, default false). (default: false) |
| `lengths` | boolean | No | -l: sequence length columns (default false). |
| `minLength` | integer | No | -L: minimum alignment length (default 0). |
| `annotateOverlaps` | boolean | No | -o: annotate overlaps (default false). |
| `sortByQuery` | boolean | No | -q: sort by query (default false). |
| `sortByReference` | boolean | No | -r: sort by reference (default false). |
| `las` | boolean | No | -g: longest ascending subset, deprecated (default false). |
| `weightedLas` | boolean | No | -G: weighted longest ascending subset, hidden (default false). |
| `tabular` | boolean | No | -T: tab-delimited output (default false). |
| `warnings` | boolean | No | -w: warnings column, hidden (default false). |
| `randomSeed` | integer | No | Seed of the rand() that breaks ties of -g / -G (the program uses the current time); default = the current time. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What show-coords writes to standard output. |
| `messages` | string | What it writes to standard error (warnings). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Delta file text cannot be null or empty |
| 1001 | ERROR: Percentages must be within the range [0.0. 100.0] |

## Examples

### Example 1: Coordinates of three substitutions

**User Prompt:**
> Show the coordinates of the alignments of this delta file.

**Expected Tool Call:**
```json
{"tool": "show_coords", "arguments": {"delta": "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n"}}
```

**Response:**
```json
{"text": "ref.fa qry.fa\nNUCMER\n\n    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]\n=====================================================================================\n       6      140  |        1      135  |      135      135  |    97.78  | ref\tqry\n", "messages": ""}
```

show-coords d1.delta: 135 of 135 bases aligned, 97.78 % identity.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
