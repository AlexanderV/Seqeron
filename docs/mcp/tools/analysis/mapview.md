# mapview

MUMmer mapview: FIG maps of show-coords -rl alignments (percent identity along the reference, links between contig matches, genes), byte for byte.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mapview` |
| **Method ID** | `MummerMapview.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 mapview (perl script): turns the output of show-coords -rl (nucmer or promer; or of mgaps with mgapsReference = -g) into FIG 3.2 files that show the alignments along the reference - the percent identity (50 - 100) on the vertical axis, a red match line per alignment, a coloured link line between the alignments of one query contig (one of three colours; no link when the next match of the contig is more than maxDistance = -d (default 50 000) bases away) - partitioned into numberOfFiles = -n (default 10) files of an equal share of the reference (a multi-reference input gets one file per group of references that fit a share), named <prefix>_<n>.fig (prefix = -p; default 'NUCMER_graph' / 'PROMER_graph': the program name in the second line of the coords file), with a kilobase scale, the percent identity ticks and a legend. With a UTR file and a CDS file (GFF: sequence id, source, exon type single-exon / initial-exon / internal-exon / last-exon, start, end, score, strand, frame, gene name) the genes are drawn above the map (5' UTR, CDS, 3' UTR, introns; overlapping genes on levels), x1 / x2 (-x1 / -x2) bound the displayed range, printQueryNames (-I) / printGeneNames (-Ir) write the names. format (-f) 'pdf' or 'ps' makes the script convert every FIG file with fig2dev: the result then lists the tmp<prefix>_<n>.fig texts the script feeds to it (figs), the fig2dev command lines (commands) and the names of the converted files (outputs) - the conversion itself is done by fig2dev, which is not part of this library. verbose (-v) lists the files and their ranges on standard output (stdout). A UTR file without a CDS file, or GFF files whose first column names no reference of the coords file, stop the script (exitCode 1, the message on messages). An mgaps input (-g) also writes the file <coordsName>coords (side): the clusters in the format of show-coords, sorted by position. The messages of a wrong option (x1 > x2, a negative bound) are printed and the script goes on, as the program does. Byte for byte equal to the installed mapview on 263 + runs (nucmer and promer coords, one or several references, mgaps input, GFF genes on both strands, every option) and 800 more random runs.

## Core Documentation Reference

- Source: [MummerMapview.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerMapview.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_mapview_oracle.py](../../../../scripts/alignment/generate_mapview_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coords` | string | Yes | The coords file (show-coords -rl output, or mgaps output with mgapsReference). |
| `utrGff` | string | No | UTR coordinates in GFF format (needs the CDS file too). |
| `cdsGff` | string | No | CDS coordinates in GFF format. |
| `maxDistance` | integer | No | -d: the largest distance in base pairs between linked matches (default 50 000). |
| `format` | string | No | -f: output format 'fig' (default), 'pdf' or 'ps' (the last two through fig2dev: see the description). |
| `magnification` | number | No | -m: magnification handed to fig2dev (default none). |
| `numberOfFiles` | integer | No | -n: number of output files the reference is partitioned into (default 10). |
| `prefix` | string | No | -p: prefix of the output files (default 'NUCMER_graph' / 'PROMER_graph' by the coords file). |
| `verbose` | boolean | No | -v: list the output files and their ranges (default false). |
| `x1` | integer | No | -x1: lower coordinate bound of the display (default none). |
| `x2` | integer | No | -x2: upper coordinate bound of the display (default none). |
| `mgapsReference` | string | No | -g: the input is mgaps output; the reference sequence id (as in the first column of the GFF files). (default: ) |
| `printQueryNames` | boolean | No | -I: display the names of the query sequences (default false). |
| `printGeneNames` | boolean | No | -Ir: display the names of the reference genes (default false). |
| `coordsName` | string | No | The file name of the coords file as the script receives it (the mgaps side file is <name>coords; default coords). (default: coords) |
| `utrName` | string | No | The file name of the UTR file; when it equals cdsName the exon lists are not split (default utr.gff). |
| `cdsName` | string | No | The file name of the CDS file (default cds.gff). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `figs` | object | The FIG files by name (with pdf / ps the temporary tmp<prefix>_<n>.fig files the script feeds to fig2dev). |
| `commands` | array | The fig2dev command lines (pdf / ps). |
| `outputs` | array | The names of the files fig2dev writes. |
| `stdout` | string | Standard output of the script. |
| `messages` | string | Standard error of the script. |
| `exitCode` | integer | The exit status. |
| `side` | object | The side file of an mgaps input (<coordsName>coords). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coords text cannot be null or empty |
| 1001 | format must be 'fig', 'pdf' or 'ps' |

## Examples

### Example 1: One alignment on one reference

**User Prompt:**
> Draw the alignments of this coords file into one FIG file.

**Expected Tool Call:**
```json
{"tool": "mapview", "arguments": {"coords": "ref.fa qry.fa\nNUCMER\n\n    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  | [TAGS]\n==========================================================================================================\n      11      250  |        1      240  |      240      240  |    98.76  |      260      240  | ref\tqry\n", "numberOfFiles": 1, "prefix": "demo"}}
```

**Response:**
```json
{"figs": {"demo_0.fig": "#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nMultiple\n-2\n1200 2\n2 1 0 10 11 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 0 3050 17 3050\n2 1 2 1 0 7 60 0 -1 4.000 0 0 -1 0 0 2\n\t0 5250 17 5250\n2 1 2 1 0 7 60 0 -1 4.000 0 0 -1 0 0 2\n\t0 4250 17 4250\n2 1 2 1 0 7 60 0 -1 4.000 0 0 -1 0 0 2\n\t0 3250 17 3250\n2 1 0 1 0 7 60 0 -1 0.000 0 0 -1 0 0 2\n\t0 5250 0 5300\n4 0 0 100 0 0 8 0.0000 4 135 405 \t 0 5400 0k\\001\n4 0 0 50 0 0 14 0.0000 4 135 450 0 3000 0\\001\n4 0 0 50 0 0 14 0.0000 4 135 810 -233 3000 259\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t-200 5270 50%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t37 5270 50%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t-200 4270 75%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t37 4270 75%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t-200 3270 100%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t37 3270 100%\\001\n2 1 0 2 4 0 40 0 -1 0.000 0 0 -1 0 0 2\n\t0 3330 17 3330\n2 1 0 5 20 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 0 3150 17 3150\n4 0 0 100 0 0 8 0.0000 4 135 405 \t100 6100  Legend \\001\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 70 6300 99 6300\n2 1 0 5 27 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 100 6300 200 6300\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 200 6300 230 6300\n4 0 0 100 0 0 8 0.0000 4 135 405 \t300 6330  5' utr \\001\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 70 6450 99 6450\n2 1 0 5 2 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 100 6450 200 6450\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 200 6450 230 6450\n4 0 0 100 0 0 8 0.0000 4 135 405 \t300 6480  cds \\001\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 70 6600 99 6600\n2 1 0 5 1 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 100 6600 200 6600\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 200 6600 230 6600\n4 0 0 100 0 0 8 0.0000 4 135 405 \t300 6630  3' utr \\001\n2 1 0 2 4 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t100 6750 200 6750\n4 0 0 100 0 0 8 0.0000 4 135 405 \t300 6780  match found by NUCMER \\001\n"}, "commands": [], "outputs": [], "stdout": "", "messages": "", "exitCode": 0, "side": {}}
```

mapview -n 1 -p demo x.coords: the reference line, the grid at 50 / 75 / 100 percent identity, the red match line at 98 percent, the kilobase scale and the legend.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
