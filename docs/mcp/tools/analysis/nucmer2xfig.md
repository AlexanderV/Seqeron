# nucmer2xfig

MUMmer nucmer2xfig: an xfig plot of a show-coords -c -l file (contigs against a reference chromosome).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `nucmer2xfig` |
| **Method ID** | `MummerNucmer2Xfig.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

nucmer2xfig of MUMmer 3.23 (perl script, Steven Salzberg 2001): an xfig file (FIG 3.2) plotting a coords file of 'show-coords -c -l' (19 white-space separated fields per alignment: s1 e1 | s2 e2 | len1 len2 | %idy | lenR lenQ | covR covQ | idR idQ; no header lines, sorted by query id and start) of the contigs of one genome against a reference chromosome: the reference as a red line (colour 4) with its id, x-axis tick marks every 250 000 bases and y-axis tick marks for 50 - 100 percent identity, every alignment as a green line (colour 2, x scaled by 0.005, y = 20 x (identity - 50)), every contig - the stretch its alignments span, split where they are more than twice the query length apart - as a blue line (colour 1) at 100 percent with a vertical label (4.7124 rad) when the line is wider than the label. Perl's numeric conversion of the fields is reproduced (text is 0) and printf %.0f rounds half to even. Byte for byte equal to the installed script on 260 + runs including malformed input (the output of an empty file is just the reference line with empty fields).

## Core Documentation Reference

- Source: [MummerNucmer2Xfig.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerNucmer2Xfig.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_nucmer2xfig_oracle.py](../../../../scripts/alignment/generate_nucmer2xfig_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coords` | string | Yes | The text of the coords file (show-coords -c -l without its header). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The xfig file nucmer2xfig writes to standard output. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coords text cannot be null |

## Examples

### Example 1: Two alignments of one contig

**User Prompt:**
> Plot this coords file as xfig.

**Expected Tool Call:**
```json
{"tool": "nucmer2xfig", "arguments": {"coords": "  2551   2577 | 240 266 | 27 27 | 96.30 | 30000 1424 | 0.09 1.90 | 2R ctg1\n  2600   9000 | 300 6700 | 6401 6401 | 99.00 | 30000 1424 | 21.33 99.00 | 2R ctg1\n"}}
```

**Response:**
```json
{"text": "#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nSingle\n-2\n1200 2\n2 1 0 2 2 7 100 0 -1 0.000 0 0 -1 0 0 2\n\t 12 926 13 926\n2 1 0 2 2 7 100 0 -1 0.000 0 0 -1 0 0 2\n\t 13 980 45 980\n2 1 0 2 1 7 100 0 -1 0.000 0 0 -1 0 0 2\n\t 13 1000 45 1000\n4 0 0 100 0 0 12 0.0000 4 135 405 150 0 2R\\001\n2 1 0 2 4 7 100 0 -1 0.000 0 0 -1 0 0 2\n\t 0 0 150 0\n4 0 0 100 0 0 5 0.0000 4 135 405 -150 10 50\\001\n2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n-50 0 0 0\n4 0 0 100 0 0 5 0.0000 4 135 405 -150 210 60\\001\n2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n-50 200 0 200\n4 0 0 100 0 0 5 0.0000 4 135 405 -150 410 70\\001\n2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n-50 400 0 400\n4 0 0 100 0 0 5 0.0000 4 135 405 -150 610 80\\001\n2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n-50 600 0 600\n4 0 0 100 0 0 5 0.0000 4 135 405 -150 810 90\\001\n2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n-50 800 0 800\n4 0 0 100 0 0 5 0.0000 4 135 405 -150 1010 100\\001\n2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n-50 1000 0 1000\n"}
```

nucmer2xfig c.coords: the two alignments are green lines (identity 96.30 and 99.00), the contig a blue line at 100 percent, the reference a red line with percent identity tick marks.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
