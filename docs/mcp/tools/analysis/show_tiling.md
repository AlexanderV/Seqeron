# show_tiling

MUMmer 3.23 show-tiling: places the aligned query contigs on the reference and prints the tiling path.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `show_tiling` |
| **Method ID** | `MummerShowTiling.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 show-tiling: places the aligned query contigs of a delta file on the reference. For every contig and reference the longest consistent subset of its alignments (maxGap -g, default 1000 for nucmer files / -1 for promer) is taken; a contig is tiled when its coverage is at least minCoverage (-v, default 95) and better than the second best by minCoverageDifference (-V, default 10) and its identity at least minIdentity (-i, default 90), or, when it hangs over the reference ends, by the same rule on the overhang-corrected coverage; contigs shorter than minContigLength (-l) are ignored. The tiled contigs are linked along each reference (circular = -c lets them span the origin) and printed as '>refId N bases' followed by 'start end gapToNext length coverage identity strand id' (the tiling path), or with printAlignments (-a) as the tab-delimited alignment rows, or with xml (-x) as the XML link evidence. contigList (-t), unusable (-u) and pseudoMolecule (-p, needs the query FASTA) return the program's files. randomRepeats (-R) places repetitive contigs at a random copy (rand() seeded by randomSeed; the program uses the current time; implies minCoverageDifference 0). Byte for byte equal to the installed show-tiling on 350 + runs.

## Core Documentation Reference

- Source: [MummerShowTiling.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerShowTiling.cs)
- Algorithm doc: [Nucmer_And_Delta_Programs.md](../../../algorithms/Alignment/Nucmer_And_Delta_Programs.md)
- Oracle: [generate_mummer_show_oracle.py](../../../../scripts/alignment/generate_mummer_show_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The tiling path (or the alignment rows / XML) written to standard output. |
| `contigList` | string | The -t file (when requested). |
| `unusable` | string | The -u file (when requested). |
| `pseudoMolecule` | string | The -p file (when requested). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Delta file text cannot be null or empty |
| 1001 | ERROR: Percentages must be within the range [0.0, 100.0] |

## Examples

### Example 1: Two contigs tiled on one reference

**User Prompt:**
> Tile the query contigs on the reference.

**Expected Tool Call:**
```json
{"tool": "show_tiling", "arguments": {"delta": "ref2.fa qry2.fa\nNUCMER\n>r q1 600 320\n1 320 1 320 0 0 0\n0\n>r q2 600 320\n281 600 320 1 0 0 0\n0\n"}}
```

**Response:**
```json
{"text": ">r 600 bases\n1\t320\t-40\t320\t100.00\t100.00\t+\tq1\n281\t600\t0\t320\t100.00\t100.00\t-\tq2\n", "contigList": null, "unusable": null, "pseudoMolecule": null}
```

show-tiling d2.delta: q1 forward at 1-320, q2 reverse at 281-600 (gap -40 = overlap).

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
