# show_diff

MUMmer 3.23 show-diff: breaks, gaps, jumps, inversions, duplications and sequence jumps of a delta file.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `show_diff` |
| **Method ID** | `MummerShowDiff.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 show-diff: the structural differences between the sequences of a delta file read off the many-to-many mapping (delta-filter -m) along every reference record (referenceDiff, -r, the default) or query record (queryDiff, -q; not both): BRK breaks (unaligned stretches), GAP (gap sizes in both sequences and their difference), JMP (a jump within a sequence), INV (an inversion), DUP (a duplication) and SEQ (a jump to another sequence, with the two neighbour ids); columns [SEQ] [TYPE] [S1] [E1] [LEN 1] (+ extra fields per type); amos (-f) writes the AMOS layout with the coordinates repeated. header = false is -H. The ties of the longest-increasing-subset choice use the program's rand() (seed 1). Byte for byte equal to the installed show-diff on 300 + runs.

## Core Documentation Reference

- Source: [MummerShowDiff.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerShowDiff.cs)
- Algorithm doc: [Nucmer_And_Delta_Programs.md](../../../algorithms/Alignment/Nucmer_And_Delta_Programs.md)
- Oracle: [generate_mummer_show_oracle.py](../../../../scripts/alignment/generate_mummer_show_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `delta` | string | Yes | Delta file text. |
| `amos` | boolean | No | -f: AMOS output (default false). |
| `referenceDiff` | boolean | No | -r: differences of the reference (default; the default when neither -r nor -q is set). |
| `queryDiff` | boolean | No | -q: differences of the query (cannot be combined with -r). (default: false) |
| `header` | boolean | No | Print the header (default true; false = -H). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What show-diff writes to standard output. |
| `messages` | string | What it writes to standard error. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Delta file text cannot be null or empty |
| 1001 | ERROR: -r and -q options cannot be combined |

## Examples

### Example 1: A jump between two query contigs

**User Prompt:**
> Classify the differences along the reference.

**Expected Tool Call:**
```json
{"tool": "show_diff", "arguments": {"delta": "ref2.fa qry2.fa\nNUCMER\n>r q1 600 320\n1 320 1 320 0 0 0\n0\n>r q2 600 320\n281 600 320 1 0 0 0\n0\n", "referenceDiff": true}}
```

**Response:**
```json
{"text": "ref2.fa qry2.fa\nNUCMER\n\n[SEQ]\t[TYPE]\t[S1]\t[E1]\t[LEN 1]\nr\tSEQ\t321\t280\t-40\tq1\tq2\n", "messages": ""}
```

show-diff -r d2.delta: the reference is covered by q1 (1-320) and by the inverted q2 (281-600); the end of q1 and the start of q2 overlap by 40 bases, reported as a SEQ jump of length -40 between the two contigs.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
