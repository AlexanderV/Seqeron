# run_mummer1

MUMmer run-mummer1 script: mummer -mum -l 20 -r | gaps | mummer-annotate; returns .out, .gaps, .align, .errorsgaps.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `run_mummer1` |
| **Method ID** | `MummerRunPrograms.RunMummer1` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 run-mummer1, the 'severely antiquated' script that runs the mummer 1 suite: mummer -mum -l 20 -r (reverse-complement matches only; the first output line is dropped) | gaps | mummer-annotate. reference and query are single-sequence FASTA texts; referenceName is the file name printed in the gaps file. Returns the files of the script: out (matches), gaps, align (the gaps with their alignments) and errorsGaps (the gaps with the number of errors). Without reverse (the script's fourth argument -r) the script passes an empty argument to mummer, which fails, and the tool refuses; a query with a second record makes gaps loop forever (refused), no match makes gaps crash (refused). Byte for byte equal to the installed script on 55 + runs.

## Core Documentation Reference

- Source: [MummerRunPrograms.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerRunPrograms.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_mummer_gaps_oracle.py](../../../../scripts/alignment/generate_mummer_gaps_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference FASTA text (one sequence). |
| `query` | string | Yes | Query FASTA text (one sequence). |
| `referenceName` | string | No | The reference file name as the script receives it (printed in the gaps file). (default: reference.fa) |
| `prefix` | string | No | The prefix of the output files (names the files in the result). (default: out) |
| `reverse` | boolean | No | The fourth argument -r (default true; without it the script fails). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `out` | string | The matches (prefix.out). |
| `gaps` | string | The gaps file (prefix.gaps). |
| `align` | string | The alignments of the gaps (prefix.align). |
| `errorsGaps` | string | The gaps file with error counts (prefix.errorsgaps). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | /usr/bin/mummer: cannot open file "" or file "" is empty (without -r ...) |
| 1001 | gaps: no matches (the program dereferences a null array and crashes). |

## Examples

### Example 1: Reverse matches of a query

**User Prompt:**
> Run the mummer 1 suite.

**Expected Tool Call:**
```json
{"tool": "run_mummer1", "arguments": {"reference": ">r\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n", "query": ">q\nacgt\n", "reverse": false}}
```

**Response:**
```json
{}
```

Without -r the script passes an empty argument to mummer, which fails; the tool refuses (an error).

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
