# mash_screen

Mash 2.3 'mash screen' for nucleotide queries, byte for byte.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mash_screen` |
| **Method ID** | `MashCommands.Screen` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Mash 2.3 'mash screen'.

## Core Documentation Reference

- Algorithm: [Mash_Commands.md](../../algorithms/Alignment/Mash_Commands.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `queries` | array | Yes | The query files as FASTA texts (records start with '>': the first word of the header is the sequence id, the rest the comment). |
| `mixtures` | array | Yes | The mixture files as FASTA texts (reads or contigs). |
| `queryNames` | array<string> | No | Names of the query files as the program prints them (default query1, ...). |
| `k` | integer | No | k-mer size of the query sketch (-k of mash sketch, 1 - 32; default 21, or 9 with protein). |
| `sketchSize` | integer | No | Sketch size of the query sketch (-s, default 1000). (default: 1000) |
| `seed` | integer | No | Hash seed of the query sketch (-S, 0 - 4294967295, default 42). (default: 42) |
| `noncanonical` | boolean | No | The query sketch was made with -n (strand preserved). (default: false) |
| `individual` | boolean | No | The query sketch was made with -i (one reference per sequence). (default: false) |
| `winnerTakesAll` | boolean | No | -w: winner-takes-all strategy for identity estimates. (default: false) |
| `minIdentity` | number | No | -i: minimum identity to report (-1 - 1, default 0). (default: 0) |
| `maxPValue` | number | No | -v: maximum p-value to report (0 - 1, default 1). (default: 1) |
| `protein` | boolean | No | -a: the amino-acid alphabet (A-Z except BJOUXZ; implies -n and, without k, -k 9): k-mers of the 20 letters, 64-bit hashes when 20^k > 2^32, p-values over the k-mer space 20^k. (default: false) |
| `alphabet` | string | No | -z: the alphabet the hashes are based on (case ignored unless preserveCase; k-mers with other characters are ignored; implies -n; ASCII only: the program has no defined result for a byte of 0x80 or more). protein takes precedence. |
| `preserveCase` | boolean | No | -Z: preserve the case of the letters in the k-mers and the alphabet (a letter whose case is not in the alphabet is skipped). (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text mash screen writes. |

## Errors

| Code | Message |
|------|---------|
| 1001 | ERROR: Did not find sequence records in inputs |
| 1001 | -v must be in 0..1 (mash) |

## Examples

### Example 1: a query contained in a mixture of reads

**User Prompt:**
> Screen a 40-base genome against reads that contain it.

**Expected Tool Call:**
```json
{"tool": "mash_screen", "arguments": {"queries": [">g\nACGTTGCATGCCGATAGCTAGGATCCATGCATGCAAGT\n"], "mixtures": [">r1\nACGTTGCATGCCGATAGCTAGGATCC\n>r2\nGATAGCTAGGATCCATGCATGCAAGT\n"], "queryNames": ["g.fa"], "k": 8, "sketchSize": 100}}
```

**Response:**
```json
{"text": "1\t30/30\t1\t6.59538e-101\tg.fa\tg \n"}
```

## See Also

- [mash_dist](mash_dist.md) - Distances of queries to references
