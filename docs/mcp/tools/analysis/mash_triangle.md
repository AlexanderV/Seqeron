# mash_triangle

Mash lower-triangular distance matrix or edge list (mash triangle).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mash_triangle` |
| **Method ID** | `MashCommands.Triangle` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Mash 2.3 'mash triangle': the distance of each input to every other input as a lower-triangular matrix in relaxed Phylip format (a first line with a TAB and the number of inputs, then one line per input: its name and the distances to the inputs before it), or with edgeList = true (-E; implied by a p-value or distance filter) the edge list 'name TAB name TAB distance TAB p-value TAB shared/total' of the pairs that pass maxPValue (-v) and maxDistance (-d). Sketches as mash_dist (a file without a record of at least k symbols is refused). VERIFICATION: 220 / 220 random sets (2 - 6 files, k 4 - 32, sketch size 5 - 2000, -n, -i, -E, -v, -d) equal the installed mash 2.3 byte for byte.

## Core Documentation Reference

- Algorithm: [Mash_Commands.md](../../algorithms/Alignment/Mash_Commands.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | The input files (FASTA texts), at least two for a matrix. |
| `names` | array<string> | No | Names of the inputs as the program prints them (default input1, ...). |
| `k` | integer | No | k-mer size (-k, 1 - 32, default 21). (default: 21) |
| `sketchSize` | integer | No | Sketch size (-s, default 1000). (default: 1000) |
| `seed` | integer | No | Hash seed (-S, 0 - 4294967295, default 42). (default: 42) |
| `noncanonical` | boolean | No | -n: preserve the strand (no canonical k-mers). (default: false) |
| `individual` | boolean | No | -i: sketch individual sequences instead of whole files. (default: false) |
| `edgeList` | boolean | No | -E: output the edge list instead of the matrix. (default: false) |
| `maxPValue` | number | No | -v: maximum p-value to report in the edge list (0 - 1, default 1). (default: 1.0) |
| `maxDistance` | number | No | -d: maximum distance to report in the edge list (0 - 1, default 1). (default: 1.0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text mash triangle writes. |

## Errors

| Code | Message |
|------|---------|
| 1001 | k must be in 1..32 (mash -k) |
| 1001 | seed must be in 0..4294967295 (mash -S) |
| 1001 | -v must be in 0..1 (mash) |

## Examples

### Example 1: Three genomes

**User Prompt:**
> Give the lower-triangular Mash distance matrix of g1, g2 and g3 (k = 9, sketch size 50).

**Expected Tool Call:**
```json
{"tool": "mash_triangle", "arguments": {"sequences": [">g1 first genome\nGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCAT\nGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGGTTAACG\n", ">g2 close\nGGATGACAGTCTACTCTGCTCACTCCAATCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATAACACTGGGTAGGATACGGCGGAGGGCACGTCAATACTGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTCGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCGTTCCAGATAGCGTTTCTGTTTCGGTGTGGGTGCTAATCGACTATGCTACAGCGGTTAACG\n", ">g3 distant\nGGATCATAGTCTACACTTCTCACTCCAACCCCGGCCCCCGCGTCCGAGGACACGATGCTTCAGAGCATGTATACCACTGCATAGGATACGGCGGATGGCACGTCAATAGGGTTCAATTCCCTACTGCATGCTCGTTAGGTTCATCCGGATTGAGAGGGTAGGCATGGGTGGAGGTGTTGGCTCGAGATCTGGAACTCCCATACAGAGCCCAATGTACCGAGTGTAGAGAGGCGCTTGTCCTTCCATTTATCGCGTCTCTTTCGGTATAGGTGCTAATTGACTATGCTACTGAGGTTAACG\n"], "names": ["g1.fa", "g2.fa", "g3.fa"], "k": 9, "sketchSize": 50}}
```

**Response:**
```json
{"text": "\t3\ng1.fa\ng2.fa\t0.0394246\ng3.fa\t0.212171\t0.359853\n"}
```

Equal to mash triangle -k 9 -s 50 g1.fa g2.fa g3.fa.

## See Also

- [mash_dist](mash_dist.md) - Distances of queries to references
- [mash_bounds](mash_bounds.md) - Error-bounds table
