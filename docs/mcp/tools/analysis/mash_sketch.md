# mash_sketch

Mash 2.3 'mash sketch -o x.msh': the sketch file of the inputs, byte for byte as the program writes it (120 random runs equal: genomes and read sets, k 4 - 32, sketch sizes 1 - 3000, seeds, -n, -i, -r -m -c -g with all the files of a read set in ONE sketch, -I / -C, -l). The file is a Cap'n Proto message of MinHash.capnp built as capnp::MallocMessageBuilder builds it (segments of 1024, 1024, 2048, ... words, far pointers with landing pads); it is returned as base64 and, prefixed with 'msh:', it is accepted by mash_dist / mash_triangle / mash_screen / mash_info / mash_paste as an input. One sketch per file (its name is the file name, its comment the first record's name and comment) or, with individual, one per sequence. -l (a list of file names) is the caller's list of inputs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mash_sketch` |
| **Method ID** | `MashCommands.Sketch` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Mash 2.3 'mash sketch -o x.msh': the sketch file of the inputs, byte for byte as the program writes it (120 random runs equal: genomes and read sets, k 4 - 32, sketch sizes 1 - 3000, seeds, -n, -i, -r -m -c -g with all the files of a read set in ONE sketch, -I / -C, -l). The file is a Cap'n Proto message of MinHash.capnp built as capnp::MallocMessageBuilder builds it (segments of 1024, 1024, 2048, ... words, far pointers with landing pads); it is returned as base64 and, prefixed with 'msh:', it is accepted by mash_dist / mash_triangle / mash_screen / mash_info / mash_paste as an input. One sketch per file (its name is the file name, its comment the first record's name and comment) or, with individual, one per sequence. -l (a list of file names) is the caller's list of inputs.

## Core Documentation Reference

- Algorithm: [Mash_Commands.md](../../../algorithms/Alignment/Mash_Commands.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | The input files as FASTA texts (records start with '>': the first word of the header is the sequence id, the rest the comment). |
| `names` | array<string> | No | Names of the inputs as the program stores them (the file names; default input1, ...). |
| `k` | integer | No | k-mer size (-k, 1 - 32, default 21). (default: 21) |
| `sketchSize` | integer | No | Sketch size (-s, default 1000). (default: 1000) |
| `seed` | integer | No | Hash seed (-S, 0 - 4294967295, default 42; 42 is stored in referenceListOld, another seed in referenceList). (default: 42) |
| `noncanonical` | boolean | No | -n: preserve the strand (no canonical k-mers). (default: false) |
| `individual` | boolean | No | -i: sketch individual sequences instead of whole files. (default: false) |
| `reads` | boolean | No | -r: the inputs are read sets (noise filter -m, target coverage -c, genome size -g); the hash counts are stored. (default: false) |
| `minCopies` | integer | No | -m: minimum copies of each k-mer required to pass the noise filter of a read set (default 1; only with reads). |
| `targetCoverage` | number | No | -c: target coverage of a read set (default 0 = read everything; only with reads). |
| `genomeSize` | integer | No | -g: genome size used as the sketch length of a read set (default 0 = estimated; only with reads). |
| `id` | string | No | -I: the id (name) of the first sketch. |
| `comment` | string | No | -C: the comment of the first sketch. |
| `bloomFilter` | string | No | -b: a Bloom filter of this size (raw bytes or with K/M/G/T, e.g. 500M) filters out the k-mers seen once, instead of the exact filter -m (cannot be combined with minCopies > 1; implies reads; copies cannot be counted beyond 2; some unique k-mers pass erroneously; the filter has one hash function over bytes x 8 bits as the program builds it). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `msh` | string | msh: + base64 of the .msh file. |
| `bytes` | integer | Size of the file. |
| `sketches` | integer | Number of sketches. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (the program's message for incompatible options or a sketch file that is not a .msh). |

## Examples

### Example 1: Sketch one genome

**User Prompt:**
> Sketch this genome with Mash (k 8, 5 hashes).

**Expected Tool Call:**
```json
{"tool": "mash_sketch", "arguments": {"sequences": [">s1\nACGTACGTAGCTAGCTAGGATCGATCGATCGTAGCTAGCTAGTCGATCGATGCTAGCTAGC\n"], "names": ["a.fa"], "k": 8, "sketchSize": 5}}
```

**Response:**
```json
{"msh": "msh:...", "bytes": 216, "sketches": 1}
```

## See Also

- [mash_dist](mash_dist.md) - distances
- [mash_screen](mash_screen.md) - containment
