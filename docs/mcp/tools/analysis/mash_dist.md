# mash_dist

Mash distances of queries to references (mash dist) with the -t, -v, -d, -i, -n options.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mash_dist` |
| **Method ID** | `MashCommands.Dist` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Mash 2.3 'mash dist': MinHash (bottom-s) sketches of the reference and query files (canonical k-mers, MurmurHash3, the sketch of Ondov et al. 2016: KmerAnalyzer.CreateMinHashSketch) and, for every query in order and every reference, the line 'reference TAB query TAB distance TAB p-value TAB shared/total' of the pairs whose p-value is at most maxPValue (-v) and distance at most maxDistance (-d); table = true (-t) writes the '#query' header and a row of distances per query (a pair that fails a filter leaves its field blank); individual = true (-i) sketches every sequence on its own, named by its id (sequences shorter than k are skipped). A file without a record of at least k symbols is refused as the program refuses it (exit status 1: 'All fasta records in input files were shorter than the k-mer size'). Numbers are printed as C++ cout does (6 significant digits). The alphabet options -a (protein: amino acids, implies -n -k 9), -z (own alphabet, implies -n) and -Z (case kept) as in the program (k-mer space alphabetSize^k in the p-value, 64-bit hashes when alphabetSize^k > 2^32). VERIFICATION: 220 / 220 random sets (1 - 4 references x 1 - 3 queries of 1 - 3 records, k 4 - 32, sketch size 5 - 2000, seeds, -n, -i, -t, -v, -d) and 300 / 300 sets with -a / -z / -Z equal the installed mash 2.3 byte for byte.

## Core Documentation Reference

- Algorithm: [Mash_Commands.md](../../algorithms/Alignment/Mash_Commands.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `references` | array<string> | Yes | The reference files as FASTA texts (records start with '>': the first word of the header is the sequence id, the rest the comment). |
| `queries` | array<string> | Yes | The query files (FASTA texts). |
| `referenceNames` | array<string> | No | Names of the reference files as the program prints them (default reference1, ...). |
| `queryNames` | array<string> | No | Names of the query files (default query1, ...). |
| `k` | integer | No | k-mer size (-k, 1 - 32; default 21, or 9 with protein). |
| `sketchSize` | integer | No | Sketch size (-s, default 1000). (default: 1000) |
| `seed` | integer | No | Hash seed (-S, 0 - 4294967295, default 42). (default: 42) |
| `noncanonical` | boolean | No | -n: preserve the strand (no canonical k-mers). (default: false) |
| `individual` | boolean | No | -i: sketch individual sequences instead of whole files. (default: false) |
| `table` | boolean | No | -t: table output. (default: false) |
| `maxPValue` | number | No | -v: maximum p-value to report (0 - 1, default 1). (default: 1.0) |
| `maxDistance` | number | No | -d: maximum distance to report (0 - 1, default 1). (default: 1.0) |
| `reads` | boolean | No | -r: the inputs are read sets (noise filter -m, target coverage -c, genome size -g); incompatible with individual. (default: false) |
| `minCopies` | integer | No | -m: minimum copies of each k-mer required to pass the noise filter of a read set (default 1; only with reads). |
| `targetCoverage` | number | No | -c: target coverage of a read set: sketching concludes when it is reached (default 0 = read everything; only with reads). |
| `genomeSize` | integer | No | -g: genome size used as the sketch length of a read set (default 0 = estimated; only with reads). |
| `showComments` | boolean | No | -C: show the comment fields with the names (name:comment for dist, the comment instead of the id for triangle; the comment of a file sketch is the first record's name and comment, '[N seqs] ... [...]' for several). (default: false) |
| `bloomFilter` | string | No | -b: a Bloom filter of this size (raw bytes or with K/M/G/T, e.g. 500M) filters out the k-mers seen once, instead of the exact filter -m (cannot be combined with minCopies > 1; implies reads; copies cannot be counted beyond 2; some unique k-mers pass erroneously; the filter has one hash function over bytes x 8 bits as the program builds it). |
| `protein` | boolean | No | -a: the amino-acid alphabet (A-Z except BJOUXZ; implies -n and, without k, -k 9): k-mers of the 20 letters, 64-bit hashes when 20^k > 2^32, p-values over the k-mer space 20^k. (default: false) |
| `alphabet` | string | No | -z: the alphabet the hashes are based on (case ignored unless preserveCase; k-mers with other characters are ignored; implies -n; ASCII only: the program has no defined result for a byte of 0x80 or more). protein takes precedence. |
| `preserveCase` | boolean | No | -Z: preserve the case of the letters in the k-mers and the alphabet (a letter whose case is not in the alphabet is skipped). (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text mash dist writes. |

## Errors

| Code | Message |
|------|---------|
| 1001 | k must be in 1..32 (mash -k) |
| 1001 | seed must be in 0..4294967295 (mash -S) |
| 1001 | -v must be in 0..1 (mash) |

## Examples

### Example 1: One reference, two queries

**User Prompt:**
> Estimate the Mash distance of g2 and g3 to g1 (k = 9, sketch size 50).

**Expected Tool Call:**
```json
{"tool": "mash_dist", "arguments": {"references": [">g1 first genome\nGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCAT\nGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGGTTAACG\n"], "queries": [">g2 close\nGGATGACAGTCTACTCTGCTCACTCCAATCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATAACACTGGGTAGGATACGGCGGAGGGCACGTCAATACTGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTCGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCGTTCCAGATAGCGTTTCTGTTTCGGTGTGGGTGCTAATCGACTATGCTACAGCGGTTAACG\n", ">g3 distant\nGGATCATAGTCTACACTTCTCACTCCAACCCCGGCCCCCGCGTCCGAGGACACGATGCTTCAGAGCATGTATACCACTGCATAGGATACGGCGGATGGCACGTCAATAGGGTTCAATTCCCTACTGCATGCTCGTTAGGTTCATCCGGATTGAGAGGGTAGGCATGGGTGGAGGTGTTGGCTCGAGATCTGGAACTCCCATACAGAGCCCAATGTACCGAGTGTAGAGAGGCGCTTGTCCTTCCATTTATCGCGTCTCTTTCGGTATAGGTGCTAATTGACTATGCTACTGAGGTTAACG\n"], "referenceNames": ["g1.fa"], "queryNames": ["g2.fa", "g3.fa"], "k": 9, "sketchSize": 50}}
```

**Response:**
```json
{"text": "g1.fa\tg2.fa\t0.0394246\t2.98721e-74\t27/50\ng1.fa\tg3.fa\t0.212171\t2.41195e-08\t4/50\n"}
```

Equal to mash dist -k 9 -s 50 on the three files (the reference sketched with mash sketch -k 9 -s 50).

## See Also

- [mash_triangle](mash_triangle.md) - Lower-triangular distance matrix
- [mash_bounds](mash_bounds.md) - Error-bounds table
