# mash_taxscreen

Mash 2.3 'mash taxscreen': the Kraken-style taxonomic report of query references contained in a pool, byte for byte.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mash_taxscreen` |
| **Method ID** | `MashTaxScreen.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Mash 2.3 'mash taxscreen', byte for byte (CommandTaxScreen.cpp + taxdb.hpp): the Kraken-style taxonomic report of how well query references are contained in a pool of sequences. The pass over the pool is that of mash_screen (reads or contigs; 6-frame translation for amino-acid queries); the taxonomy ID of every query reference is the one of the mapping file (lines 'taxid SEP reference-name') or the number after the word 'taxid' in its comment, the hashes of references that share a hash are assigned to their lowest common ancestor in the taxonomy (names.dmp / nodes.dmp, the NCBI taxonomy dump), and the report lists from the root: percent of the hashes, number of contained hashes in the clade, in the taxon, number of hashes of the references in the clade and the taxon, rank, taxonomy ID and the name indented by two blanks per level, children by decreasing clade count. The report is the program's, including its order of equal counts (std::sort of libstdc++) and the effect of its unordered_map growing while the clade counts are summed. Reads the dump from the texts given (no download: a taxonomy of your own in the dump format works). A node whose parent is missing from nodes.dmp is refused (the program uses an uninitialised pointer there: valgrind-verified). VERIFICATION: 300 / 300 random runs (synthetic taxonomies of 4 - 300 taxa in the NCBI dump format: random, deep and star-shaped trees with up to 45 equal children, nucleotide and amino-acid queries, mapping files, taxonomy IDs in comments, IDs missing from the dump) equal the installed mash 2.3.

## Core Documentation Reference

- Source: [MashTaxScreen.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MashTaxScreen.cs)
- Algorithm doc: [Mash_Commands.md](../../../algorithms/Alignment/Mash_Commands.md)
- Oracle: [generate_mash_taxscreen_oracle.py](../../../../scripts/alignment/generate_mash_taxscreen_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `queries` | array<string> | Yes | The query files as FASTA texts, or one sketch file 'msh:' + base64 (the output of mash_sketch); the reference names are the file names (or the sequence ids with individual). |
| `pools` | array<string> | Yes | The pool files as FASTA texts (reads or contigs, nucleotides). |
| `namesDump` | string | Yes | The content of names.dmp (tax_id \| name \| unique name \| name class \|, one taxon per line). |
| `nodesDump` | string | Yes | The content of nodes.dmp (tax_id \| parent tax_id \| rank \| ...). |
| `mappingFile` | string | No | -m: the mapping file ('taxonomy ID, one separator character, reference name' per line); null = the IDs come from the comments of the references ('taxid N'). |
| `queryNames` | array<string> | No | Names of the query files (default query1, ...). |
| `k` | integer | No | k-mer size of the query sketch (-k of mash sketch, 1 - 32; default 21, or 9 with protein). |
| `sketchSize` | integer | No | Sketch size of the query sketch (-s, default 1000). (default: 1000) |
| `seed` | integer | No | Hash seed of the query sketch (-S, 0 - 4294967295, default 42). (default: 42) |
| `noncanonical` | boolean | No | The query sketch was made with -n (strand preserved). (default: false) |
| `individual` | boolean | No | The query sketch was made with -i (one reference per sequence). (default: false) |
| `protein` | boolean | No | -a: the amino-acid alphabet of the query sketch (implies -n and, without k, -k 9; the pool is translated in six frames). (default: false) |
| `alphabet` | string | No | -z: the alphabet of the query sketch (implies -n). |
| `preserveCase` | boolean | No | -Z: the case of the letters is preserved in the query sketch. (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The report: the header line and one line per taxon from the root down. |

## Errors

| Code | Message |
|------|---------|
| 1001 | A node whose parent is missing from nodes.dmp (the program uses an uninitialised pointer: refused) |
| 1001 | At least one input is required |
| 1001 | Input has no FASTA record |

## Examples

### Example 1: Two references of one genus screened against a read of the first

**User Prompt:**
> Screen this read against the two reference genomes and report the taxonomy.

**Expected Tool Call:**
```json
{"tool": "mash_taxscreen", "arguments": {"queries": [">refA taxid 10\nGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTT\n", ">refB taxid 11\nGGATCACAGTCTACACTGCTCACTCCAACCCAGAGTATGTATACCACTGGGTAGGATACG\n"], "pools": [">read1\nGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTT\n"], "namesDump": "1\t|\troot\t|\t\t|\tscientific name\t|\n5\t|\tGenusX\t|\t\t|\tscientific name\t|\n10\t|\tGenusX alpha\t|\t\t|\tscientific name\t|\n11\t|\tGenusX beta\t|\t\t|\tscientific name\t|\n", "nodesDump": "1\t|\t1\t|\tno rank\t|\n5\t|\t1\t|\tgenus\t|\n10\t|\t5\t|\tspecies\t|\n11\t|\t5\t|\tspecies\t|\n", "queryNames": ["a.fa", "b.fa"], "k": 9, "sketchSize": 100}}
```

**Response:**
```json
{"text": "%\thashes\ttaxHashes\thashesDB\ttaxHashesDB\ttaxID\trank\tname\n100.0000\t52\t0\t81\t0\tno rank\t1\troot\n100.0000\t52\t23\t81\t23\tgenus\t5\t  GenusX\n55.7692\t29\t29\t29\t29\tspecies\t10\t    GenusX alpha\n"}
```

Both references share the first 30 bases, so their common hashes (23) are assigned to the genus; the read holds all of refA (29 more hashes).

## See Also

- [mash_screen](mash_screen.md) - the containment screen that taxscreen builds on
- [mash_sketch](mash_sketch.md) - the query sketch file
