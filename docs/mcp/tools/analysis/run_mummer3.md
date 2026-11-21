# run_mummer3

MUMmer run-mummer3: mummer | mgaps | combineMUMs, the files .out .gaps .align .errorsgaps, byte for byte.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `run_mummer3` |
| **Method ID** | `MummerRunPrograms.RunMummer3` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 run-mummer3, the script that runs the basic mummer 3 suite: mummer -mumreference -b -l 20 (maximal matches unique in the reference, both strands) | mgaps -l 100 -f .12 -s 600 (clusters) | combineMUMs -x -e .10 -W prefix.errorsgaps (extension of the clusters; no cover files). reference is a FASTA text (only its first record is aligned), query the multi-FASTA text; referenceName / queryName are the file names as the script receives them (printed in the messages of mummer). Returns the files of the script - out (the matches; the reverse matches are relative to the reverse complement of the query), gaps (the clusters), align (the alignments of the gaps: see combine_mums) and errorsGaps - the script's standard output (the three progress lines) and standard error (messages: the lines of mummer and the summary of combineMUMs). An empty input file or a query without a sequence stops the script (sh -e): refused. Byte for byte equal to the installed script on 216 + runs (files and messages) and 1 600 more random runs.

## Core Documentation Reference

- Source: [MummerRunPrograms.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerRunPrograms.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_combinemums_oracle.py](../../../../scripts/alignment/generate_combinemums_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference FASTA text (only the first record is aligned). |
| `query` | string | Yes | Multi-FASTA text of the query records. |
| `referenceName` | string | No | The reference file name as the script receives it. (default: reference.fa) |
| `queryName` | string | No | The query file name as the script receives it. (default: query.fa) |
| `prefix` | string | No | The prefix of the output files (names the files in the result). (default: out) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `out` | string | The matches (prefix.out). |
| `gaps` | string | The clusters (prefix.gaps). |
| `align` | string | The alignments of the gaps (prefix.align). |
| `errorsGaps` | string | The gaps with the number of errors (prefix.errorsgaps). |
| `stdout` | string | The standard output of the script. |
| `messages` | string | The standard error of the script. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | mummer: empty sequence in multiple fasta file (the script stops) |

## Examples

### Example 1: One query with a substitution and two indels

**User Prompt:**
> Run the mummer 3 suite on this pair.

**Expected Tool Call:**
```json
{"tool": "run_mummer3", "arguments": {"reference": ">ref demo\nggcccagtccagatcctcggaagtcccaactaacgaataagtagatccttctaaatagta\ngtatacgaatgtaatcaacatattaccgcacccaactttgattgggtcataaacaaacat\ncatgtttctggaggcggtttgtgagaccgacaaacgcgatagagtaatcttataaatcag\ngcaaatactcggtcctttctacgggtcgcggccgcgtggtcataggcagttttatgtttc\ntttgaggttaatgggaggca\n", "query": ">qry demo\nagatcctcggaagtcccaactaacgaataagtagatccttctaaatagtagtatacgaat\nataatcaacatattaccgcacccaactttgattgggtcataaacaaacatcatgtttctg\ngaggcggttttgagaccgacaaacgcgatagagtaatcttataaatcaggcaaatactcg\nggtcctttctacgggtcgcggccgcgtggtcataggcagttttatgtttctttgaggtta\n", "referenceName": "ref.fa", "queryName": "qry.fa", "prefix": "pf"}}
```

**Response:**
```json
{"out": "> qry\n      11         1        60\n      72        62        69\n     142       131        51\n     191       181        60\n> qry Reverse\n", "gaps": "> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     142      131     51    none      1      0\n     193      183     58      -2      0      1\n> qry Reverse\n", "align": "> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     Errors = 1\nA: gtatacgaatgtaatcaacat\nB: gtatacgaatataatcaacat\n   ==========^==========\n     142      131     51    none      1      0\n     Errors = 1\nA: gaggcggtttgtgagaccgac\nB: gaggcggttt-tgagaccgac\n   ==========^==========\n     193      183     58      -2      0      1\n     Errors = 1\nA: aaatactcgg-tcctttctac\nB: aaatactcgggtcctttctac\n   ==========^==========\nRegion:         11 .. 250                1 .. 240              3 / 240         1.25%\n> qry Reverse\n", "errorsGaps": "> qry\n      11        1     60    none      -      -       -\n      72       62     69    none      1      1       1\n     142      131     51    none      1      0       1\n     193      183     58      -2      0      1       1\n> qry Reverse\n", "stdout": "Find MUMs\nDetermine gaps\nAlign gaps\n", "messages": "# reading input file \"ref.fa\" of length 260\n# construct suffix tree for sequence of length 260\n# (maximum reference length is 536870908)\n# (maximum query length is 4294967295)\n# reading input file \"qry.fa\" of length 240\n# matching query-file \"qry.fa\"\n# against subject-file \"ref.fa\"\n           Ref len = 260\n            acgt's = 260\n        Non acgt's = 0\n Number of matches = 1\nSum of match bases = 240\n   Avg match bases = 240\n"}
```

run-mummer3 ref.fa qry.fa pf: four MUMs form one cluster; the three gaps between them are aligned and the region (240 bases of the query) has 3 errors.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
