# combine_mums

MUMmer combineMUMs: extends the clusters of mgaps and aligns the gaps between their matches (edit distance), byte for byte.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `combine_mums` |
| **Method ID** | `MummerCombineMums.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 combineMUMs: takes the clusters of exact matches of mgaps (gaps = its output: '> id' / '> id Reverse' lines and 'start1 start2 length ...' lines; a '#' line ends a cluster) and extends every cluster off both ends and fills the gaps between its matches by edit-distance alignments (Myers / Landau-Vishkin diagonals; branch point score 0.272 per match and -0.728 per error in single precision, or 1 - (1 - errorRate) with errorRate = -e; give up 200 bases past the best score, at most 500 errors, extensions in steps of 10 000 bases). reference = a FASTA text (only its first record is read), query = the multi-FASTA texts of the query records; the n-th forward block of the gaps file belongs to the n-th query record and the reverse block after it to its reverse complement. Output (text) = the gaps file with, under every gap line, '     Errors = n' and the alignment (rows A: reference, B: query and a mark row with the 10 matching bases of each flanking match, 60 columns per block) or '   No extension', and a 'Region:  refStart .. refEnd  qryStart .. qryEnd  errors / length  percent' line after every combined region (a region on the reverse strand shows its query coordinates relative to the original query). onlyDifferences (-D) writes only the positions and characters of the differences instead; showDifferences (-S) every difference of the extensions and gaps; verbose (-v 1) every extension with its alignment; nucleotidesOnly (-n) lets only a, c, g, t match (others become '2' in the reference and '1' in the query); consecutiveNonAcgt (-N) breaks matches at that many consecutive non-ACGT characters; tagFromFasta (-t) tags the matches with the first word of the FASTA headers; queryTag (-q, default 'Query') / referenceTag (-r, default 'Ref') are suffixes of the cover file names. errorsFile = the file -W names (default witherrors.gaps): the gap lines with the number of errors appended. Without noCoverFiles (-x) the files '<reference base name><referenceTag>.cover' and '<query base name><queryTag>.cover' (regions and gaps with their N counts) are returned too; messages = standard error (a summary of the matches and the coverage of both sequences). A failed assertion of the program (no reference record, fewer query records than forward blocks, a line of the gaps file that is not three numbers) and a coordinate beyond a sequence (the program reads outside its buffers) are errors; the coverage report of a reference whose length is 7 modulo 16 reads the byte after its realloc'd buffer (valgrind: invalid read; the value printed depends on the allocator) and is refused unless noCoverFiles is set. Byte for byte equal to the installed combineMUMs on 240 + runs (every option, both strands, N / IUPAC characters, extensions of 24 kb) and 3 000 more random runs.

## Core Documentation Reference

- Source: [MummerCombineMums.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerCombineMums.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_combinemums_oracle.py](../../../../scripts/alignment/generate_combinemums_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference FASTA text (only the first record is read). |
| `query` | string | Yes | Multi-FASTA text of the query records. |
| `gaps` | string | Yes | The gaps file (output of mgaps). |
| `onlyDifferences` | boolean | No | -D: write only the positions and characters of the differences (default false). |
| `errorRate` | number | No | -e: the error rate cutoff of the extensions, e.g. 0.02 (default: the built-in branch point score). |
| `nucleotidesOnly` | boolean | No | -n: only a, c, g, t match (default false). |
| `consecutiveNonAcgt` | integer | No | -N: break matches at this many consecutive non-ACGT characters (default 0 = never). |
| `queryTag` | string | No | -q: suffix of the query tag (default Query). |
| `referenceTag` | string | No | -r: suffix of the reference tag (default Ref). |
| `showDifferences` | boolean | No | -S: show every difference of the extensions and gaps (default false). |
| `tagFromFasta` | boolean | No | -t: tag the matches with the first word of the FASTA headers (default false). |
| `verbose` | integer | No | -v: verbose level, above 0 writes every extension with its alignment (default 0). |
| `noCoverFiles` | boolean | No | -x: do not write the cover files (default false). |
| `errorsFileName` | string | No | -W: name of the errors file (default witherrors.gaps). |
| `referenceName` | string | No | The reference file name as given on the command line (its base name is the reference tag; default reference.fa). (default: reference.fa) |
| `queryName` | string | No | The query file name as given on the command line (its base name is the query tag; default query.fa). (default: query.fa) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What combineMUMs writes to standard output. |
| `errorsFile` | string | The file -W names (witherrors.gaps): the gap lines with the number of errors appended. |
| `messages` | string | What the program writes to standard error: the summary of the matches and, without -x, the coverage of both sequences. |
| `files` | object | The files the program leaves behind by name: the errors file and the .cover files. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Gaps file text cannot be null |
| 1001 | combineMUMs: Assertion `Read_String (...)' failed (no reference record, fewer query records than forward blocks) |
| 1001 | combineMUMs: Assertion `sscanf (...) == 3' failed (a line of the gaps file that is not three numbers) |

## Examples

### Example 1: Alignment of the gaps of one query

**User Prompt:**
> Extend and combine the clusters of this gaps file (no cover files).

**Expected Tool Call:**
```json
{"tool": "combine_mums", "arguments": {"reference": ">ref demo\nggcccagtccagatcctcggaagtcccaactaacgaataagtagatccttctaaatagta\ngtatacgaatgtaatcaacatattaccgcacccaactttgattgggtcataaacaaacat\ncatgtttctggaggcggtttgtgagaccgacaaacgcgatagagtaatcttataaatcag\ngcaaatactcggtcctttctacgggtcgcggccgcgtggtcataggcagttttatgtttc\ntttgaggttaatgggaggca\n", "query": ">qry demo\nagatcctcggaagtcccaactaacgaataagtagatccttctaaatagtagtatacgaat\nataatcaacatattaccgcacccaactttgattgggtcataaacaaacatcatgtttctg\ngaggcggttttgagaccgacaaacgcgatagagtaatcttataaatcaggcaaatactcg\nggtcctttctacgggtcgcggccgcgtggtcataggcagttttatgtttctttgaggtta\n", "gaps": "> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     142      131     51    none      1      0\n     193      183     58      -2      0      1\n> qry Reverse\n", "noCoverFiles": true, "referenceName": "ref.fa", "queryName": "qry.fa"}}
```

**Response:**
```json
{"text": "> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     Errors = 1\nA: gtatacgaatgtaatcaacat\nB: gtatacgaatataatcaacat\n   ==========^==========\n     142      131     51    none      1      0\n     Errors = 1\nA: gaggcggtttgtgagaccgac\nB: gaggcggttt-tgagaccgac\n   ==========^==========\n     193      183     58      -2      0      1\n     Errors = 1\nA: aaatactcgg-tcctttctac\nB: aaatactcgggtcctttctac\n   ==========^==========\nRegion:         11 .. 250                1 .. 240              3 / 240         1.25%\n> qry Reverse\n", "errorsFile": "> qry\n      11        1     60    none      -      -       -\n      72       62     69    none      1      1       1\n     142      131     51    none      1      0       1\n     193      183     58      -2      0      1       1\n> qry Reverse\n", "messages": "           Ref len = 260\n            acgt's = 260\n        Non acgt's = 0\n Number of matches = 1\nSum of match bases = 240\n   Avg match bases = 240\n", "files": {"witherrors.gaps": "> qry\n      11        1     60    none      -      -       -\n      72       62     69    none      1      1       1\n     142      131     51    none      1      0       1\n     193      183     58      -2      0      1       1\n> qry Reverse\n"}}
```

combineMUMs -x ref.fa qry.fa pf.gaps: the three gaps of the cluster are aligned (one substitution, one extra base in the query, one extra base in the reference), the region is extended to the ends of the query (240 bases, 3 errors).

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
