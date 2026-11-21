# promer

MUMmer 3.23 promer: amino-acid alignments between two DNA files (six-frame translation, mummer, mgaps, postpro), the .delta file byte for byte.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `promer` |
| **Method ID** | `PromerProgram.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 promer (PROtein MUMmer 3.07): amino-acid alignments between two multi-FASTA DNA files, the whole pipeline: prepro (both files translated in all six reading frames; a stretch of at most maskLength amino acids between stop codons is masked) | mummer (anchors in amino acids: matchMode 'mumreference' (default), 'mum', 'maxmatch'; minMatch = -l, default 6) | mgaps (clusters: maxGap = -g (30), minCluster = -c (20), diagonalFactor = -d (0.11)) | postpro (the affine Smith-Waterman extension of the clusters in amino-acid space with the BLOSUM matrix of matrix (-x: 1 = BLOSUM 45, 2 = BLOSUM 62 default, 3 = BLOSUM 80; 0 and other values are refused with a warning), give-up distance breakLength = -b (60 amino acids); extend / optimize as --[no]extend / --[no]optimize). Returns the .delta file: first line = the two file names, then PROMER, per record pair '>refId qryId refLen qryLen' and per alignment 'sR eR sQ eQ errors similarityErrors stopCodons' in nucleotide coordinates of the original DNA (a reverse-strand alignment has sR > eR / sQ > eQ) followed by the delta values - amino-acid distances, one delta = three nucleotides - and 0. With delta = false (--nodelta) the .cluster text is returned instead; coords (-o) adds the PROmer 1.1 coords file (show-coords -r of the delta). The prefix (-p) names the files of the result and the temporary files in the messages. Messages are the progress lines of the script, the warnings of the programs ('Forcing unrecognized DNA char to N', invalid break length / matrix) and mummer's reading lines. Byte for byte equal to the installed promer on 320 recorded and 1 920 further random runs (delta, cluster, coords, the temporary files and the messages; every option). The inputs on which the pipeline aborts (a file without sequence, a non-positive minMatch) are refused with an error.

## Core Documentation Reference

- Source: [PromerProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/PromerProgram.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_promer_oracle.py](../../../../scripts/alignment/generate_promer_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference FASTA text (multi-record DNA). |
| `query` | string | Yes | Query FASTA text (multi-record DNA). |
| `matchMode` | string | No | Anchor matches: 'mumreference' (default, --mumreference / --mumcand), 'mum' (--mum) or 'maxmatch' (--maxmatch). |
| `minMatch` | integer | No | Minimum length of a single match in amino acids (-l, default 6). (default: 6) |
| `maxGap` | integer | No | Maximum gap between two adjacent matches of a cluster in amino acids (-g, default 30). (default: 30) |
| `minCluster` | integer | No | Minimum length of a cluster of matches in amino acids (-c, default 20). (default: 20) |
| `diagonalFactor` | number | No | Diagonal difference as a fraction of the gap length (-d, default 0.11). (default: 0.11) |
| `breakLength` | integer | No | Distance an extension tries to extend poor scoring regions before giving up, in amino acids (-b, default 60; outside 1 - 10000 ignored with a warning). (default: 60) |
| `matrix` | integer | No | Alignment matrix (-x): 1 = BLOSUM 45, 2 = BLOSUM 62 (default), 3 = BLOSUM 80. |
| `maskLength` | integer | No | Maximum bookend masking length in amino acids (-m, default 8; negative = invalid, ignored with a warning). (default: 8) |
| `extend` | boolean | No | Extend the clusters (--extend, default true). (default: true) |
| `delta` | boolean | No | Create the delta (--delta, default true); false returns the cluster file (--nodelta). (default: true) |
| `optimize` | boolean | No | Alignment score optimisation (--optimize, default true); false = --nooptimize. (default: true) |
| `coords` | boolean | No | Also return the PROmer 1.1 coords file (-o / --coords: show-coords -r of the delta; default false; not with delta = false). (default: false) |
| `prefix` | string | No | Output prefix (-p, default out): names the files in the messages. (default: out) |
| `referenceName` | string | No | Reference file name printed in the first line of the delta (promer prints the absolute path). (default: reference.fa) |
| `queryName` | string | No | Query file name printed in the first line of the delta. (default: query.fa) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `delta` | string | Text of the .delta file (null with delta = false). |
| `cluster` | string | Text of the .cluster file (only with delta = false). |
| `messages` | string | The progress lines of the script, the warnings of the programs and mummer's reading lines. |
| `coords` | string | Text of the coords file (only with coords = true). |
| `aminoReference` | string | The temporary .aaref file (the reference in six frames). |
| `aminoQuery` | string | The temporary .aaqry file (the query in six frames). |
| `mgaps` | string | The temporary .mgaps file (the clusters). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | matchMode must be 'mumreference', 'mum' or 'maxmatch' |
| 1001 | ERROR: mummer and/or mgaps returned non-zero |
| 1001 | ERROR: Could not parse input from 'stdin'. |

## Examples

### Example 1: A gene with substitutions and an insertion

**User Prompt:**
> Align the two genes at the protein level.

**Expected Tool Call:**
```json
{"tool": "promer", "arguments": {"reference": ">ref gene\nACGTTGCATATAAAAACACTGCACGGATTTGCGGAGAGAATTCGATTCCCGTCGAATGGATCGCACACCAGAAGGGTCCAGGAGATTTCTATTATCGATTCTTTGCAGCACATCACGGAGGATTGATGCATGTACACGGCCATCCTCTAGCGAATCAAGGATTCAAGGACATGTTATATTGGGCATATATGGATCTGAACTTATCCATGGGATCCTCTGGATC\n", "query": ">qry\nTTTACAAGAATACGGCCAGAATCTGCGGAGAGTGGTCAATCCCGGTCGAATGGATAGCCCACCAGAAAGGTCCGGGAGACTTTTACTATCGATTCTTCGCGGGACACCACGGCGGTCTCATGCACGTACACGGACACCCACTTGCAAACCAGAAAAAGGGATTCAAAGATATGCTTTACTGGGCTTACATGGACCTTAACCTGTCTATGGGCAGTTCGACGTA\n", "minMatch": 5, "minCluster": 10, "referenceName": "ref.fa", "queryName": "qry.fa"}}
```

**Response:**
```json
{"delta": "ref.fa qry.fa\nPROMER\n>ref qry 223 223\n9 218 3 218 4 4 0\n-51\n-1\n0\n217 11 217 5 22 15 8\n-21\n-1\n0\n"}
```

promer -l 5 -c 10 ref.fa qry.fa: two alignments, forward in frame +3 (an insertion of two residues in the query: delta -51, -1) and on the reverse strand; the delta values count amino acids.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
