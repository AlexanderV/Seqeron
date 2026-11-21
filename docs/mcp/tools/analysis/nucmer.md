# nucmer

MUMmer 3.23 nucmer: the whole pipeline (prenuc | mummer | mgaps | postnuc), the .delta file byte for byte.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `nucmer` |
| **Method ID** | `GenomicAnalyzer.Nucmer` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 nucmer, the whole pipeline: prenuc (the reference records joined with an x after each) | mummer (anchors: matchMode 'mum' = unique in reference and query, 'mumreference' = unique in the reference (default), 'maxmatch' = all; minMatch = -l; only A/C/G/T match; both strands unless forwardOnly / reverseOnly) | mgaps (clusters: maxGap = -g, minCluster = -c, diagonalDifference = -D, diagonalFactor = -d) | postnuc (the banded Smith-Waterman extension of the clusters: gaps between the matches of a cluster, then outward off both ends until the score has not improved for breakLength = -b bases; clusters that meet fuse; with simplify (default) clusters shadowed by an earlier alignment are dropped; extend / optimize / banded as --[no]extend / --[no]optimize (nooptimize forces the extension to the sequence ends) / --banded (hard diagonal band = diagonalDifference)). Returns the text of the .delta file - the alignments between every reference record and every query record as the distances between insertions and deletions (first line = the two file names, then NUCMER, then per record pair '>refId qryId refLen qryLen' and per alignment 'sR eR sQ eQ errors simErrors nonAlpha' followed by the delta values and 0; a reverse-strand alignment has sQ > eQ) - byte for byte what the installed nucmer 3.23 writes (434 runs: multi-record inputs, both strands, repeats, N runs, every option family, the abort paths). With delta = false (--nodelta) the .cluster text is returned instead. A break length outside 1 - 10000 is ignored with a warning (messages), as the program does. The inputs on which the programs abort (empty file, no header, a file whose only record is empty, an unparsable header) are refused with an error.

## Core Documentation Reference

- Source: [GenomicAnalyzer.Nucmer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.Nucmer.cs)
- Algorithm doc: [Nucmer_And_Delta_Programs.md](../../../algorithms/Alignment/Nucmer_And_Delta_Programs.md)
- Oracle: [generate_nucmer_oracle.py](../../../../scripts/alignment/generate_nucmer_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference FASTA text (multi-record). |
| `query` | string | Yes | Query FASTA text (multi-record). |
| `matchMode` | string | No | Anchor matches: 'mumreference' (default, --mumreference / --mumcand), 'mum' (--mum) or 'maxmatch' (--maxmatch). |
| `minMatch` | integer | No | Minimum length of a single match (-l, default 20). (default: 20) |
| `maxGap` | integer | No | Maximum gap between two adjacent matches of a cluster (-g, default 90). (default: 90) |
| `minCluster` | integer | No | Minimum length of a cluster of matches (-c, default 65). (default: 65) |
| `diagonalDifference` | integer | No | Maximum diagonal difference of two adjacent anchors (-D, default 5). (default: 5) |
| `diagonalFactor` | number | No | Diagonal difference as a fraction of the gap length (-d, default 0.12). (default: 0.12) |
| `breakLength` | integer | No | Distance an extension tries to extend poor scoring regions before giving up (-b, default 200; outside 1 - 10000 ignored with a warning). (default: 200) |
| `banded` | boolean | No | Enforce absolute banding of the dynamic programming matrix with diagonalDifference (--banded, default false). (default: false) |
| `forwardOnly` | boolean | No | Use only the forward strand of the query (-f, default false). (default: false) |
| `reverseOnly` | boolean | No | Use only the reverse complement of the query (-r, default false); both flags = both strands. (default: false) |
| `extend` | boolean | No | Extend the clusters (--extend, default true). (default: true) |
| `delta` | boolean | No | Create the delta (--delta, default true); false returns the cluster file (--nodelta). (default: true) |
| `optimize` | boolean | No | Alignment score optimisation (--optimize, default true); false = --nooptimize: an extension that reaches the end of a sequence is not backtracked. (default: true) |
| `simplify` | boolean | No | Remove shadowed alignments (--simplify, default true); false = --nosimplify (align a sequence to itself to find repeats). (default: true) |
| `referenceName` | string | No | Reference file name printed in the first line of the output (nucmer prints the absolute path). (default: reference.fa) |
| `queryName` | string | No | Query file name printed in the first line of the output. (default: query.fa) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `delta` | string | Text of the .delta file (null with delta = false). |
| `cluster` | string | Text of the .cluster file (only with delta = false). |
| `messages` | string | What the programs write to standard error (warnings). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | matchMode must be 'mumreference', 'mum' or 'maxmatch' |
| 1001 | mummer: empty sequence in multiple fasta file (and the other aborts of the programs) |

## Examples

### Example 1: Three substitutions between two sequences

**User Prompt:**
> Align these two sequences with nucmer, minimum match length 15.

**Expected Tool Call:**
```json
{"tool": "nucmer", "arguments": {"reference": ">ref\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "query": ">qry\nAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "minMatch": 15, "referenceName": "/tmp/b09nuc/ref.fa", "queryName": "/tmp/b09nuc/qry.fa"}}
```

**Response:**
```json
{"delta": "/tmp/b09nuc/ref.fa /tmp/b09nuc/qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n", "cluster": null, "messages": ""}
```

nucmer -l 15 ref.fa qry.fa writes exactly this out.delta: one alignment 6-140 / 1-135 with 3 errors, 3 similarity errors and no indel (the delta list is just the closing 0).

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
