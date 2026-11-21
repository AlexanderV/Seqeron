# delta2blocks

MUMmer delta2blocks (Debian package): the alignments of a nucmer delta as lines with a CIGAR of M / X / I / D.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `delta2blocks` |
| **Method ID** | `MummerDelta2Blocks.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

delta2blocks of the Debian MUMmer 3.23 package: the alignments of a nucmer delta file between one reference record and one query record as one line each, 'refId qryId sR eR sQ eQ CIGAR', the CIGAR written with M (matches), X (mismatches), I (an extra base in the reference, a positive delta) and D (an extra base in the query, a negative delta); the coordinates are those of the delta file (a reverse-strand alignment has sR > eR / sQ > eQ and is compared on the reverse complement). sortByReference (-r) / sortByQuery (-q) order the alignments. The sequences are read from the FASTA texts named in the delta header. A promer delta file is refused (the program reads the sequences through a null pointer and crashes). Byte for byte equal to the installed delta2blocks on 226 + runs.

## Core Documentation Reference

- Source: [MummerGaps.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerGaps.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_mummer_gaps_oracle.py](../../../../scripts/alignment/generate_mummer_gaps_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `delta` | string | Yes | Delta file text (nucmer). |
| `referenceFasta` | string | Yes | Reference FASTA text. |
| `queryFasta` | string | Yes | Query FASTA text. |
| `referenceId` | string | Yes | Reference record id. |
| `queryId` | string | Yes | Query record id. |
| `sortByQuery` | boolean | No | -q: sort the alignments by the query start (default false). |
| `sortByReference` | boolean | No | -r: sort the alignments by the reference start (default false). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What delta2blocks writes to standard output. |
| `messages` | string | What it writes to standard error. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Delta file text cannot be null or empty |
| 1001 | ERROR: Could not find any alignments for refId and qryId |
| 1001 | delta2blocks: the program reads the sequences of a promer alignment through a null pointer and crashes (SIGSEGV). |

## Examples

### Example 1: Three substitutions

**User Prompt:**
> Give the CIGAR of the alignment.

**Expected Tool Call:**
```json
{"tool": "delta2blocks", "arguments": {"delta": "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n", "referenceFasta": ">ref\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "queryFasta": ">qry\nAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "referenceId": "ref", "queryId": "qry"}}
```

**Response:**
```json
{"text": "ref qry 6 140 1 135 65M2X1M1X66M\n", "messages": ""}
```

delta2blocks d.delta ref qry: 65 matches, 2 mismatches, 1 match, 1 mismatch, 66 matches.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
