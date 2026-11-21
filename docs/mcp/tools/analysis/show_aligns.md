# show_aligns

MUMmer 3.23 show-aligns: the alignments of one reference / query pair drawn as three rows.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `show_aligns` |
| **Method ID** | `MummerShowAligns.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 show-aligns: the alignments between one reference record and one query record of a nucmer delta file drawn as three rows (reference, query, a ^ under every difference; an indel is a . in the other row) in blocks of the screen width (-w, default 60; 11 or less is replaced with a warning), each row prefixed with the 1-based forward-strand position of its first base; '-- BEGIN alignment [ +1 sR - eR | +1 sQ - eQ ]' / END lines frame every alignment (a reverse-strand alignment shows -1). The alignments are in file order unless sortByReference (-r) / sortByQuery (-q). The program reads the sequences from the FASTA files named in the delta header: pass their text. Byte for byte equal to the installed show-aligns on 300 + runs; a promer (amino-acid) delta file is refused with NotSupportedException (registered A3-9d).

## Core Documentation Reference

- Source: [MummerShowAligns.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerShowAligns.cs)
- Algorithm doc: [Nucmer_And_Delta_Programs.md](../../../algorithms/Alignment/Nucmer_And_Delta_Programs.md)
- Oracle: [generate_mummer_show_oracle.py](../../../../scripts/alignment/generate_mummer_show_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `delta` | string | Yes | Delta file text (nucmer). |
| `referenceFasta` | string | Yes | Reference FASTA text (the file named in the delta header). |
| `queryFasta` | string | Yes | Query FASTA text (the file named in the delta header). |
| `referenceId` | string | Yes | Reference record id. |
| `queryId` | string | Yes | Query record id. |
| `sortByQuery` | boolean | No | -q: sort the alignments by the query start (default false). |
| `sortByReference` | boolean | No | -r: sort the alignments by the reference start (default false). |
| `screenWidth` | integer | No | -w: screen width (default 60). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What show-aligns writes to standard output. |
| `messages` | string | What it writes to standard error (warnings). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Delta file text cannot be null or empty |
| 1001 | ERROR: Could not find any alignments for <ref> and <qry> |

## Examples

### Example 1: The three substitutions drawn

**User Prompt:**
> Show the alignment between ref and qry in 50 columns.

**Expected Tool Call:**
```json
{"tool": "show_aligns", "arguments": {"delta": "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n", "referenceFasta": ">ref\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "queryFasta": ">qry\nAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "referenceId": "ref", "queryId": "qry", "screenWidth": 50}}
```

**Response:**
```json
{"text": "ref.fa qry.fa\n\n==================================================\n-- Alignments between ref and qry\n\n-- BEGIN alignment [ +1 6 - 140 | +1 1 - 135 ]\n\n\n6          agacaattacataacatacacgtcagcacgaaacttgtt\n1          agacaattacataacatacacgtcagcacgaaacttgtt\n                                                  \n\n45         ggcccagtgtgaatcgcttaagggttaagtaagtgtgat\n40         ggcccagtgtgaatcgcttaagggttggggaagtgtgat\n                                     ^^ ^         \n\n84         gcatacgcctttacttgctgtgtccaccccatcggactg\n79         gcatacgcctttacttgctgtgtccaccccatcggactg\n                                                  \n\n123        gcatttttattacactca\n118        gcatttttattacactca\n                             \n\n\n--   END alignment [ +1 6 - 140 | +1 1 - 135 ]\n\n==================================================\n", "messages": ""}
```

show-aligns -w 50 d1.delta ref qry: a ^ marks each of the three substituted bases.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
