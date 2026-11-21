# delta2maf

MUMmer delta2maf (Debian package): the alignments of a nucmer delta file as MAF.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `delta2maf` |
| **Method ID** | `MummerDelta2Maf.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

delta2maf of the Debian MUMmer 3.23 package: the alignments of a nucmer delta file as a MAF file: '##maf version=1 scoring=single_cov2', per alignment 'a score=n' (n = |end - start| in the reference) and two 's' lines (source id, 0-based start, aligned length, strand, source length, the gapped sequence text in lower case; a reverse-strand alignment is written on the reverse complement with strand '-'), and '##eof maf' (no final newline). referenceId / queryId select the alignments of one reference record, one query record or one pair (a query id alone needs an empty referenceId); all alignments when both are empty. sortByReference (-r) / sortByQuery (-q) order the alignments. The reference and query files named in the delta header are read into one map by id (the first record of an id wins). An id without alignments is an error (the program prints the maf header and exits with 1); a promer delta file is refused (the program reads its amino-acid deltas as nucleotide deltas and runs past the sequences). Byte for byte equal to the installed delta2maf on 315 + runs.

## Core Documentation Reference

- Source: [MummerDelta2Maf.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerDelta2Maf.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_delta2maf_oracle.py](../../../../scripts/alignment/generate_delta2maf_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `delta` | string | Yes | Delta file text (nucmer). |
| `referenceFasta` | string | Yes | Reference FASTA text. |
| `queryFasta` | string | Yes | Query FASTA text. |
| `referenceId` | string | No | Only the alignments of this reference record (default: any). |
| `queryId` | string | No | Only the alignments of this query record (default: any; both ids must match when both are given). |
| `sortByQuery` | boolean | No | -q: sort the alignments by the query start (default false). |
| `sortByReference` | boolean | No | -r: sort the alignments by the reference start (default false). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What delta2maf writes to standard output. |
| `messages` | string | What it writes to standard error. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Delta file text cannot be null or empty |
| 1001 | ERROR: Could not find any alignments for refId and qryId |
| 1001 | delta2maf of a promer delta file: ... |

## Examples

### Example 1: One alignment as MAF

**User Prompt:**
> Write the alignment as MAF.

**Expected Tool Call:**
```json
{"tool": "delta2maf", "arguments": {"delta": "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n", "referenceFasta": ">ref\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "queryFasta": ">qry\nAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n"}}
```

**Response:**
```json
{"text": "##maf version=1 scoring=single_cov2\na score=134\ns ref 5 135 + 140 agacaattacataacatacacgtcagcacgaaacttgttggcccagtgtgaatcgcttaagggttaagtaagtgtgatgcatacgcctttacttgctgtgtccaccccatcggactggcatttttattacactca\ns qry 0 135 + 135 agacaattacataacatacacgtcagcacgaaacttgttggcccagtgtgaatcgcttaagggttggggaagtgtgatgcatacgcctttacttgctgtgtccaccccatcggactggcatttttattacactca\n\n##eof maf", "messages": ""}
```

delta2maf d.delta: one block, the reference starts at 0-based 5 (the alignment covers 6 - 140), three differences are visible in the aligned text.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
