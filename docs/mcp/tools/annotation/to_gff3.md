# to_gff3

Serialize gene annotations to GFF3 lines (with header).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `to_gff3` |
| **Method ID** | `GenomeAnnotator.ToGff3` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Serializes gene annotations into GFF3 text (Sequence Ontology GFF3 spec v1.26). The first line is the
`##gff-version 3` directive, followed by one tab-delimited 9-column line per annotation: `seqId`,
source (`.` when undefined), type, **1-based** start (`Start + 1`; an empty interval `[s, s)` becomes the
spec's zero-length feature `s..s`), end, score (`.` when undefined), strand, phase (accumulated 5'→3'
over the `CDS` segments sharing a `geneId`, `.` for other types), and column 9: `ID=geneId`,
`product=…`, then `attributes` (except `translation`) and `attributeValues` (written comma-separated);
`ID`/`product` are never repeated from the dictionaries and empty values are omitted. Percent-encoding
follows the spec: tab, newline, CR, `%` and control characters everywhere; `; = & ,` in column-9 tags
and values; every seqid character outside `[a-zA-Z0-9.:^*$@!+_?-|]`. A strand other than `+ - . ?`, a
negative start, an end before the start or the interval `[0, 0)` is rejected. This is the inverse of
[`parse_gff3`](parse_gff3.md).

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L752](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L752)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `annotations` | array | Yes | — | Gene annotations to serialize (≥1) |
| `seqId` | string | No | seq1 | Sequence identifier for column 1 |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `lines` | array | GFF3 lines; the first is `##gff-version 3` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Annotations cannot be null or empty |

## Examples

### Example 1: Single CDS to GFF3

A `CDS` annotation at 0-based start 99 → line 2 uses 1-based start 100, computed phase 0.

**Response:**
```json
{ "lines": [ "##gff-version 3", "chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=gene1;product=hypothetical protein;frame=1" ] }
```

## Performance

- **Time Complexity:** O(annotations)
- **Space Complexity:** O(annotations)

## See Also

- [parse_gff3](parse_gff3.md) - Inverse: parse GFF3 into features
- [predict_genes](predict_genes.md) - Produce annotations to serialize
