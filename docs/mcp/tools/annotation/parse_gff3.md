# parse_gff3

Parse a GFF3 annotation document into structured features.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `parse_gff3` |
| **Method ID** | `GenomeAnnotator.ParseGff3` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses a GFF3 document (Sequence Ontology GFF3 spec v1.26; lines separated by newlines, a trailing CR
is dropped) into structured genomic features. Blank lines and lines starting with `#` (comments and `##`
directives) are skipped, and reading stops at `##FASTA` or at a line beginning with `>` (implied
`##FASTA`). A feature line is accepted only when it follows the spec's column rules — the same rules
GenomeTools `gt gff3validator` enforces — and is skipped otherwise: exactly 9 tab-separated columns,
start/end positive integers with start ≤ end, score `.` or a number, strand one of `+ - . ?`, phase `.`
or `0`/`1`/`2`, column 9 `.` or `tag=value` pairs with non-empty tag and value and no repeated tag.
Column 6 (`score`) and column 8 (`phase`) map to `null` for `.`. Columns 1–3 and every tag/value are
percent-decoded; values keep their spaces. `attributes` maps each tag to its decoded value text
(multi-values joined by `,`); `attributeValues` keeps the value list, split on the unescaped `,`
(`Parent=a,b` → `["a","b"]`, `Note=a%2Cb` → `["a,b"]`). The feature id is the `ID` attribute,
falling back to `feature_{n}`. Coordinates stay 1-based inclusive.

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L565](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L565)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gff3Text` | string | Yes | GFF3 document text (lines separated by newlines) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `features` | array | `{ featureId, type, start, end, strand, score, phase, attributes, seqid, source, attributeValues }` per feature |

## Errors

| Code | Message |
|------|---------|
| 1001 | GFF3 text cannot be null or empty |

## Examples

### Example 1: Fully specified CDS line

`chr1  ENSEMBL  CDS  1000  2000  95.5  -  2  ID=cds1;Name=TestCDS` →

**Response:**
```json
{ "features": [ { "featureId": "cds1", "type": "CDS", "start": 1000, "end": 2000, "score": 95.5, "strand": "-", "phase": 2, "attributes": { "ID": "cds1", "Name": "TestCDS" }, "seqid": "chr1", "source": "ENSEMBL", "attributeValues": { "ID": ["cds1"], "Name": ["TestCDS"] } } ] }
```

### Example 2: Undefined score/phase

`chr1  src  gene  1  100  .  +  .  ID=gene1` yields `score: null` and `phase: null`.

### Example 3: Multi-valued attribute (spec canonical gene, line 8)

`ctg123  .  exon  1050  1500  .  +  .  ID=exon00002;Parent=mRNA00001,mRNA00002` yields
`attributeValues.Parent = ["mRNA00001", "mRNA00002"]` (and `attributes.Parent = "mRNA00001,mRNA00002"`).

## Performance

- **Time Complexity:** O(lines)
- **Space Complexity:** O(features)

## See Also

- [to_gff3](to_gff3.md) - Serialize annotations back to GFF3
- [predict_genes](predict_genes.md) - ORF-based gene prediction
