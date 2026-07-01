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

Parses a GFF3 document (lines separated by newlines) into structured genomic features. Blank lines
and lines starting with `#` (comments and `##` directives) are skipped, as are lines with fewer than
9 tab-separated columns or with non-numeric start/end/score/phase fields. Column 6 (`score`) and
column 8 (`phase`) map to `null` when they carry the `.` undefined-field placeholder. The feature id
is taken from the column-9 `ID` attribute, falling back to `feature_{n}`. Column-9 attributes are
split on `;` into a key→value map with URL-decoded values.

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L439](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L439)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gff3Text` | string | Yes | GFF3 document text (lines separated by newlines) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `features` | array | `{ featureId, type, start, end, strand, score, phase, attributes }` per feature |

## Errors

| Code | Message |
|------|---------|
| 1001 | GFF3 text cannot be null or empty |

## Examples

### Example 1: Fully specified CDS line

`chr1  ENSEMBL  CDS  1000  2000  95.5  -  2  ID=cds1;Name=TestCDS` →

**Response:**
```json
{ "features": [ { "featureId": "cds1", "type": "CDS", "start": 1000, "end": 2000, "score": 95.5, "strand": "-", "phase": 2, "attributes": { "ID": "cds1", "Name": "TestCDS" } } ] }
```

### Example 2: Undefined score/phase

`chr1  src  gene  1  100  .  +  .  ID=gene1` yields `score: null` and `phase: null`.

## Performance

- **Time Complexity:** O(lines)
- **Space Complexity:** O(features)

## See Also

- [to_gff3](to_gff3.md) - Serialize annotations back to GFF3
- [predict_genes](predict_genes.md) - ORF-based gene prediction
