# gff_statistics

Calculate statistics for GFF/GTF annotations.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `gff_statistics` |
| **Method ID** | `GffParser.CalculateStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates comprehensive statistics for GFF/GTF annotation files. Returns feature type counts, sequence IDs (chromosomes/contigs), annotation sources, and specific counts for genes and exons.

## Core Documentation Reference

- Source: [GffParser.cs#L377](../../../../Seqeron.Genomics/GffParser.cs#L377)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | GFF/GTF format content to analyze |
| `format` | string | No | Format: 'gff3', 'gtf', 'gff2', or 'auto' (default) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalFeatures` | integer | Total number of features |
| `featureTypeCounts` | object | Count by feature type |
| `sequenceIds` | array | List of unique sequence IDs |
| `sources` | array | List of unique annotation sources |
| `geneCount` | integer | Number of gene features |
| `exonCount` | integer | Number of exon features |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Get annotation statistics

**User Prompt:**
> What features are in this annotation file?

**Expected Tool Call:**
```json
{
  "tool": "gff_statistics",
  "arguments": {
    "content": "##gff-version 3\nchr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tID=gene1\nchr1\tHAVANA\texon\t100\t200\t.\t+\t.\tID=exon1\nchr1\tHAVANA\texon\t300\t500\t.\t+\t.\tID=exon2"
  }
}
```

**Response:**
```json
{
  "totalFeatures": 3,
  "featureTypeCounts": {
    "gene": 1,
    "exon": 2
  },
  "sequenceIds": ["chr1"],
  "sources": ["HAVANA"],
  "geneCount": 1,
  "exonCount": 2
}
```

## Performance

- **Time Complexity:** O(n) where n is the number of features
- **Space Complexity:** O(u) where u is the number of unique feature types/seqids/sources

## See Also

- [gff_parse](gff_parse.md) - Parse GFF/GTF format
- [gff_filter](gff_filter.md) - Filter GFF records
