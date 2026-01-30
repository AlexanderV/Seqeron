# gff_filter

Filter GFF/GTF records by feature type, sequence ID, or genomic region.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `gff_filter` |
| **Method ID** | `GffParser.Filter*` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Filters GFF/GTF records using various criteria: feature type (gene, exon, CDS), sequence ID (chromosome), or genomic region coordinates. Multiple filters can be combined. Returns matching records with pass statistics.

## Core Documentation Reference

- Source: [GffParser.cs#L215](../../../../Seqeron.Genomics/GffParser.cs#L215)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | GFF/GTF format content to filter |
| `featureType` | string | No | Filter by feature type (e.g., 'gene', 'exon', 'CDS') |
| `seqid` | string | No | Filter by sequence ID (chromosome/contig name) |
| `regionStart` | integer | No | Filter by region start position (requires seqid and regionEnd) |
| `regionEnd` | integer | No | Filter by region end position (requires seqid and regionStart) |
| `format` | string | No | Format: 'gff3', 'gtf', 'gff2', or 'auto' (default) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of filtered GFF records |
| `passedCount` | integer | Number of records that passed filters |
| `totalCount` | integer | Total number of input records |
| `passedPercentage` | number | Percentage of records that passed |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Filter by feature type

**User Prompt:**
> Show me only the genes from this annotation

**Expected Tool Call:**
```json
{
  "tool": "gff_filter",
  "arguments": {
    "content": "##gff-version 3\nchr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tID=gene1\nchr1\tHAVANA\texon\t100\t200\t.\t+\t.\tID=exon1",
    "featureType": "gene"
  }
}
```

**Response:**
```json
{
  "records": [
    {
      "seqid": "chr1",
      "source": "HAVANA",
      "type": "gene",
      "start": 100,
      "end": 500,
      "length": 401,
      "strand": "+",
      "attributes": { "ID": "gene1" }
    }
  ],
  "passedCount": 1,
  "totalCount": 2,
  "passedPercentage": 50.0
}
```

### Example 2: Filter by chromosome

**Expected Tool Call:**
```json
{
  "tool": "gff_filter",
  "arguments": {
    "content": "...",
    "seqid": "chr1"
  }
}
```

### Example 3: Filter by region

**Expected Tool Call:**
```json
{
  "tool": "gff_filter",
  "arguments": {
    "content": "...",
    "seqid": "chr1",
    "regionStart": 1000,
    "regionEnd": 5000
  }
}
```

### Example 4: Combine filters

**User Prompt:**
> Show me all exons on chromosome 1

**Expected Tool Call:**
```json
{
  "tool": "gff_filter",
  "arguments": {
    "content": "...",
    "featureType": "exon",
    "seqid": "chr1"
  }
}
```

## Performance

- **Time Complexity:** O(n) where n is the number of records
- **Space Complexity:** O(m) where m is the number of matching records

## See Also

- [gff_parse](gff_parse.md) - Parse GFF/GTF format
- [gff_statistics](gff_statistics.md) - Calculate annotation statistics
