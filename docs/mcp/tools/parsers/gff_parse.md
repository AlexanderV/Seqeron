# gff_parse

Parse GFF3/GTF format content into feature records.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `gff_parse` |
| **Method ID** | `GffParser.Parse` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses GFF3/GTF format content into structured feature records. Supports GFF3, GTF (Gene Transfer Format), and GFF2 formats commonly used for gene annotations and genomic features. Returns sequence ID, source, feature type, coordinates, strand, phase, and all attributes.

## Core Documentation Reference

- Source: [GffParser.cs#L69](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/GffParser.cs#L69)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | GFF/GTF format content to parse |
| `format` | string | No | Format: 'gff3', 'gtf', 'gff2', or 'auto' (default) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of GFF feature records |
| `count` | integer | Number of records parsed |

### GffRecordResult Schema

| Field | Type | Description |
|-------|------|-------------|
| `seqid` | string | Sequence identifier (chromosome/contig) |
| `source` | string | Source of the feature annotation |
| `type` | string | Feature type (gene, exon, CDS, etc.) |
| `start` | integer | Start position (1-based) |
| `end` | integer | End position (1-based, inclusive) |
| `length` | integer | Feature length in base pairs |
| `score` | number | Feature score (nullable) |
| `strand` | string | Strand: '+', '-', or '.' |
| `phase` | integer | Reading frame phase for CDS (nullable) |
| `attributes` | object | Key-value pairs of attributes |
| `geneName` | string | Extracted gene name (nullable) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |
| 1002 | Invalid format |

## Examples

### Example 1: Parse GFF3 content

**User Prompt:**
> Parse this GFF3 annotation file

**Expected Tool Call:**
```json
{
  "tool": "gff_parse",
  "arguments": {
    "content": "##gff-version 3\nchr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=TestGene"
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
      "score": null,
      "strand": "+",
      "phase": null,
      "attributes": { "ID": "gene1", "Name": "TestGene" },
      "geneName": "TestGene"
    }
  ],
  "count": 1
}
```

### Example 2: Parse GTF format

**Expected Tool Call:**
```json
{
  "tool": "gff_parse",
  "arguments": {
    "content": "chr1\tensembl\texon\t100\t200\t.\t+\t.\tgene_id \"ENSG00001\"; transcript_id \"ENST00001\";",
    "format": "gtf"
  }
}
```

**Response:**
```json
{
  "records": [
    {
      "seqid": "chr1",
      "source": "ensembl",
      "type": "exon",
      "start": 100,
      "end": 200,
      "length": 101,
      "score": null,
      "strand": "+",
      "phase": null,
      "attributes": { "gene_id": "ENSG00001", "transcript_id": "ENST00001" },
      "geneName": "ENSG00001"
    }
  ],
  "count": 1
}
```

## Performance

- **Time Complexity:** O(n) where n is the number of lines
- **Space Complexity:** O(n) for storing parsed records

## See Also

- [gff_statistics](gff_statistics.md) - Calculate annotation statistics
- [gff_filter](gff_filter.md) - Filter GFF records
