# genbank_statistics

Calculate statistics for GenBank records.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `genbank_statistics` |
| **Method ID** | `GenBankParser.Statistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates comprehensive statistics for GenBank records. Returns counts of records, features by type, total sequence length, gene and CDS counts, and lists of molecule types and division codes.

## Core Documentation Reference

- Source: [GenBankParser.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/GenBankParser.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | GenBank format content to analyze |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `recordCount` | integer | Number of GenBank records |
| `totalFeatures` | integer | Total number of features |
| `totalSequenceLength` | integer | Total sequence length in bp |
| `featureTypeCounts` | object | Count by feature type |
| `moleculeTypes` | array | List of molecule types (DNA, RNA, etc.) |
| `divisions` | array | List of GenBank division codes |
| `geneCount` | integer | Number of gene features |
| `cdsCount` | integer | Number of CDS features |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Get GenBank statistics

**User Prompt:**
> What's in this GenBank file?

**Expected Tool Call:**
```json
{
  "tool": "genbank_statistics",
  "arguments": {
    "content": "LOCUS..."
  }
}
```

**Response:**
```json
{
  "recordCount": 1,
  "totalFeatures": 10,
  "totalSequenceLength": 5000,
  "featureTypeCounts": {
    "source": 1,
    "gene": 3,
    "CDS": 3,
    "mRNA": 3
  },
  "moleculeTypes": ["DNA"],
  "divisions": ["BCT"],
  "geneCount": 3,
  "cdsCount": 3
}
```

## Performance

- **Time Complexity:** O(n) where n is the number of features
- **Space Complexity:** O(u) where u is the number of unique feature types

## See Also

- [genbank_parse](genbank_parse.md) - Parse full GenBank records
- [genbank_features](genbank_features.md) - Extract features by type
