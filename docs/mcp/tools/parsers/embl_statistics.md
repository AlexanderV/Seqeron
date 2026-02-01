# embl_statistics

Calculate statistics for EMBL records.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `embl_statistics` |
| **Method ID** | `EmblParser.Statistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates comprehensive statistics for EMBL records. Returns counts of records, features by type, total sequence length, gene and CDS counts, and lists of molecule types and division codes.

## Core Documentation Reference

- Source: [EmblParser.cs](../../../../src/Seqeron/Seqeron.Genomics/EmblParser.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | EMBL format content to analyze |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `recordCount` | integer | Number of EMBL records |
| `totalFeatures` | integer | Total number of features |
| `totalSequenceLength` | integer | Total sequence length in bp |
| `featureTypeCounts` | object | Count by feature type |
| `moleculeTypes` | array | List of molecule types (DNA, RNA, etc.) |
| `divisions` | array | List of taxonomic division codes |
| `geneCount` | integer | Number of gene features |
| `cdsCount` | integer | Number of CDS features |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Get EMBL statistics

**Expected Tool Call:**
```json
{
  "tool": "embl_statistics",
  "arguments": {
    "content": "ID..."
  }
}
```

**Response:**
```json
{
  "recordCount": 1,
  "totalFeatures": 3,
  "totalSequenceLength": 100,
  "featureTypeCounts": { "gene": 1, "CDS": 1 },
  "moleculeTypes": ["DNA"],
  "divisions": ["HUM"],
  "geneCount": 1,
  "cdsCount": 1
}
```

## Performance

- **Time Complexity:** O(n) where n is the number of features
- **Space Complexity:** O(u) where u is the number of unique feature types

## See Also

- [embl_parse](embl_parse.md) - Parse full EMBL records
- [embl_features](embl_features.md) - Extract features by type
