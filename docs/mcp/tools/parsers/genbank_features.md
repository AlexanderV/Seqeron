# genbank_features

Extract features from GenBank records by feature type.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `genbank_features` |
| **Method ID** | `GenBankParser.GetFeatures` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Extracts sequence features from GenBank records. Can filter by feature type (gene, CDS, mRNA, exon, etc.). Returns feature locations, qualifiers, and counts by type.

## Core Documentation Reference

- Source: [GenBankParser.cs#L531](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/GenBankParser.cs#L531)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | GenBank format content to parse |
| `featureType` | string | No | Feature type to extract (e.g., 'gene', 'CDS', 'mRNA') |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `features` | array | List of extracted features |
| `count` | integer | Number of features |
| `featureTypeCounts` | object | Count by feature type |

### GenBankFeatureResult Schema

| Field | Type | Description |
|-------|------|-------------|
| `key` | string | Feature type (gene, CDS, etc.) |
| `start` | integer | Start position (1-based) |
| `end` | integer | End position (1-based) |
| `isComplement` | boolean | True if on complement strand |
| `isJoin` | boolean | True if location is a join |
| `rawLocation` | string | Original location string |
| `qualifiers` | object | Feature qualifiers (gene, product, etc.) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Extract all genes

**User Prompt:**
> Show me all genes in this GenBank file

**Expected Tool Call:**
```json
{
  "tool": "genbank_features",
  "arguments": {
    "content": "LOCUS...",
    "featureType": "gene"
  }
}
```

**Response:**
```json
{
  "features": [
    {
      "key": "gene",
      "start": 10,
      "end": 50,
      "isComplement": false,
      "isJoin": false,
      "rawLocation": "10..50",
      "qualifiers": { "gene": "testGene" }
    }
  ],
  "count": 1,
  "featureTypeCounts": { "gene": 1 }
}
```

### Example 2: Extract CDS features

**Expected Tool Call:**
```json
{
  "tool": "genbank_features",
  "arguments": {
    "content": "...",
    "featureType": "CDS"
  }
}
```

## Performance

- **Time Complexity:** O(n) where n is the number of features
- **Space Complexity:** O(m) where m is the number of matching features

## See Also

- [genbank_parse](genbank_parse.md) - Parse full GenBank records
- [genbank_statistics](genbank_statistics.md) - Calculate record statistics
