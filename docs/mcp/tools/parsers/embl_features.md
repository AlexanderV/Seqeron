# embl_features

Extract features from EMBL records by feature type.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `embl_features` |
| **Method ID** | `EmblParser.GetFeatures` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Extracts sequence features from EMBL records. Can filter by feature type (gene, CDS, mRNA, exon, etc.). Returns feature locations, qualifiers, and counts by type.

## Core Documentation Reference

- Source: [EmblParser.cs#L575](../../../../src/Seqeron/Seqeron.Genomics/EmblParser.cs#L575)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | EMBL format content to parse |
| `featureType` | string | No | Feature type to extract (e.g., 'gene', 'CDS', 'mRNA') |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `features` | array | List of extracted features |
| `count` | integer | Number of features |
| `featureTypeCounts` | object | Count by feature type |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Extract all genes

**Expected Tool Call:**
```json
{
  "tool": "embl_features",
  "arguments": {
    "content": "ID...",
    "featureType": "gene"
  }
}
```

## Performance

- **Time Complexity:** O(n) where n is the number of features
- **Space Complexity:** O(m) where m is the number of matching features

## See Also

- [embl_parse](embl_parse.md) - Parse full EMBL records
- [embl_statistics](embl_statistics.md) - Calculate record statistics
