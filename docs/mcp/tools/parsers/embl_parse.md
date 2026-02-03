# embl_parse

Parse EMBL flat file format into structured records.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `embl_parse` |
| **Method ID** | `EmblParser.Parse` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses EMBL flat file format (.embl, .dat) into structured records. EMBL format uses two-letter line type prefixes. Returns accession, description, organism, keywords, features with qualifiers, and the nucleotide sequence.

## Core Documentation Reference

- Source: [EmblParser.cs#L114](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/EmblParser.cs#L114)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | EMBL format content to parse |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of EMBL records |
| `count` | integer | Number of records parsed |

### EmblRecordResult Schema

| Field | Type | Description |
|-------|------|-------------|
| `accession` | string | Accession number |
| `sequenceVersion` | string | Sequence version |
| `dataClass` | string | Data class (STD, CON, etc.) |
| `moleculeType` | string | Molecule type (DNA, RNA, etc.) |
| `topology` | string | Topology (linear, circular) |
| `taxonomicDivision` | string | Taxonomic division code |
| `sequenceLength` | integer | Declared sequence length |
| `description` | string | Sequence description |
| `keywords` | array | List of keywords |
| `organism` | string | Source organism name |
| `organismClassification` | array | Taxonomic classification |
| `features` | array | List of sequence features |
| `actualSequenceLength` | integer | Actual parsed sequence length |
| `sequence` | string | Nucleotide sequence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Parse EMBL record

**User Prompt:**
> Parse this EMBL file

**Expected Tool Call:**
```json
{
  "tool": "embl_parse",
  "arguments": {
    "content": "ID   TEST001; SV 1; linear; DNA; STD; HUM; 100 BP.\nXX\nAC   TEST001;\n..."
  }
}
```

**Response:**
```json
{
  "records": [
    {
      "accession": "TEST001",
      "sequenceLength": 100,
      "moleculeType": "DNA",
      "topology": "linear",
      "description": "Test sequence.",
      "features": [...],
      "sequence": "ACGTACGT..."
    }
  ],
  "count": 1
}
```

## Performance

- **Time Complexity:** O(n) where n is the content size
- **Space Complexity:** O(n) for storing parsed records

## See Also

- [embl_features](embl_features.md) - Extract features by type
- [embl_statistics](embl_statistics.md) - Calculate record statistics
- [genbank_parse](genbank_parse.md) - Parse GenBank format (similar structure)
