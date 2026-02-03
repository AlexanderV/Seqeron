# genbank_parse

Parse GenBank flat file format into structured records.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `genbank_parse` |
| **Method ID** | `GenBankParser.Parse` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses GenBank flat file format (.gb, .gbk) into structured records. Returns locus information, definition, accession, version, keywords, organism taxonomy, features with qualifiers, and the nucleotide sequence.

## Core Documentation Reference

- Source: [GenBankParser.cs#L88](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/GenBankParser.cs#L88)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | GenBank format content to parse |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of GenBank records |
| `count` | integer | Number of records parsed |

### GenBankRecordResult Schema

| Field | Type | Description |
|-------|------|-------------|
| `locus` | string | Locus name/identifier |
| `sequenceLength` | integer | Declared sequence length |
| `moleculeType` | string | Molecule type (DNA, RNA, etc.) |
| `topology` | string | Topology (linear, circular) |
| `division` | string | GenBank division code |
| `date` | string | Modification date (nullable) |
| `definition` | string | Sequence definition |
| `accession` | string | Accession number |
| `version` | string | Version identifier |
| `keywords` | array | List of keywords |
| `organism` | string | Source organism name |
| `taxonomy` | string | Taxonomic classification |
| `features` | array | List of sequence features |
| `actualSequenceLength` | integer | Actual parsed sequence length |
| `sequence` | string | Nucleotide sequence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Parse GenBank record

**User Prompt:**
> Parse this GenBank file

**Expected Tool Call:**
```json
{
  "tool": "genbank_parse",
  "arguments": {
    "content": "LOCUS       TEST001    100 bp    DNA     linear   BCT 01-JAN-2024\nDEFINITION  Test sequence.\n..."
  }
}
```

**Response:**
```json
{
  "records": [
    {
      "locus": "TEST001",
      "sequenceLength": 100,
      "moleculeType": "DNA",
      "topology": "linear",
      "division": "BCT",
      "definition": "Test sequence.",
      "accession": "TEST001",
      "features": [...],
      "sequence": "ATGGTGCTGT..."
    }
  ],
  "count": 1
}
```

## Performance

- **Time Complexity:** O(n) where n is the content size
- **Space Complexity:** O(n) for storing parsed records

## See Also

- [genbank_features](genbank_features.md) - Extract features by type
- [genbank_statistics](genbank_statistics.md) - Calculate record statistics
