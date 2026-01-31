# genbank_extract_sequence

Extract a subsequence from a GenBank record based on a feature location string.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `genbank_extract_sequence` |
| **Method ID** | `GenBankParser.ExtractSequence` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Extracts a subsequence from a GenBank record using a feature location string. Handles various location formats:

- **Simple ranges:** `100..200`
- **Single positions:** `150`
- **Complement (reverse strand):** `complement(100..200)` - returns reverse complement
- **Join (multiple segments):** `join(100..200,300..400)` - concatenates segments
- **Complex locations:** `complement(join(100..200,300..400))`

Returns the extracted sequence along with metadata about the location.

## Core Documentation Reference

- Source: [GenBankParser.cs](../../../../Seqeron.Genomics/GenBankParser.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | GenBank format content |
| `locationString` | string | Yes | Feature location string (e.g., '100..200', 'complement(100..200)', 'join(100..200,300..400)') |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | string | Extracted sequence (uppercase) |
| `length` | integer | Length of extracted sequence |
| `isComplement` | boolean | True if sequence is from complement strand |
| `isJoin` | boolean | True if sequence spans multiple joined segments |
| `location` | string | Original location string |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |
| 1002 | Location string cannot be null or empty |
| 1003 | No GenBank records found in content |

## Examples

### Example 1: Extract simple range

**Expected Tool Call:**
```json
{
  "tool": "genbank_extract_sequence",
  "arguments": {
    "content": "LOCUS       TEST001...\nORIGIN\n        1 acgtacgtac...\n//",
    "locationString": "1..10"
  }
}
```

**Response:**
```json
{
  "sequence": "ACGTACGTAC",
  "length": 10,
  "isComplement": false,
  "isJoin": false,
  "location": "1..10"
}
```

### Example 2: Extract complement sequence

**Expected Tool Call:**
```json
{
  "tool": "genbank_extract_sequence",
  "arguments": {
    "content": "LOCUS       TEST001...\nORIGIN\n        1 acgtacgtac...\n//",
    "locationString": "complement(1..10)"
  }
}
```

**Response:**
```json
{
  "sequence": "GTACGTACGT",
  "length": 10,
  "isComplement": true,
  "isJoin": false,
  "location": "complement(1..10)"
}
```

### Example 3: Extract joined segments

**Expected Tool Call:**
```json
{
  "tool": "genbank_extract_sequence",
  "arguments": {
    "content": "LOCUS       TEST001...\nORIGIN\n        1 acgtacgtac gtacgtacgt...\n//",
    "locationString": "join(1..10,21..30)"
  }
}
```

**Response:**
```json
{
  "sequence": "ACGTACGTACACGTACGTAC",
  "length": 20,
  "isComplement": false,
  "isJoin": true,
  "location": "join(1..10,21..30)"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(m) where m is extracted sequence length

## See Also

- [genbank_parse](genbank_parse.md) - Parse full GenBank records
- [genbank_parse_location](genbank_parse_location.md) - Parse location strings
- [genbank_features](genbank_features.md) - Extract features
