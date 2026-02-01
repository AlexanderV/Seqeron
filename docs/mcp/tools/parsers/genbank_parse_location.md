# genbank_parse_location

Parse a GenBank feature location string into its components.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `genbank_parse_location` |
| **Method ID** | `GenBankParser.ParseLocation` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses GenBank/EMBL feature location strings into structured components. Handles various location formats:

- **Simple ranges:** `100..200`
- **Single positions:** `150`
- **Complement (reverse strand):** `complement(100..200)`
- **Join (multiple segments):** `join(100..200,300..400)`
- **Complex locations:** `complement(join(100..200,300..400))`

Returns the overall span, individual parts, and whether the location is on the complement strand or represents joined segments.

## Core Documentation Reference

- Source: [GenBankParser.cs#L471](../../../../src/Seqeron/Seqeron.Genomics/GenBankParser.cs#L471)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `locationString` | string | Yes | Feature location string (e.g., '100..200', 'complement(100..200)') |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `start` | integer | Overall start position (minimum of all parts) |
| `end` | integer | Overall end position (maximum of all parts) |
| `length` | integer | Total span length (end - start + 1) |
| `isComplement` | boolean | True if on complement (reverse) strand |
| `isJoin` | boolean | True if location has multiple joined segments |
| `parts` | array | Individual location segments |
| `parts[].start` | integer | Segment start position |
| `parts[].end` | integer | Segment end position |
| `parts[].length` | integer | Segment length |
| `rawLocation` | string | Original location string |

## Errors

| Code | Message |
|------|---------|
| 1001 | Location string cannot be null or empty |

## Examples

### Example 1: Simple range

**Expected Tool Call:**
```json
{
  "tool": "genbank_parse_location",
  "arguments": {
    "locationString": "100..200"
  }
}
```

**Response:**
```json
{
  "start": 100,
  "end": 200,
  "length": 101,
  "isComplement": false,
  "isJoin": false,
  "parts": [{ "start": 100, "end": 200, "length": 101 }],
  "rawLocation": "100..200"
}
```

### Example 2: Complement location

**Expected Tool Call:**
```json
{
  "tool": "genbank_parse_location",
  "arguments": {
    "locationString": "complement(100..200)"
  }
}
```

**Response:**
```json
{
  "start": 100,
  "end": 200,
  "length": 101,
  "isComplement": true,
  "isJoin": false,
  "parts": [{ "start": 100, "end": 200, "length": 101 }],
  "rawLocation": "complement(100..200)"
}
```

### Example 3: Join location (multi-exon)

**Expected Tool Call:**
```json
{
  "tool": "genbank_parse_location",
  "arguments": {
    "locationString": "join(100..200,300..400)"
  }
}
```

**Response:**
```json
{
  "start": 100,
  "end": 400,
  "length": 301,
  "isComplement": false,
  "isJoin": true,
  "parts": [
    { "start": 100, "end": 200, "length": 101 },
    { "start": 300, "end": 400, "length": 101 }
  ],
  "rawLocation": "join(100..200,300..400)"
}
```

## Performance

- **Time Complexity:** O(n) where n is location string length
- **Space Complexity:** O(p) where p is number of parts

## See Also

- [genbank_parse](genbank_parse.md) - Parse full GenBank records
- [genbank_features](genbank_features.md) - Extract features
- [genbank_extract_sequence](genbank_extract_sequence.md) - Extract sequence by location
