# fastq_detect_encoding

Detect quality score encoding from FASTQ quality string.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_detect_encoding` |
| **Method ID** | `FastqParser.DetectEncoding` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Automatically detects the quality score encoding used in a FASTQ quality string. Returns either 'Phred33' (Sanger/Illumina 1.8+) or 'Phred64' (Illumina 1.3-1.7) based on the ASCII characters present.

Detection logic:
- Characters below '@' (ASCII 64) indicate Phred+33
- Characters above 'I' (ASCII 73) indicate Phred+64
- Ambiguous ranges default to Phred+33 (most common modern format)

## Core Documentation Reference

- Source: [FastqParser.cs#L148](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/FastqParser.cs#L148)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `qualityString` | string | Yes | Quality string from FASTQ record |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `encoding` | string | Detected encoding: "Phred33" or "Phred64" |
| `offset` | integer | ASCII offset: 33 or 64 |

## Errors

| Code | Message |
|------|---------|
| 1001 | Quality string cannot be null or empty |

## Examples

### Example 1: Detect Phred+33 encoding

**Expected Tool Call:**
```json
{
  "tool": "fastq_detect_encoding",
  "arguments": {
    "qualityString": "!\"#$%&'()*+,-./0123456789:"
  }
}
```

**Response:**
```json
{
  "encoding": "Phred33",
  "offset": 33
}
```

### Example 2: Detect Phred+64 encoding

**Expected Tool Call:**
```json
{
  "tool": "fastq_detect_encoding",
  "arguments": {
    "qualityString": "efghijklmnopqrstuvwxyz"
  }
}
```

**Response:**
```json
{
  "encoding": "Phred64",
  "offset": 64
}
```

## Performance

- **Time Complexity:** O(n) where n is quality string length
- **Space Complexity:** O(1)

## See Also

- [fastq_parse](fastq_parse.md) - Parse FASTQ format
- [fastq_encode_quality](fastq_encode_quality.md) - Encode quality scores
