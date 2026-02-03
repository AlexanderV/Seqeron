# fastq_encode_quality

Encode Phred quality scores to a quality string.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_encode_quality` |
| **Method ID** | `FastqParser.EncodeQualityScores` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Converts a list of Phred quality scores (integers) into an ASCII quality string suitable for FASTQ format. Scores are clamped to 0-41 range before encoding.

Encoding formula: `ASCII character = score + offset`
- Phred+33: offset = 33 (scores 0-41 map to '!' through 'J')
- Phred+64: offset = 64 (scores 0-41 map to '@' through 'i')

## Core Documentation Reference

- Source: [FastqParser.cs#L189](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/FastqParser.cs#L189)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `scores` | integer[] | Yes | - | List of Phred quality scores (0-41) |
| `encoding` | string | No | `"phred33"` | Quality encoding: `"phred33"` or `"phred64"` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `qualityString` | string | Encoded quality string |
| `length` | integer | Length of quality string |

## Errors

| Code | Message |
|------|---------|
| 1001 | Scores cannot be null or empty |
| 1002 | Invalid encoding. Use 'phred33' or 'phred64' |

## Examples

### Example 1: Encode to Phred+33

**Expected Tool Call:**
```json
{
  "tool": "fastq_encode_quality",
  "arguments": {
    "scores": [0, 10, 20, 30, 40],
    "encoding": "phred33"
  }
}
```

**Response:**
```json
{
  "qualityString": "!+5?I",
  "length": 5
}
```

### Example 2: Encode to Phred+64

**Expected Tool Call:**
```json
{
  "tool": "fastq_encode_quality",
  "arguments": {
    "scores": [0, 10, 20, 30, 40],
    "encoding": "phred64"
  }
}
```

**Response:**
```json
{
  "qualityString": "@JT^h",
  "length": 5
}
```

## Performance

- **Time Complexity:** O(n) where n is number of scores
- **Space Complexity:** O(n) for output string

## See Also

- [fastq_detect_encoding](fastq_detect_encoding.md) - Detect encoding
- [fastq_phred_to_error](fastq_phred_to_error.md) - Convert to error probability
