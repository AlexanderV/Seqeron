# fastq_format

Format a single FASTQ record to string format.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_format` |
| **Method ID** | `FastqParser.ToFastqString` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Converts sequence data into a properly formatted FASTQ string. The output follows the standard FASTQ format with four lines per record:

1. Header line starting with '@' followed by ID and optional description
2. Sequence line (nucleotide sequence)
3. Separator line with '+'
4. Quality line (same length as sequence)

## Core Documentation Reference

- Source: [FastqParser.cs#L474](../../../../src/Seqeron/Seqeron.Genomics/FastqParser.cs#L474)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `id` | string | Yes | Sequence identifier |
| `sequence` | string | Yes | DNA sequence |
| `qualityString` | string | Yes | Quality string (same length as sequence) |
| `description` | string | No | Optional sequence description |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fastq` | string | Formatted FASTQ string |

## Errors

| Code | Message |
|------|---------|
| 1001 | ID cannot be null or empty |
| 1002 | Sequence cannot be null or empty |
| 1003 | Quality string cannot be null or empty |
| 1004 | Sequence and quality string must have the same length |

## Examples

### Example 1: Format with description

**Expected Tool Call:**
```json
{
  "tool": "fastq_format",
  "arguments": {
    "id": "read001",
    "sequence": "ATGCATGC",
    "qualityString": "IIIIIIII",
    "description": "Sample read"
  }
}
```

**Response:**
```json
{
  "fastq": "@read001 Sample read\nATGCATGC\n+\nIIIIIIII"
}
```

### Example 2: Format without description

**Expected Tool Call:**
```json
{
  "tool": "fastq_format",
  "arguments": {
    "id": "read001",
    "sequence": "ATGC",
    "qualityString": "IIII"
  }
}
```

**Response:**
```json
{
  "fastq": "@read001\nATGC\n+\nIIII"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(n) for output string

## See Also

- [fastq_write](fastq_write.md) - Write FASTQ to file
- [fastq_parse](fastq_parse.md) - Parse FASTQ format
