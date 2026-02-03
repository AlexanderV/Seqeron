# fastq_write

Write FASTQ records to a file.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_write` |
| **Method ID** | `FastqParser.WriteToFile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Writes FASTQ records to a file at the specified path. The function parses the input content and writes it in standard FASTQ format. If the file already exists, it will be overwritten.

The output file is written in UTF-8 encoding with standard FASTQ format:
- Each record has 4 lines: header, sequence, separator (+), quality
- Records are separated by newlines

## Core Documentation Reference

- Source: [FastqParser.cs#L437](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/FastqParser.cs#L437)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `filePath` | string | Yes | - | File path to write FASTQ output |
| `content` | string | Yes | - | FASTQ format content to write |
| `encoding` | string | No | `"auto"` | Quality encoding for parsing: `"phred33"`, `"phred64"`, or `"auto"` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `filePath` | string | Path to the written file |
| `recordsWritten` | integer | Number of records written |
| `totalBases` | integer | Total bases written |

## Errors

| Code | Message |
|------|---------|
| 1001 | File path cannot be null or empty |
| 1002 | Content cannot be null or empty |
| 1003 | Invalid encoding. Use 'phred33', 'phred64', or 'auto' |

## Examples

### Example 1: Write FASTQ file

**Expected Tool Call:**
```json
{
  "tool": "fastq_write",
  "arguments": {
    "filePath": "/data/output.fastq",
    "content": "@read1\nATGCATGC\n+\nIIIIIIII\n@read2\nGGGCCC\n+\nHHHHHH"
  }
}
```

**Response:**
```json
{
  "filePath": "/data/output.fastq",
  "recordsWritten": 2,
  "totalBases": 14
}
```

### Example 2: Write with explicit encoding

**Expected Tool Call:**
```json
{
  "tool": "fastq_write",
  "arguments": {
    "filePath": "/data/output.fq",
    "content": "@read1\nATGC\n+\nIIII",
    "encoding": "phred33"
  }
}
```

**Response:**
```json
{
  "filePath": "/data/output.fq",
  "recordsWritten": 1,
  "totalBases": 4
}
```

## Performance

- **Time Complexity:** O(n) where n is total content length
- **Space Complexity:** O(n) for parsing records

## See Also

- [fastq_format](fastq_format.md) - Format to string (no file I/O)
- [fastq_parse](fastq_parse.md) - Parse FASTQ format
