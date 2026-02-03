# fasta_write

Write sequence to a FASTA file.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fasta_write` |
| **Method ID** | `FastaParser.WriteFile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Writes a DNA sequence to a FASTA format file. Creates a new file or overwrites an existing file at the specified path. Supports configurable line width for sequence wrapping.

## Core Documentation Reference

- Source: [FastaParser.cs#L98](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/FastaParser.cs#L98)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `filePath` | string | Yes | File path to write FASTA output |
| `id` | string | Yes | Sequence identifier |
| `sequence` | string | Yes | DNA sequence |
| `description` | string | No | Optional sequence description |
| `lineWidth` | integer | No | Line width for sequence wrapping (default: 80) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `filePath` | string | Path where file was written |
| `entriesWritten` | integer | Number of sequences written |
| `totalBases` | integer | Total number of bases written |

## Errors

| Code | Message |
|------|---------|
| 1001 | File path cannot be null or empty |
| 1002 | ID cannot be null or empty |
| 1003 | Sequence cannot be null or empty |
| 1004 | Line width must be at least 1 |
| 3001 | Cannot write to file (I/O error) |

## Examples

### Example 1: Write sequence to file

**User Prompt:**
> Save the sequence ATGCATGC with ID "seq1" to output.fasta

**Expected Tool Call:**
```json
{
  "tool": "fasta_write",
  "arguments": {
    "filePath": "/data/output.fasta",
    "id": "seq1",
    "sequence": "ATGCATGC",
    "description": "Human gene"
  }
}
```

**Response:**
```json
{
  "filePath": "/data/output.fasta",
  "entriesWritten": 1,
  "totalBases": 8
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(n) for output buffer

## See Also

- [fasta_parse](fasta_parse.md) - Parse FASTA content
- [fasta_format](fasta_format.md) - Format to FASTA string
