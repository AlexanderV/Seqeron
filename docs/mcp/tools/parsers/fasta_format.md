# fasta_format

Format sequence to FASTA string.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fasta_format` |
| **Method ID** | `FastaParser.ToFasta` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Formats a DNA sequence into FASTA format string. Takes a sequence ID, the DNA sequence, and an optional description. Supports configurable line width for sequence wrapping.

## Core Documentation Reference

- Source: [FastaParser.cs#L78](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/FastaParser.cs#L78)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `id` | string | Yes | Sequence identifier |
| `sequence` | string | Yes | DNA sequence |
| `description` | string | No | Optional sequence description |
| `lineWidth` | integer | No | Line width for sequence wrapping (default: 80) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fasta` | string | Formatted FASTA string |

## Errors

| Code | Message |
|------|---------|
| 1001 | ID cannot be null or empty |
| 1002 | Sequence cannot be null or empty |
| 1003 | Line width must be at least 1 |

## Examples

### Example 1: Format with description

**User Prompt:**
> Create a FASTA entry for sequence "ATGCATGC" with ID "seq1" and description "Human gene"

**Expected Tool Call:**
```json
{
  "tool": "fasta_format",
  "arguments": {
    "id": "seq1",
    "sequence": "ATGCATGC",
    "description": "Human gene"
  }
}
```

**Response:**
```json
{
  "fasta": ">seq1 Human gene\nATGCATGC"
}
```

### Example 2: Format without description

**User Prompt:**
> Convert sequence GGGGCCCC with ID NM_001 to FASTA format

**Expected Tool Call:**
```json
{
  "tool": "fasta_format",
  "arguments": {
    "id": "NM_001",
    "sequence": "GGGGCCCC"
  }
}
```

**Response:**
```json
{
  "fasta": ">NM_001\nGGGGCCCC"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(n) for output string

## See Also

- [fasta_parse](fasta_parse.md) - Parse FASTA content
- [fasta_write](fasta_write.md) - Write FASTA to file
