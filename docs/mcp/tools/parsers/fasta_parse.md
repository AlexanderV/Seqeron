# fasta_parse

Parse FASTA format content into sequence entries.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fasta_parse` |
| **Method ID** | `FastaParser.Parse` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses FASTA format string into individual sequence entries. Each entry includes the sequence ID (first word of header), optional description (rest of header), the DNA sequence, and its length. FASTA is the most common format for biological sequences.

## Core Documentation Reference

- Source: [FastaParser.cs#L17](../../../../src/Seqeron/Seqeron.Genomics/FastaParser.cs#L17)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | FASTA format content to parse |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entries` | array | List of parsed FASTA entries |
| `entries[].id` | string | Sequence identifier (first word of header) |
| `entries[].description` | string? | Optional description (rest of header) |
| `entries[].sequence` | string | DNA sequence |
| `entries[].length` | integer | Sequence length in base pairs |
| `count` | integer | Number of entries parsed |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Parse multiple sequences

**User Prompt:**
> Parse this FASTA file content with two sequences

**Expected Tool Call:**
```json
{
  "tool": "fasta_parse",
  "arguments": {
    "content": ">seq1 Human gene\nATGCATGC\n>seq2 Mouse gene\nGGGCCC"
  }
}
```

**Response:**
```json
{
  "entries": [
    { "id": "seq1", "description": "Human gene", "sequence": "ATGCATGC", "length": 8 },
    { "id": "seq2", "description": "Mouse gene", "sequence": "GGGCCC", "length": 6 }
  ],
  "count": 2
}
```

### Example 2: Parse single sequence without description

**User Prompt:**
> Parse this FASTA sequence: >NM_001 followed by ATGCGATCGATCG

**Expected Tool Call:**
```json
{
  "tool": "fasta_parse",
  "arguments": {
    "content": ">NM_001\nATGCGATCGATCG"
  }
}
```

**Response:**
```json
{
  "entries": [
    { "id": "NM_001", "description": null, "sequence": "ATGCGATCGATCG", "length": 13 }
  ],
  "count": 1
}
```

## Performance

- **Time Complexity:** O(n) where n is total content length
- **Space Complexity:** O(n) for storing parsed entries

## See Also

- [fasta_format](fasta_format.md) - Format sequences to FASTA string
- [fasta_write](fasta_write.md) - Write sequences to FASTA file
- [fastq_parse](fastq_parse.md) - Parse FASTQ format
