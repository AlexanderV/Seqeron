# vcf_write

Write VCF records to a file.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_write` |
| **Method ID** | `VcfParser.WriteToFile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Writes VCF records to a file at the specified path. Preserves header metadata and creates properly formatted VCF output. If the file exists, it will be overwritten.

## Core Documentation Reference

- Source: [VcfParser.cs#L702](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/VcfParser.cs#L702)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `filePath` | string | Yes | File path to write VCF output |
| `content` | string | Yes | VCF format content to write |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `filePath` | string | Path to the written file |
| `recordsWritten` | integer | Number of records written |

## Examples

### Example 1: Write VCF file

**Expected Tool Call:**
```json
{
  "tool": "vcf_write",
  "arguments": {
    "filePath": "/data/output.vcf",
    "content": "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t30\tPASS\t."
  }
}
```

**Response:**
```json
{
  "filePath": "/data/output.vcf",
  "recordsWritten": 1
}
```

## See Also

- [vcf_parse](vcf_parse.md) - Parse VCF format
- [vcf_filter](vcf_filter.md) - Filter VCF records
