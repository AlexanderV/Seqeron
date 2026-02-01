# vcf_has_flag

Check if a VCF INFO field flag is present.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_has_flag` |
| **Method ID** | `VcfParser.HasInfoFlag` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Checks if a specific INFO field flag is present in VCF records. Returns statistics about how many records contain the flag.

Common VCF INFO flags:
- **DB**: Variant found in dbSNP
- **H2**: Variant found in HapMap2
- **H3**: Variant found in HapMap3
- **SOMATIC**: Variant is somatic mutation

## Core Documentation Reference

- Source: [VcfParser.cs#L835](../../../../src/Seqeron/Seqeron.Genomics/VcfParser.cs#L835)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | VCF format content |
| `flag` | string | Yes | INFO field flag name to check |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `flag` | string | Flag name that was checked |
| `recordsWithFlag` | integer | Number of records with the flag |
| `totalRecords` | integer | Total number of records |
| `percentage` | number | Percentage of records with flag |

## Examples

### Example 1: Check for dbSNP flag

**Expected Tool Call:**
```json
{
  "tool": "vcf_has_flag",
  "arguments": {
    "content": "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t30\tPASS\tDB\nchr1\t200\t.\tC\tT\t40\tPASS\t.",
    "flag": "DB"
  }
}
```

**Response:**
```json
{
  "flag": "DB",
  "recordsWithFlag": 1,
  "totalRecords": 2,
  "percentage": 50.0
}
```

## See Also

- [vcf_parse](vcf_parse.md) - Parse VCF format
- [vcf_filter](vcf_filter.md) - Filter VCF records
