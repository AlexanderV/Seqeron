# vcf_statistics

Calculate statistics for VCF variants.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_statistics` |
| **Method ID** | `VcfParser.CalculateStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates comprehensive statistics for VCF variant data including counts by variant type (SNP, indel, complex), chromosome distribution, passing variants, and quality metrics.

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | VCF format content to analyze |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalVariants` | integer | Total number of variants |
| `snpCount` | integer | Number of SNPs |
| `indelCount` | integer | Number of insertions + deletions |
| `complexCount` | integer | Number of complex variants |
| `passingCount` | integer | Number of variants with PASS filter |
| `chromosomeCounts` | object | Variant counts per chromosome |
| `meanQuality` | number? | Average quality score |

## Examples

### Example 1: Calculate variant statistics

**Expected Tool Call:**
```json
{
  "tool": "vcf_statistics",
  "arguments": {
    "content": "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t30\tPASS\t.\nchr1\t200\t.\tAT\tA\t25\tPASS\t.\nchr2\t300\t.\tC\tT\t20\t.\t."
  }
}
```

**Response:**
```json
{
  "totalVariants": 3,
  "snpCount": 2,
  "indelCount": 1,
  "complexCount": 0,
  "passingCount": 2,
  "chromosomeCounts": { "chr1": 2, "chr2": 1 },
  "meanQuality": 25.0
}
```

## See Also

- [vcf_parse](vcf_parse.md) - Parse VCF content
- [vcf_filter](vcf_filter.md) - Filter variants
