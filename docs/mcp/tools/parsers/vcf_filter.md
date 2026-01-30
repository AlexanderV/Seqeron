# vcf_filter

Filter VCF variants by type, quality, chromosome, or PASS status.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_filter` |
| **Method ID** | `VcfParser.Filter*` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Filters VCF variants using multiple optional criteria. All filters can be combined. Useful for selecting specific variant types, high-quality calls, or variants from specific chromosomes.

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `content` | string | Yes | - | VCF format content to filter |
| `variantType` | string | No | - | Filter by type: 'snp', 'indel', 'insertion', 'deletion', 'complex' |
| `chrom` | string | No | - | Filter by chromosome name |
| `minQuality` | number | No | - | Minimum quality score |
| `passOnly` | boolean | No | false | Only include PASS variants |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of filtered VCF records |
| `passedCount` | integer | Number of records that passed filters |
| `totalCount` | integer | Total number of input records |
| `passedPercentage` | number | Percentage of records that passed |

## Examples

### Example 1: Filter SNPs only

**Expected Tool Call:**
```json
{
  "tool": "vcf_filter",
  "arguments": {
    "content": "...",
    "variantType": "snp"
  }
}
```

### Example 2: High-quality PASS variants

**Expected Tool Call:**
```json
{
  "tool": "vcf_filter",
  "arguments": {
    "content": "...",
    "minQuality": 30,
    "passOnly": true
  }
}
```

### Example 3: Combine filters

**Expected Tool Call:**
```json
{
  "tool": "vcf_filter",
  "arguments": {
    "content": "...",
    "variantType": "snp",
    "chrom": "chr1",
    "passOnly": true
  }
}
```

## See Also

- [vcf_parse](vcf_parse.md) - Parse VCF content
- [vcf_statistics](vcf_statistics.md) - Calculate variant statistics
