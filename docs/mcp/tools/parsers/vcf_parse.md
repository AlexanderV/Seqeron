# vcf_parse

Parse VCF (Variant Call Format) content into variant records.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_parse` |
| **Method ID** | `VcfParser.Parse` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses VCF format content into variant records. VCF is the standard format for storing genetic variation data including SNPs, insertions, deletions, and complex variants. Each record includes chromosome, position, reference/alternate alleles, quality score, filter status, and INFO fields.

## Core Documentation Reference

- Source: [VcfParser.cs#L116](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/VcfParser.cs#L116)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | VCF format content to parse |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of parsed VCF variant records |
| `records[].chrom` | string | Chromosome name |
| `records[].pos` | integer | 1-based position |
| `records[].id` | string | Variant identifier (e.g., rsID) |
| `records[].ref` | string | Reference allele |
| `records[].alt` | string[] | Alternate allele(s) |
| `records[].qual` | number? | Phred-scaled quality score |
| `records[].filter` | string[] | Filter status (PASS or filter names) |
| `records[].info` | object | INFO field key-value pairs |
| `records[].variantType` | string | Classified type: SNP, Insertion, Deletion, Complex |
| `count` | integer | Number of records parsed |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Parse VCF variants

**User Prompt:**
> Parse this VCF file with genetic variants

**Expected Tool Call:**
```json
{
  "tool": "vcf_parse",
  "arguments": {
    "content": "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\trs1\tA\tG\t30\tPASS\tDP=10"
  }
}
```

**Response:**
```json
{
  "records": [
    {
      "chrom": "chr1",
      "pos": 100,
      "id": "rs1",
      "ref": "A",
      "alt": ["G"],
      "qual": 30.0,
      "filter": ["PASS"],
      "info": { "DP": "10" },
      "variantType": "SNP"
    }
  ],
  "count": 1
}
```

## Variant Types

- **SNP**: Single nucleotide polymorphism (ref and alt same length = 1)
- **MNP**: Multiple nucleotide polymorphism (ref and alt same length > 1)
- **Insertion**: Alt longer than ref
- **Deletion**: Ref longer than alt
- **Complex**: Other structural changes

## See Also

- [vcf_statistics](vcf_statistics.md) - Calculate variant statistics
- [vcf_filter](vcf_filter.md) - Filter variants
