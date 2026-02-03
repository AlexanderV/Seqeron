# vcf_is_snp

Check if a variant is a SNP.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_is_snp` |
| **Method ID** | `VcfParser.IsSNP` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Determines if a variant is a Single Nucleotide Polymorphism (SNP). A SNP is defined as a variant where both the reference and alternate alleles are exactly 1 base pair long.

## Core Documentation Reference

- Source: [VcfParser.cs#L407](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/VcfParser.cs#L407)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `refAllele` | string | Yes | Reference allele |
| `altAllele` | string | Yes | Alternate allele |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `isSNP` | boolean | True if variant is a SNP |
| `refAllele` | string | Reference allele |
| `altAllele` | string | Alternate allele |

## Examples

### Example 1: SNP detected

**Expected Tool Call:**
```json
{
  "tool": "vcf_is_snp",
  "arguments": {
    "refAllele": "A",
    "altAllele": "G"
  }
}
```

**Response:**
```json
{
  "isSNP": true,
  "refAllele": "A",
  "altAllele": "G"
}
```

## See Also

- [vcf_classify](vcf_classify.md) - Full variant classification
- [vcf_is_indel](vcf_is_indel.md) - Check if indel
