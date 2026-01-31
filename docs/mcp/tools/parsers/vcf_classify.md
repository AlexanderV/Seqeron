# vcf_classify

Classify variant type for a VCF record.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_classify` |
| **Method ID** | `VcfParser.ClassifyVariant` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Classifies a variant based on its reference and alternate alleles. Returns the variant type along with length information.

Variant types:
- **SNP**: Single Nucleotide Polymorphism (ref and alt both 1 bp)
- **MNP**: Multi-Nucleotide Polymorphism (ref and alt same length > 1 bp)
- **Insertion**: Alt longer than ref
- **Deletion**: Ref longer than alt
- **Complex**: Multiple indels or substitutions
- **Symbolic**: Structural variants (starts with <, [, or ])
- **Unknown**: Cannot classify

## Core Documentation Reference

- Source: [VcfParser.cs#L376](../../../../Seqeron.Genomics/VcfParser.cs#L376)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `refAllele` | string | Yes | Reference allele |
| `altAllele` | string | Yes | Alternate allele |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variantType` | string | Classified variant type |
| `refLength` | integer | Length of reference allele |
| `altLength` | integer | Length of alternate allele |
| `lengthDifference` | integer | Absolute difference between allele lengths |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference allele cannot be null or empty |
| 1002 | Alternate allele cannot be null or empty |

## Examples

### Example 1: Classify SNP

**Expected Tool Call:**
```json
{
  "tool": "vcf_classify",
  "arguments": {
    "refAllele": "A",
    "altAllele": "G"
  }
}
```

**Response:**
```json
{
  "variantType": "SNP",
  "refLength": 1,
  "altLength": 1,
  "lengthDifference": 0
}
```

### Example 2: Classify insertion

**Expected Tool Call:**
```json
{
  "tool": "vcf_classify",
  "arguments": {
    "refAllele": "A",
    "altAllele": "ATG"
  }
}
```

**Response:**
```json
{
  "variantType": "Insertion",
  "refLength": 1,
  "altLength": 3,
  "lengthDifference": 2
}
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [vcf_is_snp](vcf_is_snp.md) - Check if SNP
- [vcf_is_indel](vcf_is_indel.md) - Check if indel
- [vcf_variant_length](vcf_variant_length.md) - Get variant length
