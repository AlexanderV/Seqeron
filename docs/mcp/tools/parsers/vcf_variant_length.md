# vcf_variant_length

Get the length difference of a variant.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_variant_length` |
| **Method ID** | `VcfParser.GetVariantLength` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the absolute length difference between reference and alternate alleles. For SNPs this returns 0. For indels, it returns the size of the insertion or deletion.

## Core Documentation Reference

- Source: [VcfParser.cs#L421](../../../../Seqeron.Genomics/VcfParser.cs#L421)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `refAllele` | string | Yes | Reference allele |
| `altAllele` | string | Yes | Alternate allele |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `length` | integer | Absolute length difference |
| `refLength` | integer | Length of reference allele |
| `altLength` | integer | Length of alternate allele |

## Examples

### Example 1: Insertion length

**Expected Tool Call:**
```json
{
  "tool": "vcf_variant_length",
  "arguments": {
    "refAllele": "A",
    "altAllele": "ATGC"
  }
}
```

**Response:**
```json
{
  "length": 3,
  "refLength": 1,
  "altLength": 4
}
```

### Example 2: SNP (zero length)

**Expected Tool Call:**
```json
{
  "tool": "vcf_variant_length",
  "arguments": {
    "refAllele": "A",
    "altAllele": "G"
  }
}
```

**Response:**
```json
{
  "length": 0,
  "refLength": 1,
  "altLength": 1
}
```

## See Also

- [vcf_classify](vcf_classify.md) - Full variant classification
- [vcf_is_indel](vcf_is_indel.md) - Check if indel
