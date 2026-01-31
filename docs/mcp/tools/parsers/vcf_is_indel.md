# vcf_is_indel

Check if a variant is an indel.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_is_indel` |
| **Method ID** | `VcfParser.IsIndel` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Determines if a variant is an insertion or deletion (indel). Returns detailed information about whether the variant is specifically an insertion or deletion.

## Core Documentation Reference

- Source: [VcfParser.cs#L412](../../../../Seqeron.Genomics/VcfParser.cs#L412)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `refAllele` | string | Yes | Reference allele |
| `altAllele` | string | Yes | Alternate allele |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `isIndel` | boolean | True if variant is an indel |
| `isInsertion` | boolean | True if variant is an insertion |
| `isDeletion` | boolean | True if variant is a deletion |
| `refAllele` | string | Reference allele |
| `altAllele` | string | Alternate allele |

## Examples

### Example 1: Insertion detected

**Expected Tool Call:**
```json
{
  "tool": "vcf_is_indel",
  "arguments": {
    "refAllele": "A",
    "altAllele": "ATG"
  }
}
```

**Response:**
```json
{
  "isIndel": true,
  "isInsertion": true,
  "isDeletion": false,
  "refAllele": "A",
  "altAllele": "ATG"
}
```

### Example 2: Deletion detected

**Expected Tool Call:**
```json
{
  "tool": "vcf_is_indel",
  "arguments": {
    "refAllele": "ATG",
    "altAllele": "A"
  }
}
```

**Response:**
```json
{
  "isIndel": true,
  "isInsertion": false,
  "isDeletion": true,
  "refAllele": "ATG",
  "altAllele": "A"
}
```

## See Also

- [vcf_classify](vcf_classify.md) - Full variant classification
- [vcf_is_snp](vcf_is_snp.md) - Check if SNP
- [vcf_variant_length](vcf_variant_length.md) - Get variant length
