# vcf_is_hom_alt

Check if a genotype is homozygous alternate.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_is_hom_alt` |
| **Method ID** | `VcfParser.IsHomAlt` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Determines if a genotype string represents a homozygous alternate call (e.g., 1/1, 2/2, 1|1). This indicates the sample has two copies of an alternate allele.

## Core Documentation Reference

- Source: [VcfParser.cs#L539](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/VcfParser.cs#L539)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genotype` | string | Yes | Genotype string (e.g., '0/0', '0/1', '1/1') |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `result` | boolean | True if homozygous alternate |
| `genotype` | string | Input genotype string |
| `checkType` | string | "HomozygousAlternate" |

## Examples

### Example 1: Homozygous alternate

**Expected Tool Call:**
```json
{
  "tool": "vcf_is_hom_alt",
  "arguments": {
    "genotype": "1/1"
  }
}
```

**Response:**
```json
{
  "result": true,
  "genotype": "1/1",
  "checkType": "HomozygousAlternate"
}
```

## See Also

- [vcf_is_hom_ref](vcf_is_hom_ref.md) - Check homozygous reference
- [vcf_is_het](vcf_is_het.md) - Check heterozygous
