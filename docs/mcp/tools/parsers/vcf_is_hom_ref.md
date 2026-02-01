# vcf_is_hom_ref

Check if a genotype is homozygous reference.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_is_hom_ref` |
| **Method ID** | `VcfParser.IsHomRef` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Determines if a genotype string represents a homozygous reference call (0/0 or 0|0). This indicates the sample has two copies of the reference allele.

## Core Documentation Reference

- Source: [VcfParser.cs#L529](../../../../src/Seqeron/Seqeron.Genomics/VcfParser.cs#L529)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genotype` | string | Yes | Genotype string (e.g., '0/0', '0/1', '1/1') |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `result` | boolean | True if homozygous reference |
| `genotype` | string | Input genotype string |
| `checkType` | string | "HomozygousReference" |

## Examples

### Example 1: Homozygous reference

**Expected Tool Call:**
```json
{
  "tool": "vcf_is_hom_ref",
  "arguments": {
    "genotype": "0/0"
  }
}
```

**Response:**
```json
{
  "result": true,
  "genotype": "0/0",
  "checkType": "HomozygousReference"
}
```

## See Also

- [vcf_is_hom_alt](vcf_is_hom_alt.md) - Check homozygous alternate
- [vcf_is_het](vcf_is_het.md) - Check heterozygous
