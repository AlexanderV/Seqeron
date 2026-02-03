# vcf_is_het

Check if a genotype is heterozygous.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `vcf_is_het` |
| **Method ID** | `VcfParser.IsHet` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Determines if a genotype string represents a heterozygous call (e.g., 0/1, 0|1, 1/2). This indicates the sample has two different alleles.

## Core Documentation Reference

- Source: [VcfParser.cs#L554](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/VcfParser.cs#L554)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genotype` | string | Yes | Genotype string (e.g., '0/0', '0/1', '1/1') |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `result` | boolean | True if heterozygous |
| `genotype` | string | Input genotype string |
| `checkType` | string | "Heterozygous" |

## Examples

### Example 1: Heterozygous

**Expected Tool Call:**
```json
{
  "tool": "vcf_is_het",
  "arguments": {
    "genotype": "0/1"
  }
}
```

**Response:**
```json
{
  "result": true,
  "genotype": "0/1",
  "checkType": "Heterozygous"
}
```

## See Also

- [vcf_is_hom_ref](vcf_is_hom_ref.md) - Check homozygous reference
- [vcf_is_hom_alt](vcf_is_hom_alt.md) - Check homozygous alternate
