# variants_to_vcf

Format variants as VCF v4.2 lines (header + records).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `variants_to_vcf` |
| **Method ID** | `VariantCaller.ToVcfLines` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Serializes variants into VCF v4.2 text: three header lines (`##fileformat=VCFv4.2`,
`##source=Seqeron.Genomics`, and the `#CHROM…{sampleName}` column header) followed by one record per
variant. Each record uses the given `chromosome` for CHROM, `position + 1` for the 1-based POS, `.`
for ID/QUAL/INFO, `PASS` for FILTER, `GT` for FORMAT and `0/1` for the sample genotype. Gap alleles
(`-`) are written as `.` in REF/ALT.

## Core Documentation Reference

- Source: [VariantCaller.cs#L350](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L350)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `variants` | array | Yes | — | Variants to serialize |
| `chromosome` | string | No | chr1 | CHROM value for column 1 |
| `sampleName` | string | No | SAMPLE | Sample column name in the header |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `lines` | array | VCF lines: 3 header lines then one record per variant |

## Errors

| Code | Message |
|------|---------|
| 1001 | Variants cannot be null |

## Examples

### Example 1: Single SNP

A `G>T` SNP at 0-based position 2 → VCF POS 3.

**Response:**
```json
{ "lines": [ "##fileformat=VCFv4.2", "##source=Seqeron.Genomics", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE", "chr1\t3\t.\tG\tT\t.\tPASS\t.\tGT\t0/1" ] }
```

### Example 2: Insertion

An insertion carries `.` in REF (gap allele) and POS = position + 1.

## Performance

- **Time Complexity:** O(variants)
- **Space Complexity:** O(variants)

## See Also

- [call_variants](call_variants.md) - Produce the variants to serialize
- [parse_vcf_variant](parse_vcf_variant.md) - Build a variant from VCF fields
