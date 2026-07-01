# parse_vcf_variant

Build a `Variant` record from VCF fields and classify its type.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `parse_vcf_variant` |
| **Method ID** | `VariantAnnotator.ParseVcfVariant` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assembles a variant record from the core VCF columns (CHROM, POS, ID, REF, ALT, QUAL), preserving each
field verbatim, and assigns a `type` by delegating to `VariantAnnotator.ClassifyVariant`. Classification
is purely allele-length based:

- 1 bp REF and 1 bp ALT ⇒ `SNV`
- equal-length REF/ALT longer than 1 bp ⇒ `MNV`
- 1 bp REF where ALT starts with REF ⇒ `Insertion`
- 1 bp ALT where REF starts with ALT ⇒ `Deletion`
- any other length mismatch ⇒ `Indel`
- empty REF or ALT ⇒ `Complex`

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L1469](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L1469)
- Classification: [VariantAnnotator.cs#L196](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L196)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chromosome` | string | Yes | Chromosome / contig (CHROM) |
| `position` | integer | Yes | 1-based position (POS) |
| `id` | string | Yes | Variant ID (ID column; `.` if none) |
| `reference` | string | Yes | Reference allele (REF) |
| `alternate` | string | Yes | Alternate allele (ALT) |
| `quality` | number \| null | No | Quality score (QUAL); null if unknown |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variant.chromosome` | string | Chromosome |
| `variant.position` | integer | 1-based position |
| `variant.reference` | string | Reference allele |
| `variant.alternate` | string | Alternate allele |
| `variant.type` | string | SNV, Insertion, Deletion, MNV, Indel, or Complex |
| `variant.quality` | number \| null | Quality score |
| `variant.id` | string \| null | Variant ID |

## Errors

| Code | Message |
|------|---------|
| 1001 | Chromosome cannot be null or empty |
| 1001 | Reference allele cannot be null or empty |
| 1001 | Alternate allele cannot be null or empty |

## Examples

### Example 1: Substitution classified as SNV

**User Prompt:**
> Parse the VCF line chr1 100 rs123 A G 42.

**Expected Tool Call:**
```json
{
  "tool": "parse_vcf_variant",
  "arguments": { "chromosome": "chr1", "position": 100, "id": "rs123", "reference": "A", "alternate": "G", "quality": 42.0 }
}
```

**Response:**
```json
{
  "variant": { "chromosome": "chr1", "position": 100, "reference": "A", "alternate": "G", "type": "SNV", "quality": 42.0, "id": "rs123" }
}
```

### Example 2: Insertion

**Expected Tool Call:**
```json
{
  "tool": "parse_vcf_variant",
  "arguments": { "chromosome": "chr1", "position": 100, "id": ".", "reference": "A", "alternate": "AT" }
}
```

**Response:**
```json
{
  "variant": { "chromosome": "chr1", "position": 100, "reference": "A", "alternate": "AT", "type": "Insertion", "quality": null, "id": "." }
}
```

## Performance

- **Time Complexity:** O(L) where L is the allele length
- **Space Complexity:** O(1)

## See Also

- [normalize_variant](normalize_variant.md) — left-align and trim a variant
- [classify_variant](classify_variant.md) — classify REF/ALT into a variant type
