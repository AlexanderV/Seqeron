# normalize_variant

Normalize a variant by trimming common prefixes/suffixes, then classify its type.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `normalize_variant` |
| **Method ID** | `VariantAnnotator.NormalizeVariant` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Normalizes a variant to its minimal representation: the alleles are uppercased, then their common
**suffix** is trimmed, then their common **prefix** is trimmed (advancing the position by one per
trimmed prefix base). The resulting minimal ref/alt pair is classified via
[`classify_variant`](classify_variant.md), and the normalized `chromosome`, `position`, alleles, and
`type` are returned. The optional `referenceSequence` is reserved for left-alignment and currently
unused.

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L228](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L228)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chromosome` | string | Yes | Chromosome / contig identifier |
| `position` | integer | Yes | 1-based genomic position of the reference allele |
| `reference` | string | Yes | Reference allele |
| `alternate` | string | Yes | Alternate allele |
| `referenceSequence` | string | No | Optional reference sequence (reserved for left-alignment) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variant` | object | Normalized `{ chromosome, position, reference, alternate, type, ... }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Chromosome cannot be null or empty |
| 1001 | Reference allele cannot be null or empty |
| 1001 | Alternate allele cannot be null or empty |

## Examples

### Example 1: Common suffix trimmed

`ACG`/`ATG` at position 100 → the shared `G` suffix is trimmed, leaving `C`/`T` at position 101
(SNV).

**Response:**
```json
{ "variant": { "chromosome": "chr1", "position": 101, "reference": "C", "alternate": "T", "type": "SNV" } }
```

### Example 2: Already-minimal SNV

`A`/`T` at position 100 → unchanged.

## Performance

- **Time Complexity:** O(allele length)
- **Space Complexity:** O(allele length)

## See Also

- [classify_variant](classify_variant.md) - Classify a (ref, alt) pair
- [parse_vcf_variant](parse_vcf_variant.md) - Build a variant from VCF fields
