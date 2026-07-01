# classify_variant

Classify a (ref, alt) allele pair as SNV / Insertion / Deletion / MNV / Indel / Complex.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `classify_variant` |
| **Method ID** | `VariantAnnotator.ClassifyVariant` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Classifies a reference/alternate allele pair by length comparison: equal single bases → `SNV`; equal
multi-base → `MNV`; a longer alternate that starts with the reference → `Insertion`; a longer
reference whose single-base alternate is a prefix → `Deletion`; any other length mismatch → `Indel`;
an empty allele → `Complex`.

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L196](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L196)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference allele |
| `alternate` | string | Yes | Alternate allele |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variantType` | string | `SNV` / `Insertion` / `Deletion` / `MNV` / `Indel` / `Complex` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference allele cannot be null or empty |
| 1001 | Alternate allele cannot be null or empty |

## Examples

### Example 1: SNV

`A` → `G`.

**Response:**
```json
{ "variantType": "SNV" }
```

### Example 2: Insertion

`A` → `ACGT` (alternate longer, prefixed by reference).

**Response:**
```json
{ "variantType": "Insertion" }
```

## Performance

- **Time Complexity:** O(allele length)
- **Space Complexity:** O(1)

## See Also

- [normalize_variant](normalize_variant.md) - Normalize and classify a variant
- [classify_mutation](classify_mutation.md) - Transition/transversion class of a SNP
