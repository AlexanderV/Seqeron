# variant_statistics

Compute summary variant statistics between reference and query (totals, Ti/Tv, density).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `variant_statistics` |
| **Method ID** | `VariantCaller.CalculateStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calls all variants between the reference and query (global alignment) and summarizes them: total
variant count, per-type counts (`snps`, `insertions`, `deletions`), the Ti/Tv ratio, the variant
density (variants per 1000 reference bases), and the two sequence lengths. Variant density is
`variants / referenceLength × 1000` (0 for an empty reference).

## Core Documentation Reference

- Source: [VariantCaller.cs#L229](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L229)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference DNA sequence |
| `query` | string | Yes | Query DNA sequence |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalVariants` | integer | Total called variants |
| `snps` | integer | SNP count |
| `insertions` | integer | Insertion count |
| `deletions` | integer | Deletion count |
| `tiTvRatio` | number | Transition/transversion ratio |
| `variantDensity` | number | Variants per 1000 reference bases |
| `referenceLength` | integer | Reference length |
| `queryLength` | integer | Query length |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference cannot be null or empty |
| 1001 | Query cannot be null or empty |

## Examples

### Example 1: Single SNP

`ATGC` vs `ATTC` → one SNP (G>T, a transversion), density 250 per kb.

**Response:**
```json
{ "totalVariants": 1, "snps": 1, "insertions": 0, "deletions": 0, "tiTvRatio": 0.0, "variantDensity": 250.0, "referenceLength": 4, "queryLength": 4 }
```

### Example 2: Identical sequences

`ATGCATGC` vs `ATGCATGC` → zero variants.

## Performance

- **Time Complexity:** O(n·m) global alignment
- **Space Complexity:** O(n·m)

## See Also

- [call_variants](call_variants.md) - Raw variant list
- [titv_ratio](titv_ratio.md) - Ti/Tv ratio only
