# filter_svs

Filter structural variants by quality, support, and length.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `filter_svs` |
| **Method ID** | `StructuralVariantAnalyzer.FilterSVs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Keeps structural variants whose quality is at least `minQuality`, whose supporting-read count is at least
`minSupport`, and whose length is within `[minLength, maxLength]`. All four conditions must hold.

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L1114](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L1114)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variants` | array | Yes | Structural variants to filter |
| `minQuality` | number | No | Minimum quality (default 20) |
| `minSupport` | integer | No | Minimum supporting reads (default 2) |
| `minLength` | integer | No | Minimum SV length (default 50) |
| `maxLength` | integer | No | Maximum SV length (default 100000000) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants` | array | Structural variants passing all filters |

## Errors

| Code | Message |
|------|---------|
| 1001 | variants cannot be null |

## Examples

### Example 1: Default thresholds

Only the SV with sufficient quality, support and length passes.

### Example 2: Relaxed thresholds

Lowering all minimums lets short/low-quality SVs through.

## Performance

- **Time Complexity:** O(n) for n SVs
- **Space Complexity:** O(k) for k passing SVs

## See Also

- [merge_overlapping_svs](merge_overlapping_svs.md) — merge SV calls
- [genotype_sv](genotype_sv.md) — genotype a filtered SV
