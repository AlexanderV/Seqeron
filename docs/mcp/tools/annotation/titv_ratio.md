# titv_ratio

Compute transition/transversion (Ti/Tv) ratio from a list of variants.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `titv_ratio` |
| **Method ID** | `VariantCaller.CalculateTiTvRatio` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Counts transitions (purine↔purine or pyrimidine↔pyrimidine SNPs) and transversions
(purine↔pyrimidine SNPs) across the provided variants — only `SNP`-typed entries are counted — and
returns their ratio `transitions / transversions`. When there are no transversions the ratio is `0`.
The Ti/Tv ratio is a standard sequencing-quality metric (whole-genome ≈2.0–2.1).

## Core Documentation Reference

- Source: [VariantCaller.cs#L203](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L203)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variants` | array | Yes | Variants (only SNP-typed entries are counted) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `ratio` | number | `transitions / transversions`; 0 when there are no transversions |

## Errors

| Code | Message |
|------|---------|
| 1001 | Variants cannot be null |

## Examples

### Example 1: One transition, one transversion

A>G (transition) and A>C (transversion) → 1/1 = 1.0.

**Response:**
```json
{ "ratio": 1.0 }
```

### Example 2: No transversions

A single A>G transition → ratio 0.0 (division guarded).

**Response:**
```json
{ "ratio": 0.0 }
```

## Performance

- **Time Complexity:** O(n)
- **Space Complexity:** O(1)

## See Also

- [classify_mutation](classify_mutation.md) - Classify a single SNP
- [variant_statistics](variant_statistics.md) - Summary statistics including Ti/Tv
