# identify_cnvs

Convert copy-number segments into deletion / duplication SVs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `identify_cnvs` |
| **Method ID** | `StructuralVariantAnalyzer.IdentifyCNVs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Turns copy-number segments whose copy number differs from `normalCopyNumber` and whose length is at least
`minLength` into structural variants: **Deletion** when the copy number is below baseline, **Duplication**
when above. Each CNV carries an `id` (`CNV1`, `CNV2`, …), its span/length, a quality of `|logRatio| · 50`,
and the segment probe count as supporting reads. Baseline or short segments are skipped.

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L794](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L794)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `segments` | array | Yes | Copy-number segments |
| `normalCopyNumber` | integer | No | Diploid baseline copy number (default 2) |
| `minLength` | integer | No | Minimum segment length to emit (default 10000) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants[]` | object | CNVs (id, chromosome, span, type, length, quality, supportingReads) |

## Errors

| Code | Message |
|------|---------|
| 1001 | segments cannot be null |

## Examples

### Example 1: Deletion and duplication

A CN-1 (19 kb) and a CN-4 (30 kb) segment become a Deletion and a Duplication.

### Example 2: Baseline only

A CN-2 baseline segment yields nothing:

**Response:**
```json
{ "variants": [] }
```

## Performance

- **Time Complexity:** O(n) for n segments
- **Space Complexity:** O(k) for k CNVs

## See Also

- [segment_copy_number](segment_copy_number.md) — produce the segments
- [annotate_svs](annotate_svs.md) — annotate CNVs with genes
