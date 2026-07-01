# call_variants_from_alignment

Detect variants from already-aligned reference and query strings.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `call_variants_from_alignment` |
| **Method ID** | `VariantCaller.CallVariantsFromAlignment` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Walks two equal-length aligned strings column by column (gaps encoded as `-`) and emits a variant
for each difference without re-aligning:

- a reference-side gap (`-` in the reference, base in the query) → **Insertion**,
- a query-side gap (base in the reference, `-` in the query) → **Deletion**,
- a mismatched column (both bases, different) → **SNP**.

`position` is the 0-based reference coordinate and `queryPosition` the 0-based query coordinate at
the variant. The two aligned strings must be the same length, otherwise an error is raised.

## Core Documentation Reference

- Source: [VariantCaller.cs#L41](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L41)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedReference` | string | Yes | Aligned reference sequence (gaps as `-`) |
| `alignedQuery` | string | Yes | Aligned query sequence (gaps as `-`), same length |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants` | array | `{ position, referenceAllele, alternateAllele, type, queryPosition }` per called difference |

## Errors

| Code | Message |
|------|---------|
| 1001 | Aligned reference cannot be null or empty |
| 1001 | Aligned query cannot be null or empty |
| 2001 | Aligned sequences must have the same length. |

## Examples

### Example 1: Insertion (reference gap)

`AT-GC` / `ATXGC` → an insertion of `X` at reference position 2.

**Response:**
```json
{ "variants": [ { "position": 2, "referenceAllele": "-", "alternateAllele": "X", "type": "Insertion", "queryPosition": 2 } ] }
```

### Example 2: Deletion (query gap)

`ATGC` / `AT-C` → a deletion of `G` at reference position 2.

**Response:**
```json
{ "variants": [ { "position": 2, "referenceAllele": "G", "alternateAllele": "-", "type": "Deletion", "queryPosition": 2 } ] }
```

## Performance

- **Time Complexity:** O(n) over the alignment columns
- **Space Complexity:** O(v) in the number of variants

## See Also

- [call_variants](call_variants.md) - Align then call variants
- [find_indels](find_indels.md) - Insertions and deletions only
