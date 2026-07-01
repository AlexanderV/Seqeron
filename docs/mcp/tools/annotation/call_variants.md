# call_variants

Detect SNPs and indels between two DNA sequences using global alignment.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `call_variants` |
| **Method ID** | `VariantCaller.CallVariants` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Globally aligns the query against the reference and reports every difference as a variant. Column
mismatches become `SNP`s; reference-side gaps become `Insertion`s (query has extra bases) and
query-side gaps become `Deletion`s. Each variant carries its 0-based reference `position`, the
`referenceAllele`/`alternateAllele`, its `type`, and the corresponding `queryPosition`. Identical
sequences yield no variants.

## Core Documentation Reference

- Source: [VariantCaller.cs#L27](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L27)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference DNA sequence |
| `query` | string | Yes | Query DNA sequence |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants` | array | `{ position, referenceAllele, alternateAllele, type, queryPosition }` per called difference |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference cannot be null or empty |
| 1001 | Query cannot be null or empty |

## Examples

### Example 1: Single SNP

`ATGC` vs `ATTC` → one `G>T` SNP at position 2.

**Response:**
```json
{ "variants": [ { "position": 2, "referenceAllele": "G", "alternateAllele": "T", "type": "SNP", "queryPosition": 2 } ] }
```

### Example 2: Identical sequences

**Response:**
```json
{ "variants": [] }
```

## Performance

- **Time Complexity:** O(n·m) global alignment
- **Space Complexity:** O(n·m)

## See Also

- [call_variants_from_alignment](call_variants_from_alignment.md) - Call variants from pre-aligned strings
- [find_snps](find_snps.md) - SNPs only
- [annotate_variants](annotate_variants.md) - Call variants and annotate effects
