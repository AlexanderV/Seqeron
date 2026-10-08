# find_snps_direct

Detect SNPs by direct positional comparison without alignment. Inputs must be pre-aligned and of equal length (unequal lengths are rejected); bases compare case-insensitively and gap columns are not SNPs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_snps_direct` |
| **Method ID** | `VariantCaller.FindSnpsDirect` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Compares the pre-aligned reference and query column by column (no alignment is run), reporting a
`SNP` for every column whose two bases differ (case-insensitive, VCF v4.3). It delegates to
`VariantCaller.CallVariantsFromAlignment` and keeps the SNP columns: gap (`-`) columns are indels and
are not reported, and positions are ungapped coordinates. For gap-free inputs `queryPosition` equals
`position` and the SNP count equals the Hamming distance. The inputs must have equal length (Hamming
distance is undefined otherwise); unequal lengths are rejected. Use this when the sequences are already aligned/registered; use
[`find_snps`](find_snps.md) when they may contain indels.

## Core Documentation Reference

- Source: [VariantCaller.cs#L125](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L125)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference DNA sequence |
| `query` | string | Yes | Query DNA sequence |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants` | array | `{ position, referenceAllele, alternateAllele, type, queryPosition }` per SNP |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference cannot be null or empty |
| 1001 | Query cannot be null or empty |
| 1001 | Aligned sequences must have the same length. |

## Examples

### Example 1: Single substitution

`ATGC` vs `ATTC` → one `G>T` SNP at position 2.

**Response:**
```json
{ "variants": [ { "position": 2, "referenceAllele": "G", "alternateAllele": "T", "type": "SNP", "queryPosition": 2 } ] }
```

### Example 2: Unequal lengths are rejected

`ATGCAA` vs `ATTC` → error "Aligned sequences must have the same length." (positional comparison is
undefined for unequal lengths).

### Example 3: Case-insensitive

`acgt` vs `ACGT` → no variants (soft-masked bases are matches).

## Performance

- **Time Complexity:** O(n)
- **Space Complexity:** O(k)

## See Also

- [find_snps](find_snps.md) - Alignment-based SNP detection (indel-tolerant)
- [call_variants](call_variants.md) - All variant types
