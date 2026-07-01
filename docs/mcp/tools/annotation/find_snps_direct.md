# find_snps_direct

Detect SNPs by direct positional comparison without alignment.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_snps_direct` |
| **Method ID** | `VariantCaller.FindSnpsDirect` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Compares the reference and query base-by-base at the same index (no alignment), reporting a `SNP`
for every mismatched position over the common prefix (the shorter of the two lengths). Because it is
purely positional, `queryPosition` always equals `position`, and the SNP count equals the Hamming
distance for equal-length inputs. Use this when the sequences are already aligned/registered; use
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

## Examples

### Example 1: Single substitution

`ATGC` vs `ATTC` → one `G>T` SNP at position 2.

**Response:**
```json
{ "variants": [ { "position": 2, "referenceAllele": "G", "alternateAllele": "T", "type": "SNP", "queryPosition": 2 } ] }
```

### Example 2: Unequal lengths compare the common prefix

`ATGCAA` vs `ATTC` → only the 4-base common prefix is compared, yielding one SNP at position 2.

## Performance

- **Time Complexity:** O(min(n, m))
- **Space Complexity:** O(k)

## See Also

- [find_snps](find_snps.md) - Alignment-based SNP detection (indel-tolerant)
- [call_variants](call_variants.md) - All variant types
