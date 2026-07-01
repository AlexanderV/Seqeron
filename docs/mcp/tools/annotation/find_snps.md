# find_snps

Detect only SNPs between two DNA sequences using global alignment.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_snps` |
| **Method ID** | `VariantCaller.FindSnps` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Globally aligns the query against the reference (via `call_variants`) and returns **only** the
single-nucleotide substitutions, filtering out insertions and deletions. Each SNP carries its
0-based reference `position`, the `referenceAllele`/`alternateAllele`, the `type` (always `SNP`),
and the corresponding `queryPosition`. Identical sequences yield no SNPs.

## Core Documentation Reference

- Source: [VariantCaller.cs#L117](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L117)

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

`ATGCATGC` vs `ATGAATGC` → one `C>A` SNP at position 3.

**Response:**
```json
{ "variants": [ { "position": 3, "referenceAllele": "C", "alternateAllele": "A", "type": "SNP", "queryPosition": 3 } ] }
```

### Example 2: Identical sequences

`ATGCATGC` vs `ATGCATGC` → no SNPs.

## Performance

- **Time Complexity:** O(n·m) global alignment
- **Space Complexity:** O(n·m)

## See Also

- [find_snps_direct](find_snps_direct.md) - SNPs by direct positional comparison (no alignment)
- [call_variants](call_variants.md) - All variant types
- [titv_ratio](titv_ratio.md) - Transition/transversion ratio
