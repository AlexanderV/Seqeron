# find_indels

Detect insertions and deletions (indels) between two DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_indels` |
| **Method ID** | `VariantCaller.FindIndels` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Globally aligns the query against the reference (via `call_variants`) and returns the insertions and
deletions together, filtering out substitutions. Insertions use the gap sentinel `-` as
`referenceAllele`; deletions use `-` as `alternateAllele`. Each indel carries a 0-based reference
`position` and a `type` of `Insertion` or `Deletion`.

## Core Documentation Reference

- Source: [VariantCaller.cs#L169](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L169)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference DNA sequence |
| `query` | string | Yes | Query DNA sequence |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants` | array | `{ position, referenceAllele, alternateAllele, type, queryPosition }` per indel |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference cannot be null or empty |
| 1001 | Query cannot be null or empty |

## Examples

### Example 1: Insertion

`ATGCAT` vs `ATGTCAT` → one `Insertion` of `T` at position 3.

**Response:**
```json
{ "variants": [ { "position": 3, "referenceAllele": "-", "alternateAllele": "T", "type": "Insertion", "queryPosition": 3 } ] }
```

### Example 2: Deletion

`ATGTCAT` vs `ATGCAT` → one `Deletion` of `T` at position 3.

**Response:**
```json
{ "variants": [ { "position": 3, "referenceAllele": "T", "alternateAllele": "-", "type": "Deletion", "queryPosition": 3 } ] }
```

## Performance

- **Time Complexity:** O(n·m) global alignment
- **Space Complexity:** O(n·m)

## See Also

- [find_insertions](find_insertions.md) - Insertions only
- [find_deletions](find_deletions.md) - Deletions only
- [call_variants](call_variants.md) - All variant types
