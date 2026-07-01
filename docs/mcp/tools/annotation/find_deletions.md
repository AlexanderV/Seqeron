# find_deletions

Detect only deletions between two DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_deletions` |
| **Method ID** | `VariantCaller.FindDeletions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Globally aligns the query against the reference (via `call_variants`) and returns **only** the
deletions — positions where the reference carries bases absent from the query. Each deletion uses the
deleted base(s) as its `referenceAllele`, the gap sentinel `-` as its `alternateAllele`, a 0-based
reference `position`, and `type` `Deletion`. Substitutions and insertions are filtered out.

## Core Documentation Reference

- Source: [VariantCaller.cs#L161](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L161)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference DNA sequence |
| `query` | string | Yes | Query DNA sequence |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants` | array | `{ position, referenceAllele, alternateAllele, type, queryPosition }` per deletion |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference cannot be null or empty |
| 1001 | Query cannot be null or empty |

## Examples

### Example 1: Single deleted base

`ATGTCAT` vs `ATGCAT` → one deletion of `T` at position 3 (`alternateAllele` = `-`).

**Response:**
```json
{ "variants": [ { "position": 3, "referenceAllele": "T", "alternateAllele": "-", "type": "Deletion", "queryPosition": 3 } ] }
```

### Example 2: Substitution-only input

`ATGCATGC` vs `ATGAATGC` → no deletions.

## Performance

- **Time Complexity:** O(n·m) global alignment
- **Space Complexity:** O(n·m)

## See Also

- [find_insertions](find_insertions.md) - Insertions only
- [find_indels](find_indels.md) - Insertions and deletions
- [call_variants](call_variants.md) - All variant types
