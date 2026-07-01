# find_insertions

Detect only insertions between two DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_insertions` |
| **Method ID** | `VariantCaller.FindInsertions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Globally aligns the query against the reference (via `call_variants`) and returns **only** the
insertions — positions where the query carries extra bases relative to the reference. Each insertion
uses the gap sentinel `-` as its `referenceAllele`, the inserted base(s) as its `alternateAllele`, a
0-based reference `position`, and `type` `Insertion`. Substitutions and deletions are filtered out.

## Core Documentation Reference

- Source: [VariantCaller.cs#L153](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L153)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference DNA sequence |
| `query` | string | Yes | Query DNA sequence |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants` | array | `{ position, referenceAllele, alternateAllele, type, queryPosition }` per insertion |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference cannot be null or empty |
| 1001 | Query cannot be null or empty |

## Examples

### Example 1: Single inserted base

`ATGCAT` vs `ATGTCAT` → one insertion of `T` at position 3 (`referenceAllele` = `-`).

**Response:**
```json
{ "variants": [ { "position": 3, "referenceAllele": "-", "alternateAllele": "T", "type": "Insertion", "queryPosition": 3 } ] }
```

### Example 2: Substitution-only input

`ATGCATGC` vs `ATGAATGC` → no insertions.

## Performance

- **Time Complexity:** O(n·m) global alignment
- **Space Complexity:** O(n·m)

## See Also

- [find_deletions](find_deletions.md) - Deletions only
- [find_indels](find_indels.md) - Insertions and deletions
- [call_variants](call_variants.md) - All variant types
