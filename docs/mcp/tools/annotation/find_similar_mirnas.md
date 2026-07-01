# find_similar_mirnas

Find miRNAs with seed regions similar to a query.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_similar_mirnas` |
| **Method ID** | `MiRnaAnalyzer.FindSimilarMiRnas` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the database miRNAs whose seed differs from the query seed by at most `maxMismatches` positions
(Hamming distance over the overlapping seed length). The query is excluded from its own results by name.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2679](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2679)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | object | Yes | Query miRNA record |
| `database` | array | Yes | miRNA records to search |
| `maxMismatches` | integer | No | Maximum allowed seed mismatches (default 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | array | Database miRNAs within the mismatch budget |

## Errors

| Code | Message |
|------|---------|
| 1001 | query cannot be null |
| 1001 | database cannot be null |

## Examples

### Example 1: Default 1-mismatch budget

Query seed `GAGGUAG`; a 0-mismatch and a 1-mismatch database entry are returned, a 2-mismatch entry is not:

**Response:**
```json
{ "matches": [ { "name": "same-seed" }, { "name": "one-off" } ] }
```

### Example 2: Exact seed only

With `maxMismatches = 0`, only exact seed matches are returned:

**Response:**
```json
{ "matches": [ { "name": "same-seed" } ] }
```

## Performance

- **Time Complexity:** O(n · L) for n database miRNAs and seed length L
- **Space Complexity:** O(k) for k matches

## See Also

- [group_by_seed_family](group_by_seed_family.md) — exact seed grouping
- [compare_seed_regions](compare_seed_regions.md) — pairwise seed comparison
