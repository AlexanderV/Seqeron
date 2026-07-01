# group_by_seed_family

Group miRNAs by identical seed sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `group_by_seed_family` |
| **Method ID** | `MiRnaAnalyzer.GroupBySeedFamily` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Partitions a set of miRNAs into seed families by grouping on the exact seed sequence. miRNAs in the same
family share a seed and are predicted to regulate an overlapping set of targets.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2669](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2669)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `miRnas` | array | Yes | miRNA records to group |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `families[].seedFamily` | string | Shared seed sequence |
| `families[].members` | array | miRNA records in the family |

## Errors

| Code | Message |
|------|---------|
| 1001 | miRnas cannot be null |

## Examples

### Example 1: Two seed families

let-7a and let-7c share seed `GAGGUAG`; a third miRNA has `GAGGUCG`:

**Response:**
```json
{
  "families": [
    { "seedFamily": "GAGGUAG", "members": [ { "name": "let-7a" }, { "name": "let-7c" } ] },
    { "seedFamily": "GAGGUCG", "members": [ { "name": "other" } ] }
  ]
}
```

### Example 2: Empty input

```json
{ "families": [] }
```

## Performance

- **Time Complexity:** O(n) for n miRNAs
- **Space Complexity:** O(n)

## See Also

- [compare_seed_regions](compare_seed_regions.md) — pairwise seed comparison
- [find_similar_mirnas](find_similar_mirnas.md) — near-seed matches within a Hamming budget
