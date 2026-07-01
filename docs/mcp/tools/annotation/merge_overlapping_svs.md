# merge_overlapping_svs

Merge overlapping structural variants of the same type.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `merge_overlapping_svs` |
| **Method ID** | `StructuralVariantAnalyzer.MergeOverlappingSVs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Sorts SVs by chromosome and start, then merges consecutive same-chromosome, same-type variants whose
overlap fraction (overlap length ÷ smaller SV length) is at least `overlapFraction`. A merged SV spans the
union of the two, keeps the first SV's id, sums their supporting reads, and takes the maximum quality.
Because merging is over consecutive sorted entries, only adjacent overlapping SVs of the same type combine.

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L1055](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L1055)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variants` | array | Yes | Structural variants to merge |
| `overlapFraction` | number | No | Minimum overlap fraction 0..1 (default 0.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `merged` | array | Merged structural variants |

## Errors

| Code | Message |
|------|---------|
| 1001 | variants cannot be null |

## Examples

### Example 1: Two overlapping deletions merge

del1 (1000–2000) and del2 (1500–2500) overlap by 0.5 of the min length and merge into 1000–2500 with
summed support (8) and max quality (40).

### Example 2: Different types not merged

A Deletion and a Duplication at the same span are kept separate:

**Response:**
```json
{ "merged": [ { "type": "Deletion" }, { "type": "Duplication" } ] }
```

## Performance

- **Time Complexity:** O(n log n) for n SVs
- **Space Complexity:** O(n)

## See Also

- [cluster_discordant_pairs](cluster_discordant_pairs.md) / [identify_cnvs](identify_cnvs.md) — SV sources
- [filter_svs](filter_svs.md) — filter SV calls
