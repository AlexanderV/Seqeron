# compare_seed_regions

Compare the seed regions of two miRNAs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `compare_seed_regions` |
| **Method ID** | `MiRnaAnalyzer.CompareSeedRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Compares the seed strings of two miRNA records position-by-position, returning the number of matching
positions, the mismatch count (Hamming distance over the overlapping length plus any length difference),
and `isSameFamily` — `true` only when the two seed strings are exactly equal. miRNAs sharing a seed are
predicted to regulate an overlapping set of targets (a seed family).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L124](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L124)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `miRna1` | object | Yes | First miRNA record (with seed metadata) |
| `miRna2` | object | Yes | Second miRNA record |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | integer | Positions where the seeds agree |
| `mismatches` | integer | Hamming distance plus length difference |
| `isSameFamily` | boolean | true when the seed strings are identical |

## Errors

| Code | Message |
|------|---------|
| 1001 | miRna1 cannot be null |
| 1001 | miRna2 cannot be null |
| 1001 | Both miRNAs must have non-empty sequences |

## Examples

### Example 1: Same seed family

Two miRNAs with the identical seed `GAGGUAG`:

**Response:**
```json
{ "matches": 7, "mismatches": 0, "isSameFamily": true }
```

### Example 2: One-nucleotide difference

Seeds `GAGGUAG` vs `GAGGUCG`:

**Response:**
```json
{ "matches": 6, "mismatches": 1, "isSameFamily": false }
```

## Performance

- **Time Complexity:** O(L) for seed length L
- **Space Complexity:** O(1)

## See Also

- [create_mirna](create_mirna.md) — build miRNA records with seed metadata
- [group_by_seed_family](group_by_seed_family.md) — cluster miRNAs by shared seed
