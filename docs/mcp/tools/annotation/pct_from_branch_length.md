# pct_from_branch_length

Probability of conserved targeting (PCT) from a branch-length score: explicit sigmoid row or bundled TargetScan seed row.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `pct_from_branch_length` |
| **Method ID** | `MiRnaAnalyzer.PctFromBranchLength` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

PCT = b0 + b1 / (1 + exp(-b2·Bls + b3)), rounded to 4 decimals, negative → 0, Bls = 0 → 0 (targetscan_70_BL_PCT.pl). Mode 1: supply b0-b3 (optional `siteType`; Seed6mer → 0). Mode 2: supply `seedRegion` (7 nt, nt 2-8) and `siteType` (Seed8mer, Seed7merM8 or Seed7merA1) to use the bundled TargetScan 7.0 row; a seed outside the 111 PCT families returns `pct` null / `isNa` true (TargetScan 'NA'). Exactly one mode.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2229](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2229)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `branchLengthScore` | number | Yes | Bls ≥ 0 |
| `siteType` | string | No | Seed8mer, Seed7merM8, Seed7merA1, Seed6mer, Offset6mer |
| `seedRegion` | string | No | 7-nt seed region for the bundled lookup |
| `b0` | number | No | Sigmoid b0 |
| `b1` | number | No | Sigmoid b1 |
| `b2` | number | No | Sigmoid b2 |
| `b3` | number | No | Sigmoid b3 |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pct` | number | null | PCT; null = NA |
| `isNa` | boolean | No bundled row for the seed |
| `parameters` | object | null | Sigmoid row used |

## Errors

| Code | Message |
|------|---------|
| 1001 | Provide exactly one of seedRegion (+ siteType) or the sigmoid parameters b0..b3 |
| 1001 | b0, b1, b2 and b3 must be given together |
| 1001 | siteType is required with seedRegion |
| 1001 | TargetScan has PCT rows only for Seed8mer, Seed7merM8 and Seed7merA1 |
| 1001 | Branch-length score must be a finite value ≥ 0 |

## Examples

### Example 1

let-7 (GAGGUAG) 7mer-m8 with the kit Bls 1.6379 gives PCT 0.0058.

**Request:**
```json
{"branchLengthScore": 1.6379, "siteType": "Seed7merM8", "seedRegion": "GAGGUAG"}
```

**Response (abridged):**
```json
{"pct": 0.0058, "isNa": false, "parameters": {"b0": -0.221367319203409, "b1": 1.24081845475817, "b2": 1.17370692612852, "b3": 3.41806979609985}}
```

## Performance

- **Time Complexity:** linear in the input size (alignment columns × species for the UTR bin)
- **Space Complexity:** O(n)

## See Also

- [score_targetscan_conservation](score_targetscan_conservation.md)
- [score_target_site_context_plus_plus](score_target_site_context_plus_plus.md)
