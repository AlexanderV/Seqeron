# get_targetscan_parameters

Bundled TargetScan 7.0 PCT family / PCT sigmoid row / TA / SPS lookup keyed by 7-nt seed region.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `get_targetscan_parameters` |
| **Method ID** | `MiRnaAnalyzer.TryGetTargetScanPctParameters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

One lookup over `IsTargetScanPctSeedFamily`, `TryGetTargetScanPctParameters` and `TryGetTargetScanTaSps` (T read as U, any case). `pctParameters` needs `siteType` Seed8mer/Seed7merM8/Seed7merA1; `sps` needs a siteType with an SPS column (8mer, 7mer-m8, 7mer-A1, 6mer); `ta` is per seed. `hasTaSps` is false unless the seed is 7 nt over ACGU.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2294](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2294)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seedRegion` | string | Yes | 7-nt seed region (miRNA nt 2-8) |
| `siteType` | string | No | Selects the PCT row and the SPS column |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `seedRegion` | string | Echo |
| `isPctSeedFamily` | boolean | One of the 111 PCT seed families |
| `pctParameters` | object | null | b0-b3 |
| `hasTaSps` | boolean | TA/SPS row exists |
| `ta` | number | null | log10 site abundance |
| `sps` | number | null | Seed-pairing stability (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | seedRegion cannot be null or empty |
| 1001 | Unknown siteType |

## Examples

### Example 1

let-7 seed, 7mer-m8: TA 3.393, SPS -9.25 (TA_SPS_by_seed_region.txt) and the let-7 PCT row.

**Request:**
```json
{"seedRegion": "GAGGUAG", "siteType": "Seed7merM8"}
```

**Response (abridged):**
```json
{"seedRegion": "GAGGUAG", "isPctSeedFamily": true, "pctParameters": {"b0": -0.221367319203409, "b1": 1.24081845475817, "b2": 1.17370692612852, "b3": 3.41806979609985}, "hasTaSps": true, "ta": 3.393, "sps": -9.25}
```

## Performance

- **Time Complexity:** linear in the input size (alignment columns × species for the UTR bin)
- **Space Complexity:** O(n)

## See Also

- [score_targetscan_conservation](score_targetscan_conservation.md)
- [score_target_site_context_plus_plus](score_target_site_context_plus_plus.md)
