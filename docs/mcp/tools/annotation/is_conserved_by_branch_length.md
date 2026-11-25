# is_conserved_by_branch_length

TargetScan 'Conserved' call: Bls ≥ the site-type threshold (human or mouse tree).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `is_conserved_by_branch_length` |
| **Method ID** | `MiRnaAnalyzer.IsConservedByBranchLength` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

TargetScanHuman7 thresholds: 8mer 1.8, 7mer-m8 2.8, 7mer-A1 3.6; TargetScanMouse7: 0.6, 1.8, 2.5. A 6mer or offset 6mer is never conserved (`threshold` null).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2450](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2450)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `branchLengthScore` | number | Yes | Bls ≥ 0 |
| `siteType` | string | Yes | Seed8mer, Seed7merM8, Seed7merA1, Seed6mer, Offset6mer |
| `thresholds` | string | No | TargetScanHuman7 (default) or TargetScanMouse7 |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `conserved` | boolean | Conserved call |
| `threshold` | number | null | Threshold applied |

## Errors

| Code | Message |
|------|---------|
| 1001 | Unknown siteType |
| 1001 | Unknown thresholds |
| 1001 | Branch-length score must be a finite value ≥ 0 |

## Examples

### Example 1

An 8mer with Bls 1.8 is conserved on the human thresholds.

**Request:**
```json
{"branchLengthScore": 1.8, "siteType": "Seed8mer"}
```

**Response (abridged):**
```json
{"conserved": true, "threshold": 1.8}
```

## Performance

- **Time Complexity:** linear in the input size (alignment columns × species for the UTR bin)
- **Space Complexity:** O(n)

## See Also

- [score_targetscan_conservation](score_targetscan_conservation.md)
- [score_target_site_context_plus_plus](score_target_site_context_plus_plus.md)
