# methylation_from_bisulfite

Compute per-CpG methylation levels from bisulfite sequencing reads.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `methylation_from_bisulfite` |
| **Method ID** | `EpigeneticsAnalyzer.CalculateMethylationFromBisulfite` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calls methylation levels at every reference CpG from aligned bisulfite reads, following the Bismark
convention (Krueger & Andrews 2011): at a reference CpG cytosine a read base of `C` is a methylated
(protected) call and a read base of `T` is an unmethylated (converted) call; any other base is ignored.
Each covered CpG reports `methylationLevel = methylated / (methylated + unmethylated)` and `coverage`
= the number of valid C/T calls. CpG sites with zero coverage are omitted. Read bases beyond the
reference (or at its last position) are ignored.

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L434](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L434)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `referenceSequence` | string | Yes | Reference DNA sequence (min length: 1) |
| `bisulfiteReads` | array | Yes | Reads as `{ readSequence, startPosition }`; at least one |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].position` | integer | 0-based CpG position |
| `sites[].type` | string | Always `CpG` |
| `sites[].context` | string | Up-to-3-base context at the CpG |
| `sites[].methylationLevel` | number | methylated / (methylated + unmethylated) |
| `sites[].coverage` | integer | Number of valid C/T calls |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference sequence cannot be null or empty |
| 1001 | Bisulfite reads cannot be null or empty |

## Examples

### Example 1: Two CpG sites with mixed methylation

Reference `ACGTACGT` has CpGs at positions 1 and 5. Two reads aligned at 0 — one fully methylated, one
converted at site 5:

**Response:**
```json
{
  "sites": [
    { "position": 1, "type": "CpG", "context": "CGT", "methylationLevel": 1.0, "coverage": 2 },
    { "position": 5, "type": "CpG", "context": "CGT", "methylationLevel": 0.5, "coverage": 2 }
  ]
}
```

### Example 2: Uncovered CpG excluded

A single short read covering only site 1 leaves site 5 uncovered:

**Response:**
```json
{ "sites": [ { "position": 1, "type": "CpG", "context": "CGT", "methylationLevel": 1.0, "coverage": 1 } ] }
```

## Performance

- **Time Complexity:** O(R·L) for R reads of length L
- **Space Complexity:** O(S) for S CpG sites

## See Also

- [simulate_bisulfite_conversion](simulate_bisulfite_conversion.md) — forward simulation of conversion
- [methylation_profile](methylation_profile.md) — aggregate sites into a profile
