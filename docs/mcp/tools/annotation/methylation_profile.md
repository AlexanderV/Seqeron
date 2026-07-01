# methylation_profile

Aggregate methylation sites into a global / CpG / CHG / CHH profile.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `methylation_profile` |
| **Method ID** | `EpigeneticsAnalyzer.GenerateMethylationProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Summarises a set of measured methylation sites into per-context **weighted methylation levels**
(Schultz et al. 2012): for each context the level is `Σ(level·coverage) / Σ(coverage)`, i.e. the fraction
of methylated reads over total reads (falling back to an unweighted mean when total coverage is zero).
It also reports the total number of CpG sites, a descriptive count of CpG sites with level ≥ 0.5, and a
per-position methylation list sorted by position.

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L507](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L507)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sites` | array | Yes | Methylation sites (`position`, `type`, `context`, `methylationLevel`, `coverage`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `globalMethylation` | number | Weighted mean over all sites |
| `cpGMethylation` | number | Weighted mean over CpG sites |
| `cHGMethylation` | number | Weighted mean over CHG sites |
| `cHHMethylation` | number | Weighted mean over CHH sites |
| `totalCpGSites` | integer | Number of CpG sites |
| `methylatedCpGSites` | integer | CpG sites with level ≥ 0.5 |
| `methylationByPosition` | array | Per-position `{ position, level }`, sorted by position |

## Errors

| Code | Message |
|------|---------|
| 1001 | sites cannot be null |
| 1002 | Unknown methylation type '&lt;value&gt;'. |

## Examples

### Example 1: Mixed CpG/CHG profile

Two CpG sites (levels 1.0 and 0.0, coverage 10 each) and one CHG site (level 0.5, coverage 4):

**Response:**
```json
{
  "globalMethylation": 0.5,
  "cpGMethylation": 0.5,
  "cHGMethylation": 0.5,
  "cHHMethylation": 0.0,
  "totalCpGSites": 2,
  "methylatedCpGSites": 1,
  "methylationByPosition": [
    { "position": 10, "level": 1.0 },
    { "position": 20, "level": 0.0 },
    { "position": 30, "level": 0.5 }
  ]
}
```

### Example 2: Empty input

```json
{ "tool": "methylation_profile", "arguments": { "sites": [] } }
```

**Response:**
```json
{ "globalMethylation": 0, "cpGMethylation": 0, "cHGMethylation": 0, "cHHMethylation": 0, "totalCpGSites": 0, "methylatedCpGSites": 0, "methylationByPosition": [] }
```

## Performance

- **Time Complexity:** O(n log n) (sort by position) for n sites
- **Space Complexity:** O(n)

## See Also

- [methylation_from_bisulfite](methylation_from_bisulfite.md) — produce measured sites from reads
- [find_dmrs](find_dmrs.md) — compare two methylation profiles for DMRs
