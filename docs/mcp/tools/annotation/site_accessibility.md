# site_accessibility

Estimate miRNA target-site accessibility from local secondary structure.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `site_accessibility` |
| **Method ID** | `MiRnaAnalyzer.CalculateSiteAccessibility` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Estimates how accessible a target site is to a miRNA by measuring local secondary-structure density in a
window spanning ±50 nt around the site. The structure score counts Watson-Crick (non-wobble) pairs at
least 4 nt apart; accessibility is `max(0, 1 − density·10)` where `density = structureScore / maxPairs`
and `maxPairs = W·(W−4)/2` for window width `W`. Higher values indicate a more open, targetable site.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2630](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2630)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `mRnaSequence` | string | Yes | mRNA nucleotide sequence |
| `siteStart` | integer | Yes | 0-based inclusive site start |
| `siteEnd` | integer | Yes | 0-based inclusive site end |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `accessibility` | number | `max(0, 1 − density·10)`; higher = more accessible |

## Errors

| Code | Message |
|------|---------|
| 1001 | mRNA sequence cannot be null or empty |
| 1001 | Site indices are out of range (ArgumentOutOfRangeException) |

## Examples

### Example 1: Windowed structure density

For `GAAAAUAAAC` and site (2, 7), the window covers the whole sequence: two non-wobble pairs (G0-C9,
A1-U5), `maxPairs = 30`, so accessibility = `1 − (2/30)·10`:

**Response:**
```json
{ "accessibility": 0.33333333333333337 }
```

### Example 2: Site starting at position 0

`siteStart = 0` is valid and yields the same window/value:

**Response:**
```json
{ "accessibility": 0.33333333333333337 }
```

## Performance

- **Time Complexity:** O(W²) for window width W (≤ ~100 nt)
- **Space Complexity:** O(1)

## See Also

- [analyze_target_context](analyze_target_context.md) — AU/context score around a site
- [find_mirna_target_sites](find_mirna_target_sites.md) — canonical target-site search
