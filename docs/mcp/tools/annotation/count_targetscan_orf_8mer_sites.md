# count_targetscan_orf_8mer_sites

ORF8m: 8mer sites of a seed region in an ORF as counted by targetscan_count_8mers.pl.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `count_targetscan_orf_8mer_sites` |
| **Method ID** | `MiRnaAnalyzer.CountTargetScanOrf8merSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Gap characters ('-', '.') are removed, T is read as U, and reverse_complement(seed region) + 'A' is counted case-insensitively and non-overlapping on the given strand only. Feed `count` to `score_target_site_context_plus_plus` (`orf8mCount`).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2049](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2049)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `orfSequence` | string | Yes | ORF sequence, RNA or DNA, any case, gaps allowed |
| `seedRegion` | string | Yes | 7-nt seed region |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `count` | integer | Number of 8mer sites |

## Errors

| Code | Message |
|------|---------|
| 1001 | orfSequence cannot be null |
| 1001 | seedRegion must be 7 nt (miRNA nt 2-8) |

## Examples

### Example 1

Overlap-free counting: seed AAAAAAU, site AUUUUUUA, 3 sites (an overlapping scan would give 4).

**Request:**
```json
{"orfSequence": "GGAUUUUUUAUUUUUUAccauuuuuuaATTTTTTAG", "seedRegion": "AAAAAAU"}
```

**Response (abridged):**
```json
{"count": 3}
```

## Performance

- **Time Complexity:** linear in the input size (alignment columns × species for the UTR bin)
- **Space Complexity:** O(n)

## See Also

- [score_targetscan_conservation](score_targetscan_conservation.md)
- [score_target_site_context_plus_plus](score_target_site_context_plus_plus.md)
