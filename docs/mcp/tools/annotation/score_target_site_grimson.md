# score_target_site_grimson

TargetScan 5.0 / Grimson et al. (2007) context score of an 8mer, 7mer-m8 or 7mer-A1 site with its four contributions (log2 units).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `score_target_site_grimson` |
| **Method ID** | `MiRnaAnalyzer.ScoreTargetSiteGrimson` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Faithful port of the TargetScan 5.0 context score (Grimson et al. 2007, Mol Cell 27:91; `targetscan_50_context_scores.pl`): site type (−0.310 / −0.161 / −0.099 for 8mer / 7mer-m8 / 7mer-A1) + 3′ pairing + local AU + position contributions, in log2 fold-change units (negative = repressed). `mRnaSequence` is the 3′UTR (index 0 = first nt after the stop codon); `siteStart`/`siteEnd` are 0-based inclusive (an 8mer includes its A1). `contextScore` is null when the 1-based site start is < 15 (TargetScan `too_close`; the contributions are still returned). `scored = false` for `Seed6mer`, `Offset6mer`, `Supplementary` and `Centered` (no TargetScan 5.0 coefficients). `find_mirna_target_sites` reports `max(0, −contextScore)` as its `score`.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.GrimsonContext.cs#L75](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.GrimsonContext.cs#L75)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `mRnaSequence` | string | Yes | 3′UTR / mRNA nucleotide sequence (T→U applied internally) |
| `miRnaSequence` | string | Yes | miRNA nucleotide sequence |
| `siteStart` | integer | Yes | 0-based inclusive site start |
| `siteEnd` | integer | Yes | 0-based inclusive site end |
| `siteType` | string | Yes | `Seed8mer`, `Seed7merM8` or `Seed7merA1` (others are not scored) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `scored` | boolean | False for site types without coefficients |
| `score.type` | string | Site type scored |
| `score.siteTypeContribution` | number | −0.310 / −0.161 / −0.099 |
| `score.threePrimePairingContribution` | number | Slope × 3′ pairing score + intercept |
| `score.localAuContribution` | number | Slope × local-AU fraction + intercept (3 decimals) |
| `score.positionContribution` | number | Distance-to-UTR-end term (3 decimals) |
| `score.threePrimePairingScore` | number | Raw 3′ pairing score |
| `score.localAuFraction` | number | Position-weighted A/U fraction of the 30 nt flanks |
| `score.contextScore` | number | null | Sum of the four contributions; null when too close to the stop codon |

## Errors

| Code | Message |
|------|---------|
| 1001 | mRNA sequence cannot be null or empty |
| 1001 | miRNA sequence cannot be null or empty |
| 1001 | Site indices [s, e] are out of range for sequence of length L |
| 1001 | Unknown siteType (valid values listed) |

## Examples

### Example 1: 8mer in a poly-G 3′UTR

`GGG…(100 G)CUACCUCA(100 G)` with let-7a, site 100–107, `Seed8mer`: siteTypeContribution −0.310, positionContribution round(100·0.000172 − 0.07, 3) = −0.053, localAuFraction 0.5/(H30 + Σ1/(i+2)).

**Request:**
```json
{"mRnaSequence": "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCUACCUCAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "miRnaSequence": "UGAGGUAGUAGGUUGUAUAGUU", "siteStart": 100, "siteEnd": 107, "siteType": "Seed8mer"}
```

**Response (abridged):**
```json
{"scored": true, "score": {"type": "Seed8mer", "siteTypeContribution": -0.31, "positionContribution": -0.053}}
```

## Performance

- **Time Complexity:** O(n + m) (flanks of 30 nt, fixed-size alignment)
- **Space Complexity:** O(n)

## See Also

- [find_mirna_target_sites](find_mirna_target_sites.md)
- [analyze_target_context](analyze_target_context.md)
- [find_centered_sites](find_centered_sites.md)
