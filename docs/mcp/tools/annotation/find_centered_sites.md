# find_centered_sites

Find centered miRNA target sites (Shin et al. 2010): >= 11 contiguous Watson-Crick pairs at miRNA nt 4-15 without a perfect seed match.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_centered_sites` |
| **Method ID** | `MiRnaAnalyzer.FindCenteredSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds centered sites (Shin et al. 2010, Mol Cell 38:789): 3′UTR sites that lack perfect seed pairing (miRNA nt 2–7 not all Watson–Crick) and instead pair the miRNA centre with at least 11 contiguous Watson–Crick pairs (nt 4–14, 5–15 or 4–15; G:U wobble does not count). `start`/`end` span the paired run (11 or 12 nt = `seedMatchLength`). TargetScan has no coefficients for centered sites, so `score` is 0. The additional Shin criterion “lacks 3′-compensatory pairing” is not tested. Use [find_mirna_target_sites](find_mirna_target_sites.md) for canonical seed sites.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.GrimsonContext.cs#L170](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.GrimsonContext.cs#L170)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `mRnaSequence` | string | Yes | mRNA nucleotide sequence (T→U applied internally) |
| `miRna` | object | Yes | Query miRNA (`name`, `sequence`, `seedSequence`, `seedStart`, `seedEnd`; >= 15 nt) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].start / end` | integer | 0-based inclusive span of the paired run |
| `sites[].type` | string | `Centered` |
| `sites[].seedMatchLength` | integer | 11 or 12 |
| `sites[].score` | number | 0 (no TargetScan coefficients) |
| `sites[].targetSequence / alignment / freeEnergy` |  | Duplex fields as in find_mirna_target_sites |

## Errors

| Code | Message |
|------|---------|
| 1001 | mRNA sequence cannot be null or empty |
| 1001 | miRNA cannot be null or empty |

## Examples

### Example 1: 12-nt run at nt 4-15

let-7a nt 4–15 fully paired inside a poly-G flank gives one site, start 10, end 21, seedMatchLength 12.

**Request:**
```json
{"mRnaSequence": "GGGGGGGGGGCAACCUACUACCGGGGGGGGGG", "miRna": {"name": "let-7a", "sequence": "UGAGGUAGUAGGUUGUAUAGUU", "seedSequence": "GAGGUAG", "seedStart": 1, "seedEnd": 7}}
```

**Response (abridged):**
```json
{"sites": [{"start": 10, "end": 21, "type": "Centered", "seedMatchLength": 12, "score": 0}]}
```

## Performance

- **Time Complexity:** O(n·m) for mRNA length n
- **Space Complexity:** O(1) per site

## See Also

- [find_mirna_target_sites](find_mirna_target_sites.md)
- [score_target_site_grimson](score_target_site_grimson.md)
