# compute_ta_3utr

TA_3UTR (target-site abundance) of a miRNA seed over a set of 3′UTRs: site count N and log10(N).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `compute_ta_3utr` |
| **Method ID** | `MiRnaAnalyzer.ComputeTa3Utr / MiRnaAnalyzer.CountSeedSites3Utr` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

TargetScan context++ feature TA_3UTR (Garcia et al. 2011; Agarwal et al. 2015): `ta = log10(N)` where `N` is the total number of non-overlapping canonical 8mer, 7mer-m8 and 7mer-A1 sites of the miRNA seed (nt 2–8) across the supplied 3′UTRs (bare 6mers and offset 6mers are not counted; the scan is left-to-right greedy). Returns `siteCount` (N, `CountSeedSites3Utr`) and `ta` (0 when N = 0). Pass `ta` to `score_target_site_context_plus_plus`.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L1358](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L1358)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `miRnaSequence` | string | Yes | miRNA nucleotide sequence (≥ 8 nt; seed = nt 2–8) |
| `threePrimeUtrs` | array | Yes | 3′UTR sequences (RNA or DNA, any case) that define the transcriptome |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `siteCount` | integer | Total non-overlapping 8mer + 7mer-m8 + 7mer-A1 sites N |
| `ta` | number | log10(N); 0 when N = 0 |
| `utrCount` | integer | Number of supplied UTRs |

## Errors

| Code | Message |
|------|---------|
| 1001 | miRNA sequence must have at least 8 nt (seed region nt 2-8) |
| 1001 | threePrimeUtrs cannot be null |

## Examples

### Example 1: Four sites over four UTRs

8mer; 7mer-m8 + 7mer-A1; a bare 6mer (not counted); a lower-case DNA 8mer → N = 4, TA = log10(4).

**Request:**
```json
{"miRnaSequence": "UGAGGUAGUAGGUUGUAUAGUU", "threePrimeUtrs": ["GGGCUACCUCAGGG", "GGGCUACCUCGGGGGUACCUCAGG", "GGGUACCUCGGG", "gggctacctcaggg"]}
```

**Response (abridged):**
```json
{"siteCount": 4, "ta": 0.6020599913279624, "utrCount": 4}
```

## Performance

- **Time Complexity:** O(total UTR length)
- **Space Complexity:** O(1)

## See Also

- [score_target_site_context_plus_plus](score_target_site_context_plus_plus.md)
- [find_mirna_target_sites](find_mirna_target_sites.md)
