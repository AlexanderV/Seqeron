# find_mirna_target_sites

Find miRNA target sites in an mRNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_mirna_target_sites` |
| **Method ID** | `MiRnaAnalyzer.FindTargetSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans an mRNA (DNA T→U applied internally) for canonical miRNA target sites keyed on the seed reverse
complement, classifying each into the Bartel (2009) / TargetScan site hierarchy — **8mer**, **7mer-m8**,
**7mer-A1**, **6mer**, or **offset-6mer** — and scoring each by site type and pairing. Only sites at or
above `minScore` are returned. 8mer sites (full seed match plus an A opposite miRNA position 1) score
highest.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L157](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L157)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `mRnaSequence` | string | Yes | mRNA nucleotide sequence |
| `miRna` | object | Yes | Query miRNA record |
| `minScore` | number | No | Minimum site score to report (default 0.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].start` | integer | 0-based site start |
| `sites[].end` | integer | 0-based site end |
| `sites[].targetSequence` | string | mRNA subsequence at the site |
| `sites[].miRnaName` | string | Query miRNA name |
| `sites[].type` | string | Seed8mer / Seed7merM8 / Seed7merA1 / Seed6mer / Offset6mer / Supplementary / Centered |
| `sites[].seedMatchLength` | integer | Length of the seed match |
| `sites[].score` | number | Site score |
| `sites[].freeEnergy` | number | Predicted duplex free energy |
| `sites[].alignment` | string | miRNA–target alignment string |

## Errors

| Code | Message |
|------|---------|
| 1001 | mRNA sequence cannot be null or empty |
| 1001 | miRNA cannot be null or empty |

## Examples

### Example 1: 8mer site

For let-7a (seed RC `CUACCUC`), an mRNA containing `CUACCUC` + `A` yields an 8mer site:

**Response:**
```json
{ "sites": [ { "type": "Seed8mer", "seedMatchLength": 8 } ] }
```

### Example 2: No seed match

A poly-A mRNA has no let-7a target site:

**Response:**
```json
{ "sites": [] }
```

## Performance

- **Time Complexity:** O(n) scan over the mRNA of length n (with local extension)
- **Space Complexity:** O(k) for k reported sites

## See Also

- [create_mirna](create_mirna.md) — build the query miRNA record
- [analyze_target_context](analyze_target_context.md) — AU/context score around a site
- [align_mirna_to_target](align_mirna_to_target.md) — full duplex alignment
