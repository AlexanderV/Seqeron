# find_donor_sites

Find candidate 5' (donor) splice sites.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_donor_sites` |
| **Method ID** | `SpliceSitePredictor.FindDonorSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a sequence for 5' splice donor sites at every GU/GT dinucleotide, scoring a 9-mer window against the
`MAG|GURAGU` consensus (Shapiro & Senapathy 1987). Sites scoring at least `minScore` are reported with the
0-based dinucleotide position, motif, PWM score, and confidence. Optionally includes non-canonical GC and
U12 AT/AU donor motifs.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L204](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L204)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to scan |
| `minScore` | number | No | Minimum PWM score [0,1] (default 0.5) |
| `includeNonCanonical` | boolean | No | Include GC / U12 AT/AU donors (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].position` | integer | 0-based GU/GT dinucleotide index |
| `sites[].type` | string | `Donor` |
| `sites[].motif` | string | Matched motif |
| `sites[].score` | number | PWM match score [0,1] |
| `sites[].confidence` | number | Confidence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Perfect GU consensus

`CAGGUAAGU` matches `MAG|GURAGU` at all 9 positions:

**Response:**
```json
{ "sites": [ { "position": 3, "type": "Donor", "score": 1.0, "confidence": 1.0 } ] }
```

### Example 2: No GU dinucleotide

```json
{ "sites": [] }
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) for k sites

## See Also

- [find_acceptor_sites](find_acceptor_sites.md) — 3' acceptor sites
- [predict_introns](predict_introns.md) — pair donors and acceptors into introns
- [maxent_score](maxent_score.md) — MaxEntScan-like motif score
