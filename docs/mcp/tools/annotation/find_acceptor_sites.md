# find_acceptor_sites

Find candidate 3' (acceptor) splice sites.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_acceptor_sites` |
| **Method ID** | `SpliceSitePredictor.FindAcceptorSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a sequence for 3' splice acceptor sites at every AG dinucleotide, scoring the upstream
polypyrimidine tract and `(Y)nNCAG|G` consensus (Shapiro & Senapathy 1987; Burge et al. 1999). Sites
scoring at least `minScore` are reported with the 0-based position of the G in the AG, motif, PWM score,
and confidence. Optionally includes non-canonical U12 AC acceptor motifs.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L264](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L264)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to scan |
| `minScore` | number | No | Minimum PWM score [0,1] (default 0.5) |
| `includeNonCanonical` | boolean | No | Include U12 AC acceptors (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].position` | integer | 0-based index of the G of the AG |
| `sites[].type` | string | `Acceptor` |
| `sites[].motif` | string | Matched motif |
| `sites[].score` | number | PWM match score [0,1] |
| `sites[].confidence` | number | Confidence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Strong polypyrimidine tract + CAG

`UUUUUUUUUUUUUUUUCAGGG` scores a strong acceptor at the AG:

**Response:**
```json
{ "sites": [ { "position": 18, "type": "Acceptor", "score": 0.8393 } ] }
```

### Example 2: No AG dinucleotide

```json
{ "sites": [] }
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) for k sites

## See Also

- [find_donor_sites](find_donor_sites.md) — 5' donor sites
- [find_branch_points](find_branch_points.md) — intronic branch points
- [predict_introns](predict_introns.md) — pair donors and acceptors into introns
