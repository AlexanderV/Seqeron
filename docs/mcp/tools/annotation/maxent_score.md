# maxent_score

Compute a MaxEntScan-like log-likelihood score for a splice-site motif.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `maxent_score` |
| **Method ID** | `SpliceSitePredictor.CalculateMaxEntScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scores a donor or acceptor splice motif (T→U normalised) by summing `log2(weight + 0.01)` over the
positions of the corresponding position weight matrix (donor 9-mer around GT, acceptor window around AG).
Higher scores indicate a stronger match to the consensus; the value is a MaxEntScan-style log-likelihood,
not bounded to [0,1].

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L937](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L937)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `motif` | string | Yes | Splice-site motif sequence |
| `type` | string | Yes | `Donor` or `Acceptor` (case-insensitive) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | Sum of `log2(weight + 0.01)` over motif positions |

## Errors

| Code | Message |
|------|---------|
| 1001 | Motif cannot be null or empty |
| 1001 | Type cannot be null or empty |
| 1002 | Unknown splice-site type '&lt;value&gt;'. Expected 'Donor' or 'Acceptor'. |

## Examples

### Example 1: Perfect donor consensus

`CAGGUAAGU` matches every donor PWM position (weight 1.0), so score = `9 · log2(1.01)`:

**Response:**
```json
{ "score": 0.12919763679363047 }
```

### Example 2: Non-consensus donor

Breaking the invariant positions drops the score sharply:

**Response:**
```json
{ "score": -13.187225328709957 }
```

## Performance

- **Time Complexity:** O(L) for motif length L
- **Space Complexity:** O(1)

## See Also

- [find_donor_sites](find_donor_sites.md) / [find_acceptor_sites](find_acceptor_sites.md) — site scanning
