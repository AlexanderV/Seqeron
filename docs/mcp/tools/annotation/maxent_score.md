# maxent_score

Score a splice-site motif: Donor = MaxEntScan score5 (Yeo & Burge 2004); Acceptor = PWM log-weight sum.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `maxent_score` |
| **Method ID** | `SpliceSitePredictor.CalculateMaxEntScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

- **Donor:** the Yeo & Burge (2004) MaxEntScan **score5** maximum-entropy score in bits, computed by
  `SpliceSitePredictor.ScoreDonorMaxEnt` (me2x5 model, embedded `Data/maxent_score5.txt`). The motif must be
  the 9-nt donor window (3 exon + 6 intron nt, `GT` at 1-based positions 4-5) over A/C/G/T(/U), any case.
  Non-GT dinucleotides are scored (not rejected), as in `score5.pl`. Identical to maxentpy `score5`
  bit-for-bit on all 4^9 ACGT 9-mers; `score5.pl` prints the value to 2 decimals.
- **Acceptor:** sum of `log2(weight + 0.01)` over the acceptor PWM columns (T→U normalised; legacy heuristic).

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L930](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L931)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `motif` | string | Yes | Splice-site motif sequence; Donor: exactly 9 nt of A/C/G/T(/U) |
| `type` | string | Yes | `Donor` or `Acceptor` (case-insensitive) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | Donor: MaxEntScan score5 (bits); Acceptor: sum of `log2(weight + 0.01)` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Motif cannot be null or empty |
| 1001 | Type cannot be null or empty |
| 1002 | Unknown splice-site type '&lt;value&gt;'. Expected 'Donor' or 'Acceptor'. |
| 1001 | MaxEntScan score5 requires a window of exactly 9 nt; got &lt;n&gt;. (Donor) |
| 1001 | MaxEntScan score5 requires an A/C/G/T(/U) sequence; found '&lt;c&gt;'. (Donor) |

## Examples

### Example 1: Canonical donor (MaxEntScan reference example)

`CAGGUAAGU` (= `cagGTAAGT`, the score5.pl example, printed 10.86):

**Response:**
```json
{ "score": 10.858313101356437 }
```

### Example 2: Non-GT (GC) donor

`CAGGCAAGU` (maxentpy / score5.pl: 3.10):

**Response:**
```json
{ "score": 3.104161570996358 }
```

## Performance

- **Time Complexity:** Donor O(1) (one table lookup); Acceptor O(L) for motif length L
- **Space Complexity:** O(1)

## See Also

- [find_donor_sites](find_donor_sites.md) / [find_acceptor_sites](find_acceptor_sites.md) — site scanning
