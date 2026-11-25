# maxent_score

Score a splice-site motif: Donor = MaxEntScan score5, Acceptor = MaxEntScan score3 (Yeo & Burge 2004).

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
- **Acceptor:** the Yeo & Burge (2004) MaxEntScan **score3** maximum-entropy score in bits, computed by
  `SpliceSitePredictor.ScoreAcceptorMaxEnt` (me2x3acc1-9 model, embedded `Data/maxent_score3.txt`). The motif must
  be the 23-nt acceptor window (20 intron + 3 exon nt, `AG` at 1-based positions 19-20) over A/C/G/T(/U), any case.
  Non-AG dinucleotides are scored (not rejected), as in `score3.pl`. Identical to maxentpy `score3` bit-for-bit on
  200,025 test windows (100,000 random + 100,000 AG 23-mers + all 16 dinucleotides); `score3.pl` prints 2 decimals.
- U is read as T (extension: the references accept only A/C/G/T); other `type` values (`Branch`, `U12Donor`,
  `U12Acceptor`) return 0.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L935](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L935)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `motif` | string | Yes | Splice-site motif sequence; Donor: exactly 9 nt, Acceptor: exactly 23 nt, of A/C/G/T(/U) |
| `type` | string | Yes | `Donor` or `Acceptor` (case-insensitive) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | Donor: MaxEntScan score5 (bits); Acceptor: MaxEntScan score3 (bits) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Motif cannot be null or empty |
| 1001 | Type cannot be null or empty |
| 1002 | Unknown splice-site type '&lt;value&gt;'. Expected 'Donor' or 'Acceptor'. |
| 1001 | MaxEntScan score5 requires a window of exactly 9 nt; got &lt;n&gt;. (Donor) |
| 1001 | MaxEntScan score5 requires an A/C/G/T(/U) sequence; found '&lt;c&gt;'. (Donor) |
| 1001 | MaxEntScan score3 requires a window of exactly 23 nt; got &lt;n&gt;. (Acceptor) |
| 1001 | MaxEntScan score3 requires an A/C/G/T(/U) sequence; found '&lt;c&gt;'. (Acceptor) |

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

### Example 3: Acceptor (MaxEntScan reference example)

`TTCCAAACGAACTTTTGTAGGGA` (= `ttccaaacgaacttttgtAGgga`, the score3 example, printed 2.89):

**Response:**
```json
{ "score": 2.8867730651152104 }
```

## Performance

- **Time Complexity:** O(1) (Donor: one table lookup; Acceptor: nine sub-table lookups)
- **Space Complexity:** O(1)

## See Also

- [find_donor_sites](find_donor_sites.md) / [find_acceptor_sites](find_acceptor_sites.md) — site scanning
