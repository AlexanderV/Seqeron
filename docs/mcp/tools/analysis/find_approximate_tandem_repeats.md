# find_approximate_tandem_repeats

Find approximate tandem repeats with the Tandem Repeats Finder (TRF) model.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_approximate_tandem_repeats` |
| **Method ID** | `RepeatFinder.FindApproximateTandemRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds approximate (imperfect / interrupted) tandem repeats with the **Tandem Repeats Finder** model
(Benson 1999; TRF 4.10.0 `trf File Match Mismatch Delta PM PI Minscore MaxPeriod [-l] [-r] [-f]`). Defaults are
TRF's recommended `2 7 7 80 10 50 500`. Detection follows TRF (k-tuple matches, sum-of-heads and apparent-size
criteria, random-walk distance range for d > 20, best-period list, wraparound DP alignment, narrow band for
patterns > 20, redundancy elimination); candidate distances are examined up to TRF's MAXDISTANCE and `maxPeriod`
filters the reported periods. Each item carries TRF's `.dat` columns — `start` (0-based; TRF prints 1-based),
`spanLength`, `period`, `copyNumber`, `consensusSize`, `percentMatches`, `percentIndels` (exact; TRF truncates),
`alignmentScore`, base percentages, `entropyTrf` (TRF's entropy column) and `entropy` (A/C/G/T-only) — plus the
final alignment rows and, when `flankLength > 0`, the flanking sequences. `minPeriod` is a library filter applied
after redundancy elimination. `examineUpToMaxPeriodOnly = true` (legacy library mode) examines candidate
distances only up to `maxPeriod` (faster for short periods) and requires the recommended weights/PM/PI and default
`-l/-r/-f`. Measured against the TRF binary: 99.8 % identical rows (Evidence REP-APPROX-001).

## Core Documentation Reference

- Source: [RepeatFinder.cs#L1002](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L1002)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence (case-insensitive; any non-A/C/G/T symbol never matches) |
| `minPeriod` | integer | No | Minimum reported period (default 1) |
| `maxPeriod` | integer | No | TRF MaxPeriod (default 500) |
| `minScore` | integer | No | TRF Minscore (default 50) |
| `matchWeight` | integer | No | TRF Match weight (default 2) |
| `mismatchPenalty` | integer | No | TRF Mismatch penalty (default 7) |
| `indelPenalty` | integer | No | TRF Delta (indel penalty) (default 7) |
| `matchProbability` | integer | No | TRF PM in percent (80 or 75) (default 80) |
| `indelProbability` | integer | No | TRF PI in percent (default 10) |
| `maxRepeatLength` | integer | No | TRF -l: maximum tandem-repeat length in bp (default 2000000) |
| `eliminateRedundancy` | boolean | No | Redundancy elimination (TRF -r turns it off) (default true) |
| `flankLength` | integer | No | TRF -f: flanking-sequence length on each side (0 = none; TRF -f uses 500) (default 0) |
| `examineUpToMaxPeriodOnly` | boolean | No | Examine candidate distances only up to maxPeriod (legacy mode; recommended weights required) (default false) |

## Output Schema

`items`: `start, spanLength, period, consensusSize, consensus, copyNumber, percentMatches, percentIndels, alignmentScore, percentA, percentC, percentG, percentT, entropy, entropyTrf, alignedSequence, alignedConsensus` (+ `leftFlank, rightFlank` when `flankLength > 0`)

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | examineUpToMaxPeriodOnly requires the TRF recommended weights and default -l/-r/-f |
| 1003 | TRF parameter out of range (weights >= 1, PM 80/75, PI 1-100, MaxPeriod 1-2000, -l >= 1, flankLength >= 0, minPeriod >= 1) |

## Examples

### Example 1: TRF recommended parameters (TRF row 61 124 7 9.1 7 92 0 110 14 14 26 42 1.83 TCATTGG)

**Input:** `{"sequence": "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTCATTGGTCATTGGTCANTGGTCATTGGTCATTGGTCATTNGTCATTGGTCATTGGTCATTGGTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA"}`

**Output:**

```json
{"items": [{"start": 60, "spanLength": 64, "period": 7, "consensusSize": 7, "consensus": "TCATTGG", "copyNumber": 9.142857142857142, "percentMatches": 92.98245614035088, "percentIndels": 0, "alignmentScore": 110, "percentA": 14.0625, "percentC": 14.0625, "percentG": 26.5625, "percentT": 42.1875, "entropy": 1.8424627477154065, "entropyTrf": 1.829258111162015, "alignedSequence": "TCATTGGTCATTGGTCANTGGTCATTGGTCATTGGTCATTNGTCATTGGTCATTGGTCATTGGT", "alignedConsensus": "TCATTGGTCATTGGTCATTGGTCATTGGTCATTGGTCATTGGTCATTGGTCATTGGTCATTGGT"}]}
```

### Example 2: trf 2 3 5 80 10 40 200 (TRF row 1 60 2 30.0 2 89 0 105 CA)

**Input:** `{"sequence": "CACACACACACACACACGCACACACACACCGACACACACACACACACACACACACACACAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC", "maxPeriod": 200, "minScore": 40, "mismatchPenalty": 3, "indelPenalty": 5}`

**Output:**

```json
{"items": [{"start": 0, "spanLength": 60, "period": 2, "consensusSize": 2, "consensus": "CA", "copyNumber": 30, "percentMatches": 89.65517241379311, "percentIndels": 0, "alignmentScore": 105, "percentA": 46.666666666666664, "percentC": 50, "percentG": 3.3333333333333335, "percentT": 0, "entropy": 1.1766796675107107, "entropyTrf": 1.1766796675107107, "alignedSequence": "CACACACACACACACACGCACACACACACCGACACACACACACACACACACACACACACA", "alignedConsensus": "CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"}]}
```

### Example 3: Legacy short-period scan

**Input:** `{"sequence": "CAGCAGCAGTAGCAGCAG", "maxPeriod": 6, "minScore": 20, "examineUpToMaxPeriodOnly": true}`

**Output:**

```json
{"items": [{"start": 0, "spanLength": 18, "period": 3, "consensusSize": 3, "consensus": "CAG", "copyNumber": 6, "percentMatches": 86.66666666666667, "percentIndels": 0, "alignmentScore": 27, "percentA": 33.333333333333336, "percentC": 27.77777777777778, "percentG": 33.333333333333336, "percentT": 5.555555555555555, "entropy": 1.8016366412706073, "entropyTrf": 1.8016366412706075, "alignedSequence": "CAGCAGCAGTAGCAGCAG", "alignedConsensus": "CAGCAGCAGCAGCAGCAG"}]}
```

## Performance

- O(n · MAXDISTANCE) k-tuple scan plus one wraparound DP per examined candidate (O(region · pattern), narrow band for patterns > 20).

## See Also

- [mask_approximate_tandem_repeats](mask_approximate_tandem_repeats.md)
- [tandem_repeat_bernoulli_statistics](tandem_repeat_bernoulli_statistics.md)
- [find_microsatellites](find_microsatellites.md)
- [find_tandem_repeats](find_tandem_repeats.md)
