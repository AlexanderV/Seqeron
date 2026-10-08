# mask_approximate_tandem_repeats

Mask Tandem Repeats Finder repeats (TRF -m masked sequence).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mask_approximate_tandem_repeats` |
| **Method ID** | `RepeatFinder.MaskApproximateTandemRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns TRF's masked sequence (`trf ... -m`; TRF README: "every location that occurred in a tandem repeat
changed to the letter 'N'"): every position inside a repeat reported by `find_approximate_tandem_repeats` with the
same TRF parameters becomes `N`, or — with `softMask` — lower case. Positions outside repeats keep the input
characters (TRF itself upper-cases its whole output). Same length as the input. Defaults = TRF recommended
`2 7 7 80 10 50 500`. `apparentSizeTable` / `apparentSizeTableKind` optionally replace the exact apparent-size table,
as in `find_approximate_tandem_repeats` (TRF's own table → TRF's mask bit for bit).

## Core Documentation Reference

- Source: [RepeatFinder.cs#L1039](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L1039)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence (case-insensitive; any non-A/C/G/T symbol never matches) |
| `maxPeriod` | integer | No | TRF MaxPeriod (default 500) |
| `minScore` | integer | No | TRF Minscore (default 50) |
| `matchWeight` | integer | No | TRF Match weight (default 2) |
| `mismatchPenalty` | integer | No | TRF Mismatch penalty (default 7) |
| `indelPenalty` | integer | No | TRF Delta (indel penalty) (default 7) |
| `matchProbability` | integer | No | TRF PM in percent (80 or 75) (default 80) |
| `indelProbability` | integer | No | TRF PI in percent (default 10) |
| `maxRepeatLength` | integer | No | TRF -l: maximum tandem-repeat length in bp (default 2000000) |
| `eliminateRedundancy` | boolean | No | Redundancy elimination (TRF -r turns it off) (default true) |
| `softMask` | boolean | No | Lower-case repeat positions instead of writing N (default false) |
| `apparentSizeTable` | string | No | 2001 comma/space-separated integers y(d), d = 0..2000 (entry 0 ignored; each in 0..max(d,20)−1); empty = exact table |
| `apparentSizeTableKind` | string | No | `apparent` (y(d), default) or `trfWaitingTimes` (TRF `waitdata` w(d); y = max(d,20) − w − 1) |

## Output Schema

`masked`: the masked sequence (same length as the input)

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | TRF parameter out of range (weights >= 1, PM 80/75, PI 1-100, MaxPeriod 1-2000, -l >= 1) |
| 1004 | apparentSizeTable must hold 2001 integers, each in 0..max(d,20)-1 (waiting times likewise); apparentSizeTableKind must be 'apparent' or 'trfWaitingTimes' |

## Examples

### Example 1: Hard mask = TRF -m file

**Input:** `{"sequence": "CACACACACACACACACGCACACACACACCGACACACACACACACACACACACACACACAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC"}`

**Output:**

```json
{"masked": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC"}
```

### Example 2: Soft mask

**Input:** `{"sequence": "CACACACACACACACACGCACACACACACCGACACACACACACACACACACACACACACAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC", "softMask": true}`

**Output:**

```json
{"masked": "cacacacacacacacacgcacacacacaccgacacacacacacacacacacacacacacaAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC"}
```

## Performance

- As find_approximate_tandem_repeats, plus O(n) masking.

## See Also

- [find_approximate_tandem_repeats](find_approximate_tandem_repeats.md)
- [mask_low_complexity](mask_low_complexity.md)
