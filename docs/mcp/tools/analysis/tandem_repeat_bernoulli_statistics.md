# tandem_repeat_bernoulli_statistics

Estimate the TRF Bernoulli-model PM/PI of a tandem-repeat tract between adjacent copies.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `tandem_repeat_bernoulli_statistics` |
| **Method ID** | `RepeatFinder.ComputeBernoulliStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Estimates the Tandem Repeats Finder Bernoulli-model parameters of a tandem-repeat tract (Benson 1999): **PM**
(match probability) and **PI** (indel probability) between ADJACENT copies. The tract is analysed as TRF analyses a
detected repeat — wraparound-DP alignment against tandem copies of the candidate pattern (the last `period`
bases), majority-rule consensus, realignment, then comparison of each copy with the next one ("between adjacent
copies in the sequence, not between the sequence and the consensus pattern"). Each compared column is one
Bernoulli trial (heads = match). PM/PI therefore equal TRF's `% matches / % indels` (/100) for the same region.
Flanks that do not align are ignored; a tract without two aligned copies gives zero trials and PM = PI = 0.
`meetsExpectedMatchProbability` = PM ≥ `expectedMatchProbability`.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L3528](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L3528)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `repeatTract` | string | Yes | The tandem-repeat tract (>= 2 x period symbols; case-insensitive) |
| `period` | integer | Yes | Candidate period of the tract |
| `expectedMatchProbability` | number | No | PM the tract is compared with (Benson 1999 default 0.80) (default 0.8) |

## Output Schema

`period, adjacentCopyPairs, bernoulliTrials, matches, mismatches, indels, matchProbability, indelProbability, percentMatches, percentIndels, expectedMatches, meetsExpectedMatchProbability`

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | repeatTract must hold at least 2 x period symbols |
| 1003 | period must be 1-2000; expectedMatchProbability in [0, 1] |

## Examples

### Example 1: Perfect (CA)5

**Input:** `{"repeatTract": "CACACACACA", "period": 2}`

**Output:**

```json
{"period": 2, "adjacentCopyPairs": 4, "bernoulliTrials": 8, "matches": 8, "mismatches": 0, "indels": 0, "matchProbability": 1, "indelProbability": 0, "percentMatches": 100, "percentIndels": 0, "expectedMatches": 8, "meetsExpectedMatchProbability": true}
```

### Example 2: CAG x6 with one TAG copy (TRF m=13 mm=2 ind=0)

**Input:** `{"repeatTract": "CAGCAGCAGTAGCAGCAG", "period": 3}`

**Output:**

```json
{"period": 3, "adjacentCopyPairs": 5, "bernoulliTrials": 15, "matches": 13, "mismatches": 2, "indels": 0, "matchProbability": 0.8666666666666667, "indelProbability": 0, "percentMatches": 86.66666666666667, "percentIndels": 0, "expectedMatches": 13, "meetsExpectedMatchProbability": true}
```

### Example 3: Deletion tract (TRF m=25 mm=0 ind=2)

**Input:** `{"repeatTract": "CAGCAGCAGCAGCAGAGCAGCAGCAGCAG", "period": 3}`

**Output:**

```json
{"period": 3, "adjacentCopyPairs": 9, "bernoulliTrials": 27, "matches": 25, "mismatches": 0, "indels": 2, "matchProbability": 0.9259259259259259, "indelProbability": 0.07407407407407407, "percentMatches": 92.5925925925926, "percentIndels": 7.4074074074074066, "expectedMatches": 25, "meetsExpectedMatchProbability": true}
```

## Performance

- O(tract · period) wraparound DP.

## See Also

- [find_approximate_tandem_repeats](find_approximate_tandem_repeats.md)
- [find_microsatellites](find_microsatellites.md)
