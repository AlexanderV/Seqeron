# find_degenerate_repeats

Find all maximal degenerate repeats (k-differences or k-mismatches; direct or palindromic) as Vmatch.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_degenerate_repeats` |
| **Method ID** | `RepeatFinder.FindDegenerateRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **all maximal degenerate repeats** with at most `maxDifferences` differences (Kurtz et al. 2000
REPuter; Vmatch `-e k` / `-h k`, `-p`, `-allmax`). `distance = "edit"` (default) uses unit-cost edit distance, so the
two instances may differ in length (`firstLength`, `secondLength`); `"hamming"` uses mismatches only.
`reverseComplement = true` reports palindromic repeats (instance 2 matched against the reverse complement of
instance 1; Vmatch `-p`). Edit-distance extension uses the greedy furthest-reaching fronts of Ukkonen 1985 / Myers
1986 (Vmatch `frontSEP.c`); the result is the complete maximal set (identical to Vmatch with its left-extension
shortcut disabled, and to a brute force of the definition — Evidence REP-DIRECT-001 §WP8). Non-ACGT symbols always
mismatch. Filters `maxLength`, `spacing = secondPosition − firstPosition − firstLength ≥ minSpacing`
(−2147483648 = complete Vmatch set). Sorted by (firstPosition, secondPosition, firstLength, secondLength).

`reporting = "bestPerSeed"` gives Vmatch's **default output** (no `-allmax`): for every exact seed, the best
extension by E-value, then identity, then length (Vmatch `cmpmatches`). Rows may repeat, one per seed, as in
Vmatch. `vmatchCompatible = true` reproduces stock Vmatch 2.3.1 exactly. That includes its left-extension seed
shortcut, which drops some maximal edit matches and changes best-per-seed choices. It also includes its
first-seed distance label in `-allmax` mode. The default is the complete, definition-correct output.
Cross-checked against real Vmatch with 0 differing rows (Evidence REP-DIRECT-001 §WP15).

## Core Documentation Reference

- Source: [RepeatFinder.cs#L5117](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L5117)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `minLength` | integer | No | Minimum length of each instance (Vmatch -l; > maxDifferences) (default 10) |
| `maxDifferences` | integer | No | Maximum distance k (default 1) |
| `distance` | string | No | Distance: edit (k-differences) or hamming (k-mismatches) (default edit) |
| `reverseComplement` | boolean | No | Palindromic (reverse-complement) repeats (Vmatch -p) (default false) |
| `maxLength` | integer | No | Maximum length of each instance (default 2147483647) |
| `minSpacing` | integer | No | Minimum spacing secondPosition - firstPosition - firstLength (default 1) |
| `reporting` | string | No | allMaximal (Vmatch -allmax) or bestPerSeed (Vmatch default output) (default allMaximal) |
| `vmatchCompatible` | boolean | No | Reproduce stock Vmatch 2.3.1 (seed shortcut) (default false) |

## Output Schema

`items`: `firstPosition, firstLength, secondPosition, secondLength, distance, spacing, firstCopy, secondCopy, isReverseComplement`

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | distance must be 'edit' or 'hamming' |
| 1003 | minLength must be >= 2 and > maxDifferences; maxDifferences >= 1; maxLength >= minLength |
| 1004 | reporting must be 'allMaximal' or 'bestPerSeed' |

## Examples

### Example 1: vmatch -l 6 -e 1 -allmax: insertions and deletions

**Input:** `{"sequence": "GATTACAGATTACATTTTTGATTCAGATTACA", "minLength": 6, "maxDifferences": 1, "distance": "edit", "reverseComplement": false, "maxLength": 2147483647, "minSpacing": -2147483648}`

**Output:**

```json
{"items": [{"firstPosition": 0, "firstLength": 7, "secondPosition": 6, "secondLength": 8, "distance": 1, "spacing": -1, "firstCopy": "GATTACA", "secondCopy": "AGATTACA", "isReverseComplement": false}, {"firstPosition": 0, "firstLength": 8, "secondPosition": 7, "secondLength": 8, "distance": 1, "spacing": -1, "firstCopy": "GATTACAG", "secondCopy": "GATTACAT", "isReverseComplement": false}, {"firstPosition": 0, "firstLength": 14, "secondPosition": 19, "secondLength": 13, "distance": 1, "spacing": 5, "firstCopy": "GATTACAGATTACA", "secondCopy": "GATTCAGATTACA", "isReverseComplement": false}, {"firstPosition": 5, "firstLength": 10, "secondPosition": 23, "secondLength": 9, "distance": 1, "spacing": 8, "firstCopy": "CAGATTACAT", "secondCopy": "CAGATTACA", "isReverseComplement": false}, {"firstPosition": 19, "firstLength": 6, "secondPosition": 25, "secondLength": 7, "distance": 1, "spacing": 0, "firstCopy": "GATTCA", "secondCopy": "GATTACA", "isReverseComplement": false}]}
```

### Example 2: vmatch -p -l 8 -h 1 -allmax: palindromic, Hamming

**Input:** `{"sequence": "TTGACCGTAACCCCCGTTACGGTCAACC", "minLength": 8, "maxDifferences": 1, "distance": "hamming", "reverseComplement": true, "maxLength": 2147483647, "minSpacing": -2147483648}`

**Output:**

```json
{"items": [{"firstPosition": 0, "firstLength": 12, "secondPosition": 14, "secondLength": 12, "distance": 1, "spacing": 2, "firstCopy": "TTGACCGTAACC", "secondCopy": "CGTTACGGTCAA", "isReverseComplement": true}, {"firstPosition": 13, "firstLength": 9, "secondPosition": 13, "secondLength": 9, "distance": 1, "spacing": -9, "firstCopy": "CCGTTACGG", "secondCopy": "CCGTTACGG", "isReverseComplement": true}]}
```

## Performance

- O(n log² n) seeding plus O(k²) front extension per seed and an output-sensitive maximality filter.

## See Also

- [find_approximate_direct_repeats](find_approximate_direct_repeats.md)
- [find_reverse_complement_repeats](find_reverse_complement_repeats.md)
- [find_direct_repeats](find_direct_repeats.md)
