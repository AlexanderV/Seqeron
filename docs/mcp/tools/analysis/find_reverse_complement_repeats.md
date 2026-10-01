# find_reverse_complement_repeats

Find maximal exact reverse-complement repeat pairs (MUMmer repeat-match reverse / Vmatch -p).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_reverse_complement_repeats` |
| **Method ID** | `RepeatFinder.FindReverseComplementRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **maximal exact reverse-complement repeat pairs**: copy 1 at `firstPosition` and copy 2 at
`secondPosition` (both 0-based, `firstPosition ≤ secondPosition`) with copy 2 = reverse complement of copy 1, not
extendable outward or inward — the reverse (`r`) lines of MUMmer `repeat-match` and Vmatch `-p` (repeat-match
`Start2` = `secondPosition + length`). One pair can carry several maximal lengths on different anti-diagonals.
Filters: `minLength ≤ length ≤ maxLength` (longer maximal pairs are not truncated) and
`spacing = secondPosition − firstPosition − length ≥ minSpacing` (negative admits overlap; `minSpacing =
−2147483648` with `maxLength = 2147483647` returns the complete repeat-match set). Only A/C/G/T pair
(case-insensitive). Sorted by (firstPosition, secondPosition, length).

## Core Documentation Reference

- Source: [RepeatFinder.cs#L3961](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L3961)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `minLength` | integer | No | Minimum repeat length (default 5) |
| `maxLength` | integer | No | Maximum repeat length (>= minLength) (default 50) |
| `minSpacing` | integer | No | Minimum number of bases between the copies (negative admits overlap) (default 1) |

## Output Schema

`items`: `firstPosition, secondPosition, repeatSequence, secondSequence, length, spacing`

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | minLength must be >= 2 |
| 1003 | maxLength must be >= minLength |

## Examples

### Example 1: Defaults: separated copies only

**Input:** `{"sequence": "TTGCATGCAAAAAATTTTTTTGCATGCAA"}`

**Output:**

```json
{"items": [{"firstPosition": 0, "secondPosition": 15, "repeatSequence": "TTGCATGCAAAAAA", "secondSequence": "TTTTTTGCATGCAA", "length": 14, "spacing": 1}, {"firstPosition": 8, "secondPosition": 14, "repeatSequence": "AAAAA", "secondSequence": "TTTTT", "length": 5, "spacing": 1}, {"firstPosition": 9, "secondPosition": 16, "repeatSequence": "AAAAA", "secondSequence": "TTTTT", "length": 5, "spacing": 2}]}
```

### Example 2: EcoRI sites, case-insensitive

**Input:** `{"sequence": "gaattcAAAAAgaattc", "minLength": 6}`

**Output:**

```json
{"items": [{"firstPosition": 0, "secondPosition": 11, "repeatSequence": "GAATTC", "secondSequence": "GAATTC", "length": 6, "spacing": 5}]}
```

### Example 3: Complete repeat-match set (8 19r 12 / 7 12r 6 / 16 19r 4)

**Input:** `{"sequence": "AAAAAAAACGTTGCAACGTAAAA", "minLength": 3, "maxLength": 2147483647, "minSpacing": -2147483648}`

**Output:**

```json
{"items": [{"firstPosition": 6, "secondPosition": 6, "repeatSequence": "AACGTT", "secondSequence": "AACGTT", "length": 6, "spacing": -6}, {"firstPosition": 7, "secondPosition": 7, "repeatSequence": "ACGTTGCAACGT", "secondSequence": "ACGTTGCAACGT", "length": 12, "spacing": -12}, {"firstPosition": 15, "secondPosition": 15, "repeatSequence": "ACGT", "secondSequence": "ACGT", "length": 4, "spacing": -4}]}
```

## Performance

- O(n log² n + z) (suffix array + LCP of S·#·revcomp(S); z = reported pairs). Space O(n).

## See Also

- [find_direct_repeats](find_direct_repeats.md)
- [find_inverted_repeats](find_inverted_repeats.md)
- [find_degenerate_repeats](find_degenerate_repeats.md)
