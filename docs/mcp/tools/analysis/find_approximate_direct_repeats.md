# find_approximate_direct_repeats

Find all maximal k-mismatch direct repeats (REPuter / Vmatch -h k).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_approximate_direct_repeats` |
| **Method ID** | `RepeatFinder.FindApproximateDirectRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **all maximal k-mismatch (Hamming) direct repeats** (Kurtz et al. 2000 REPuter; Vmatch `-h k`): two
copies of equal length `length` at 0-based `firstPosition < secondPosition` with at most `maxMismatches`
mismatches that cannot be extended on either side within k. Seeds are exact maximal pairs of length
⌊minLength/(k+1)⌋ (pigeonhole), extended to every maximal window. Default per-diagonal maximality (for k = 0 equal
to `find_direct_repeats`); `excludeContained = true` also drops repeats contained in a k-mismatch repeat on another
diagonal — output identical to `vmatch -l minLength -h k -allmax`. Non-ACGT symbols are mismatches. Filters
`length ≤ maxLength`, `spacing = secondPosition − firstPosition − length ≥ minSpacing`. Sorted by (firstPosition,
secondPosition, length). For edit distance or palindromic repeats use `find_degenerate_repeats`.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L4126](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L4126)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `minLength` | integer | No | Minimum repeat length (> maxMismatches) (default 10) |
| `maxMismatches` | integer | No | Maximum Hamming distance k (default 1) |
| `maxLength` | integer | No | Maximum repeat length (default 2147483647) |
| `minSpacing` | integer | No | Minimum number of bases between the copies (negative admits overlap) (default 1) |
| `excludeContained` | boolean | No | Drop repeats contained in a k-mismatch repeat on another diagonal (vmatch -allmax) (default false) |

## Output Schema

`items`: `firstPosition, secondPosition, length, mismatches, spacing, firstCopy, secondCopy`

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | minLength must be >= 2 and > maxMismatches; maxMismatches >= 0; maxLength >= minLength |

## Examples

### Example 1: vmatch -l 10 -h 2: one repeat with 2 mismatches

**Input:** `{"sequence": "ACGTTGCAAGCTTACGGGGGGACGATGCAAGCATACGG", "minLength": 10, "maxMismatches": 2, "maxLength": 2147483647, "minSpacing": -2147483648}`

**Output:**

```json
{"items": [{"firstPosition": 0, "secondPosition": 21, "length": 17, "mismatches": 2, "spacing": 4, "firstCopy": "ACGTTGCAAGCTTACGG", "secondCopy": "ACGATGCAAGCATACGG"}]}
```

### Example 2: Boundary-limited exact repeat (0 mismatches)

**Input:** `{"sequence": "GATTACAGATTACA", "minLength": 6, "maxMismatches": 1, "maxLength": 2147483647, "minSpacing": -2147483648}`

**Output:**

```json
{"items": [{"firstPosition": 0, "secondPosition": 7, "length": 7, "mismatches": 0, "spacing": 0, "firstCopy": "GATTACA", "secondCopy": "GATTACA"}]}
```

### Example 3: vmatch -l 5 -h 1 -allmax

**Input:** `{"sequence": "AAAAAAAACGTTGCAACGTAAAA", "minLength": 5, "maxMismatches": 1, "maxLength": 2147483647, "minSpacing": -2147483648, "excludeContained": true}`

**Output:**

```json
{"items": [{"firstPosition": 0, "secondPosition": 1, "length": 8, "mismatches": 1, "spacing": -7, "firstCopy": "AAAAAAAA", "secondCopy": "AAAAAAAC"}, {"firstPosition": 0, "secondPosition": 18, "length": 5, "mismatches": 1, "spacing": 13, "firstCopy": "AAAAA", "secondCopy": "TAAAA"}, {"firstPosition": 1, "secondPosition": 18, "length": 5, "mismatches": 1, "spacing": 12, "firstCopy": "AAAAA", "secondCopy": "TAAAA"}, {"firstPosition": 2, "secondPosition": 18, "length": 5, "mismatches": 1, "spacing": 11, "firstCopy": "AAAAA", "secondCopy": "TAAAA"}, {"firstPosition": 3, "secondPosition": 18, "length": 5, "mismatches": 1, "spacing": 10, "firstCopy": "AAAAA", "secondCopy": "TAAAA"}, {"firstPosition": 5, "secondPosition": 13, "length": 6, "mismatches": 1, "spacing": 2, "firstCopy": "AAACGT", "secondCopy": "CAACGT"}, {"firstPosition": 6, "secondPosition": 14, "length": 6, "mismatches": 1, "spacing": 2, "firstCopy": "AACGTT", "secondCopy": "AACGTA"}]}
```

## Performance

- O(n log² n + s·k + z) for s seeds and z reported repeats.

## See Also

- [find_direct_repeats](find_direct_repeats.md)
- [find_degenerate_repeats](find_degenerate_repeats.md)
- [find_supermaximal_repeats](find_supermaximal_repeats.md)
