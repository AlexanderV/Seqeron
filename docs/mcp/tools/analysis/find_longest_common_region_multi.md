# find_longest_common_region_multi

Longest common substring of several sequences with a minimum support (k-common substring problem).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_longest_common_region_multi` |
| **Method ID** | `GenomicAnalyzer.FindLongestCommonRegion` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The longest string that occurs in at least `minSupport` of the `N >= 2` input sequences (Gusfield 1997 §7.6; Rosalind LCSM for `minSupport = N`), computed with the canonical generalized suffix tree (`SuffixTree.FindLongestCommonSubstrings`). Every distinct longest string is returned in ordinal order, with its first 0-based occurrence in each sequence (`-1` where a sequence does not contain it: possible when `minSupport < N`) and its support. Only A/C/G/T match (case-insensitive); an `N` or another IUPAC code ends a region. The two-sequence tool with strands and mismatches is `find_longest_common_region` (SuffixTree server) / `find_longest_common_region_mismatches`.

## Core Documentation Reference

- Source: [GenomicAnalyzer.CommonRegions.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.CommonRegions.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | string[] | Yes | At least two nucleotide sequences: A/C/G/T plus IUPAC codes such as N (case-insensitive). |
| `minSupport` | integer | No | Number of sequences that must contain the substring, 2 .. N (default: all N). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Longest common strings: `{ sequence, positions[] (first occurrence per input, -1 = absent), length, support }` |
| `length` | integer | Common length of the strings (0 when there is none) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Fewer than two sequences, a sequence is null / empty or not a nucleotide sequence |
| 1002 | minSupport outside [2, N] |

## Examples

### Example 1: Rosalind LCSM sample

**Input:** `{ "sequences": ["GATTACA", "TAGACCA", "ATACA"] }`

**Response:**
```json
{ "items": [
  { "sequence": "AC", "positions": [4, 3, 2], "length": 2, "support": 3 },
  { "sequence": "CA", "positions": [5, 5, 3], "length": 2, "support": 3 },
  { "sequence": "TA", "positions": [3, 0, 1], "length": 2, "support": 3 }
], "length": 2 }
```

## See Also

- [find_common_regions](find_common_regions.md)
- [find_longest_common_region_mismatches](find_longest_common_region_mismatches.md)
