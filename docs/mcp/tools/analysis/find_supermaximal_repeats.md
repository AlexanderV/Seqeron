# find_supermaximal_repeats

Find supermaximal repeats (Gusfield 1997 §7.12.1; Vmatch -supermax).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_supermaximal_repeats` |
| **Method ID** | `RepeatFinder.FindSupermaximalRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **supermaximal repeats** (Gusfield 1997 §7.12.1; Vmatch `-supermax`): maximal repeats that never occur
as a substring of another maximal repeat — on the suffix array, lcp-intervals with only singleton children whose
suffixes have pairwise distinct left characters (Abouelhoda, Kurtz & Ohlebusch 2004). Each repeat is reported once
with every (possibly overlapping) 0-based occurrence, ascending. Only A/C/G/T match (case-insensitive); a suffix at
position 0 or preceded by a non-ACGT symbol has an undefined (distinct) left character. Ordered by first
occurrence, then length.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L5469](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L5469)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `minLength` | integer | No | Minimum repeat length (default 5) |

## Output Schema

`items`: `sequence, length, positions` (all 0-based occurrences)

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | minLength must be >= 1 |

## Examples

### Example 1: Three occurrences, one overlapping

**Input:** `{"sequence": "CAGCAGCAGTTTCAGCAG", "minLength": 3}`

**Output:**

```json
{"items": [{"sequence": "CAGCAG", "length": 6, "positions": [0, 3, 12]}]}
```

### Example 2: Two supermaximal repeats (vmatch -supermax -l 3)

**Input:** `{"sequence": "ACGTACGTTTTTTTTTACGTACGT", "minLength": 3}`

**Output:**

```json
{"items": [{"sequence": "ACGTACGT", "length": 8, "positions": [0, 16]}, {"sequence": "TTTTTTTT", "length": 8, "positions": [7, 8]}]}
```

## Performance

- O(n log² n) (shared suffix-array + LCP helpers).

## See Also

- [find_direct_repeats](find_direct_repeats.md)
- [find_repeats](find_repeats.md)
- [find_approximate_direct_repeats](find_approximate_direct_repeats.md)
