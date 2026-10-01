# find_clump_windows

Find (L, t)-clump k-mers together with the windows in which they form clumps.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_clump_windows` |
| **Method ID** | `KmerAnalyzer.FindClumpWindows` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Same clump definition as [find_clumps](find_clumps.md) (the Clump Finding Problem, Compeau & Pevzner
Ch. 1; Rosalind BA1E): windows are `sequence[i..i+L−1]` for `i ∈ [0, |sequence| − L]`, and an occurrence
at `p` counts when `i ≤ p ≤ i + L − k` (overlapping, case-insensitive). In addition to the k-mer set, the
tool reports **where** each k-mer forms a clump: the set of window starts in which it has at least `t`
occurrences, given as maximal runs `[firstWindowStart, lastWindowStart]` (0-based, inclusive). A run
covers `sequence[firstWindowStart .. lastWindowStart + L − 1]`. The entry's `firstWindowStart` is the
leftmost qualifying window, where `BetterClumpFinding` first detects the k-mer. Entries are ordered by
`firstWindowStart`, then ordinal k-mer; the k-mer set equals that of `find_clumps`. `windowSize` must be
at least `k`; a window longer than the sequence gives an empty list.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L2089](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L2089)
- Algorithm: [K-mer_Search.md](../../../algorithms/K-mer/K-mer_Search.md) §2.2, §7.3

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to scan (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `windowSize` | integer | Yes | Sliding window size L (≥ k) |
| `minOccurrences` | integer | Yes | Minimum occurrences t within a window (> 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `clumps` | object[] | One entry per clump k-mer, ordered by `firstWindowStart` then k-mer |
| `clumps[].kmer` | string | Clump-forming k-mer (upper case) |
| `clumps[].firstWindowStart` | integer | Leftmost window start with ≥ t occurrences |
| `clumps[].windowRuns` | object[] | Maximal runs `{ firstWindowStart, lastWindowStart }` of qualifying window starts (ascending, disjoint) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |
| 1003 | Window size must be at least k |
| 1003 | Minimum occurrences must be positive |

## Examples

### Example 1: Two clumps

**User Prompt:**
> Where do 2-mers occur at least 3 times in a 4-bp window of AAAACCCC?

**Expected Tool Call:**
```json
{
  "tool": "find_clump_windows",
  "arguments": { "sequence": "AAAACCCC", "k": 2, "windowSize": 4, "minOccurrences": 3 }
}
```

**Response:**
```json
{ "clumps": [
  { "kmer": "AA", "firstWindowStart": 0, "windowRuns": [{ "firstWindowStart": 0, "lastWindowStart": 0 }] },
  { "kmer": "CC", "firstWindowStart": 4, "windowRuns": [{ "firstWindowStart": 4, "lastWindowStart": 4 }] }
] }
```

### Example 2: Split runs (Rosalind BA1B sample, k=4, L=12, t=2)

**Expected Tool Call:**
```json
{
  "tool": "find_clump_windows",
  "arguments": { "sequence": "ACGTTGCATGTCGCATGATGCATGAGAGCT", "k": 4, "windowSize": 12, "minOccurrences": 2 }
}
```

**Response:**
```json
{ "clumps": [
  { "kmer": "GCAT", "firstWindowStart": 4, "windowRuns": [{ "firstWindowStart": 4, "lastWindowStart": 5 }, { "firstWindowStart": 11, "lastWindowStart": 12 }] },
  { "kmer": "CATG", "firstWindowStart": 5, "windowRuns": [{ "firstWindowStart": 5, "lastWindowStart": 6 }, { "firstWindowStart": 12, "lastWindowStart": 13 }] },
  { "kmer": "ATGA", "firstWindowStart": 13, "windowRuns": [{ "firstWindowStart": 13, "lastWindowStart": 14 }] }
] }
```
Expected values: Python brute force over every window (K-mer_Search.md §7.3), equal to the library output.

## Performance

- **Time Complexity:** O(n · k) (one sliding pass shared with `find_clumps`). **Space Complexity:** O(L) window state plus the output.

## See Also

- [find_clumps](find_clumps.md) — the clump k-mer set only
- [kmer_positions](kmer_positions.md) — all occurrence positions of one k-mer
