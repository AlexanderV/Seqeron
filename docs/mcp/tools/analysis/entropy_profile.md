# entropy_profile

Shannon entropy in sliding windows along a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `entropy_profile` |
| **Method ID** | `SequenceStatistics.CalculateEntropyProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Produces the Shannon-entropy profile of a sequence: for each sliding window of length
`windowSize` (advanced by `stepSize`) it reports `H = −Σ pᵢ log₂ pᵢ` in bits over the window's
letter frequencies (case-insensitive; every letter is its own symbol — `N` and IUPAC codes are distinct
symbols and `T`/`U` are not merged; non-letters such as gaps, digits and `*` are ignored). Windows are emitted for offsets
`0, stepSize, 2·stepSize, …` up to `length − windowSize`. Maximum entropy is `log₂ m` for m
distinct symbols (2 bits for the 4-letter DNA alphabet).

- **`k`** (optional, ≥ 1): k-mer (block) entropy instead — every overlapping k-substring of the
  upper-cased window is a symbol (canonical `SequenceComplexity.CalculateKmerEntropy`); a window shorter
  than k scores 0.
- **`includePartialWindow`**: after the complete windows, every further start `k·stepSize < length` is
  scored on the truncated trailing slice (`range(0, len(seq), step)` slicing); a slice without letters scores 0.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L2577](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L2577) (letters), [#L2319](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L2611) (partial windows), [#L2361](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L2653) (k-mers)
- Algorithm: `docs/algorithms/Statistics/Entropy_Profile.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `windowSize` | integer | No | Window size (default 50, ≥ 1) |
| `stepSize` | integer | No | Step size (default 1, ≥ 1) |
| `k` | integer | No | k-mer length (≥ 1) for block entropy; omit for letter entropy |
| `includePartialWindow` | boolean | No | `true`: also score trailing windows shorter than `windowSize` (default `false`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `values` | number[] | Per-window entropy in bits |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Window size must be at least 1 |
| 1003 | Step size must be at least 1 |
| 1003 | k must be at least 1 |

## Examples

### Example 1: Homopolymer (zero entropy)

**Input:** `{ "sequence": "AAAA", "windowSize": 2, "stepSize": 1 }`
→ three "AA" windows → **Response:** `{ "values": [0, 0, 0] }`

### Example 2: Uniform window (max entropy)

**Input:** `{ "sequence": "ATGC", "windowSize": 4 }`
→ one window with equal base frequencies → **Response:** `{ "values": [2.0] }`

### Example 3: k-mer entropy with partial windows

**Input:** `{ "sequence": "ACGTACGTAAAACCCGGT", "windowSize": 8, "stepSize": 3, "k": 2, "includePartialWindow": true }`

**Response:** `{ "values": [1.950212064914747, 2.2359263506290326, 2.128085278891394, 2.2359263506290326, 1.9219280948873625, 1.0] }`
(`scipy.stats.entropy` (base 2) of the overlapping dinucleotide counts of `s[i:i+8]`, `i in range(0, 18, 3)`.)

## Performance

- **Time Complexity:** O(n · windowSize). **Space Complexity:** O(n / stepSize).

## See Also

- [windowed_complexity](windowed_complexity.md)
- [gc_content_profile](gc_content_profile.md)
