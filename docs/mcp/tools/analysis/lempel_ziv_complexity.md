# lempel_ziv_complexity

Raw and normalized Lempel-Ziv (LZ76) complexity of a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `lempel_ziv_complexity` |
| **Method ID** | `SequenceComplexity.CalculateLempelZivComplexity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the raw **Lempel–Ziv (1976) complexity** `c` — the number of components of the exhaustive history of
the sequence (Lempel & Ziv 1976; values identical to the Kaspar–Schuster 1987 scan in antropy `_lz_complexity`,
computed via the longest-previous-factor array) — and the **normalized** value `c / (n / log_b(n))` (Zhang et al.
2009; antropy `lziv_complexity(normalize=True)`; b = number of distinct symbols, clamped to ≥ 2; equals
`compression_ratio`). Input is upper-cased; any symbol alphabet. A homopolymer of length ≥ 2 has c = 2.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L1581](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L1581)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence (any symbols; case-insensitive) |

## Output Schema

`complexity` (raw LZ76 c), `normalized` (c / (n / log_b n))

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Lempel & Ziv 1976 example 0·001·10·100·1000·101

**Input:** `{"sequence": "0001101001000101"}`

**Output:**

```json
{"complexity": 6, "normalized": 1.5}
```

### Example 2: (CAG)20 (antropy)

**Input:** `{"sequence": "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"}`

**Output:**

```json
{"complexity": 4, "normalized": 0.2484555351907228}
```

### Example 3: All distinct

**Input:** `{"sequence": "ACGT"}`

**Output:**

```json
{"complexity": 4, "normalized": 1}
```

## Performance

- O(n log² n) (suffix array + LPF).

## See Also

- [compression_ratio](compression_ratio.md)
- [windowed_complexity](windowed_complexity.md)
- [dust_score](dust_score.md)
