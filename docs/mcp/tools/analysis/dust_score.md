# dust_score

DUST low-complexity score of a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `dust_score` |
| **Method ID** | `SequenceComplexity.CalculateDustScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the DUST low-complexity score (Morgulis et al. 2006) as
`Σ_t c_t·(c_t − 1)/2` over all overlapping words `t`, divided by `ℓ − 1`, where
`ℓ = L − wordSize + 1` is the number of words (the normaliser of NCBI dustmasker and lh3/sdust). A **higher** score indicates **lower** complexity (more repeated words);
fully distinct words give 0. Counting is case-insensitive; fewer than two words
yield 0.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L361](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L361)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `wordSize` | integer | No | Word size (default 3, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | DUST score (≥ 0); higher = lower complexity |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Word size must be at least 1 |

## Examples

### Example 1: Homopolymer (low complexity)

**Input:** `{ "sequence": "AAAAA", "wordSize": 3 }`
→ AAA occurs 3× over ℓ = 3 words → `3·2/2 / (3 − 1) = 1.5` → **Response:** `{ "score": 1.5 }`

### Example 2: All-distinct words

**Input:** `{ "sequence": "ACGT", "wordSize": 3 }`
→ ACG, CGT distinct → **Response:** `{ "score": 0.0 }`

## Performance

- **Time Complexity:** O(n). **Space Complexity:** O(distinct words).

## See Also

- [mask_low_complexity](mask_low_complexity.md)
- [find_low_complexity_regions](find_low_complexity_regions.md)
