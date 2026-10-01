# longdust_score

Longdust complexity score S_L(x) of a sequence (Li & Li 2025).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `longdust_score` |
| **Method ID** | `SequenceComplexity.CalculateLongdustScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **longdust** score of a whole sequence x (Li & Li 2025, arXiv:2509.07357; lh3/longdust):
`S_L(x) = Σ_t log c_x(t)! − f(ℓ(x)/4^k)`, where `c_x(t)` counts k-mer t in x, `ℓ(x) = |x| − k + 1` (`kmerPositions`)
and `f(λ) = 4^k e^{−λ} Σ_n log(n!) λ^n/n!` is the expected `Σ log c!` of a random sequence — the k-mer
generalisation of the DUST score. Higher = lower complexity; longdust calls x low-complexity when
`S_L(x) − T·ℓ(x) > 0` (T = 0.6). k-mers containing a non-ACGT symbol are not counted but still count in ℓ.
`gcContent` enables longdust's GC correction (`-g`). Score 0 when ℓ ≤ 0.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L1148](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L1148)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (A/C/G/T plus IUPAC codes such as N; case-insensitive) |
| `k` | integer | No | k-mer length (longdust -k) (default 7) |
| `gcContent` | number | No | Genome GC fraction in (0, 1) for GC correction (-g); omit = uniform composition |

## Output Schema

`score` (S_L(x)), `kmerPositions` (ℓ(x) = max(0, |x| − k + 1))

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1003 | k must be 1-14; gcContent in (0, 1) |

## Examples

### Example 1: Homopolymer, k = 3

**Input:** `{"sequence": "AAAAAAAAAA", "k": 3}`

**Output:**

```json
{"score": 10.263913947178825, "kmerPositions": 8}
```

### Example 2: With GC correction

**Input:** `{"sequence": "AAAAAAAAAA", "k": 3, "gcContent": 0.41}`

**Output:**

```json
{"score": 10.230970927904716, "kmerPositions": 8}
```

### Example 3: No k-mer position

**Input:** `{"sequence": "ACGTAC"}`

**Output:**

```json
{"score": 0, "kmerPositions": 0}
```

## Performance

- O(n) plus O(ℓ) for the f table.

## See Also

- [find_longdust_regions](find_longdust_regions.md)
- [dust_score](dust_score.md)
- [find_low_complexity_intervals](find_low_complexity_intervals.md)
