# at_skew

Whole-sequence AT skew of a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `at_skew` |
| **Method ID** | `GcSkewCalculator.CalculateAtSkew` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the AT skew of a DNA sequence, `AT skew = (A − T) / (A + T)`, a strand-
composition asymmetry measure used alongside GC skew for replication-strand
analysis (Charneski et al. 2011; Lobry 1996). The value lies in `[−1, 1]`: `+1`
when there are no T, `−1` when there are no A. Only `A` and `T` are counted
(case-insensitive); all other symbols (including G/C) are ignored. When `A + T = 0`
the skew is defined as `0` (Biopython zero-division convention).

## Core Documentation Reference

- Source: [GcSkewCalculator.cs#L189](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L189)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `atSkew` | number | AT skew `(A − T) / (A + T)`, in [−1, 1] |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Adenine-rich

**User Prompt:**
> What is the AT skew of "AAAT"?

**Expected Tool Call:**
```json
{
  "tool": "at_skew",
  "arguments": { "sequence": "AAAT" }
}
```

**Response:**
```json
{ "atSkew": 0.5 }
```
`(A − T)/(A + T) = (3 − 1)/(3 + 1) = 0.5`.

### Example 2: G/C ignored

**User Prompt:**
> AT skew of "AAATGGGCCC"?

**Expected Tool Call:**
```json
{
  "tool": "at_skew",
  "arguments": { "sequence": "AAATGGGCCC" }
}
```

**Response:**
```json
{ "atSkew": 0.5 }
```
G and C are ignored; only A=3, T=1 count ⇒ `(3 − 1)/4 = 0.5`.

## Performance

- **Time Complexity:** O(n).
- **Space Complexity:** O(1).

## See Also

- [gc_skew](gc_skew.md) — GC skew `(G − C)/(G + C)`
- [analyze_gc_content](analyze_gc_content.md) — comprehensive GC/AT report
