# compression_ratio

Compression-based sequence complexity (normalized Lempel-Ziv).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `compression_ratio` |
| **Method ID** | `SequenceComplexity.EstimateCompressionRatio` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Estimates sequence complexity using the normalized Lempel-Ziv complexity
`c / (n / log_b(n))`, where `c` is the number of LZ76 exhaustive-history components,
`n` the sequence length, and `b` the alphabet size (clamped to ≥ 2). Lower values
indicate more repetitive / less complex sequences. Comparison is case-insensitive.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L1772](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L1772)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `ratio` | number | Normalized Lempel-Ziv complexity |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: LZ doctest

**User Prompt:**
> Compression complexity of "1001111011000010".

**Expected Tool Call:**
```json
{
  "tool": "compression_ratio",
  "arguments": { "sequence": "1001111011000010" }
}
```

**Response:**
```json
{ "ratio": 1.5 }
```
Normalized LZ complexity of the classic LZ76 doctest string is 6/(16/log₂16) = 1.5 (LZ76 c = 6: 1/0/01/1110/1100/0010; antropy doctest).

### Example 2: Tandem repeat (ACGT × 4)

**User Prompt:**
> Compression complexity of "ACGTACGTACGTACGT".

**Expected Tool Call:**
```json
{
  "tool": "compression_ratio",
  "arguments": { "sequence": "ACGTACGTACGTACGT" }
}
```

**Response:**
```json
{ "ratio": 0.625 }
```
The repetitive ACGT×4 has low normalized LZ complexity (A/C/G/T/ACGTACGTACGT: c = 5 → 5/(16/log₄16) = 0.625).

## Performance

- **Time Complexity:** O(n²) worst case (Kaspar–Schuster scan of the LZ76 exhaustive history).

## See Also

- [dust_score](dust_score.md)
- [windowed_complexity](windowed_complexity.md)
