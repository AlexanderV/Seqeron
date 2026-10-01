# find_low_complexity_intervals

SDUST low-complexity intervals (lh3/sdust, dustmasker -outfmt interval).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_low_complexity_intervals` |
| **Method ID** | `SequenceComplexity.FindLowComplexityIntervals` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the **SDUST** low-complexity intervals (Morgulis et al. 2006; symmetric DUST) as 0-based half-open
`[start, end)` pairs in ascending order — the native output of lh3/sdust (`name start end`) and, with `end − 1`,
of NCBI dustmasker `-outfmt interval`. These are exactly the positions `mask_low_complexity` masks. N and other
IUPAC codes split the input into independently scanned A/C/G/T runs (sdust's contract). `linker` reproduces
dustmasker `-linker` (intervals separated by fewer than `linker` unmasked bases are merged; default 1 = sdust).

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L854](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L854)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (A/C/G/T plus IUPAC codes such as N; case-insensitive) |
| `windowSize` | integer | No | SDUST window length in bases (default 64) |
| `threshold` | number | No | DUST score threshold; intervals scoring strictly above it are reported (2.0 = dustmasker level 20) (default 2.0) |
| `linker` | integer | No | dustmasker linker (default 1) |

## Output Schema

`items`: `start, end` (0-based half-open), `length`

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1003 | windowSize >= 3; threshold finite >= 0; linker 1-32 |

## Examples

### Example 1: sdust -w 64 -t 20

**Input:** `{"sequence": "ACGTGCATGCAAAAAAAAAAAAAAAAGCTAGCATCGACTGCAGCACACACACACACACAGATCGATCGTACGGTGCATGACAAAAAAAAAAAAACT"}`

**Output:**

```json
{"items": [{"start": 10, "end": 26, "length": 16}, {"start": 43, "end": 59, "length": 16}, {"start": 81, "end": 94, "length": 13}]}
```

### Example 2: dustmasker -linker 32

**Input:** `{"sequence": "ACGTGCATGCAAAAAAAAAAAAAAAAGCTAGCATCGACTGCAGCACACACACACACACAGATCGATCGTACGGTGCATGACAAAAAAAAAAAAACT", "linker": 32}`

**Output:**

```json
{"items": [{"start": 10, "end": 94, "length": 84}]}
```

### Example 3: N splits the scan

**Input:** `{"sequence": "ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA"}`

**Output:**

```json
{"items": [{"start": 6, "end": 18, "length": 12}, {"start": 23, "end": 38, "length": 15}]}
```

## Performance

- O(n · windowSize) worst case (SDUST perfect-interval list).

## See Also

- [mask_low_complexity](mask_low_complexity.md)
- [dust_score](dust_score.md)
- [find_longdust_regions](find_longdust_regions.md)
